"""Minimal ChEMBL REST client.

Standard library only, so data retrieval carries no third-party dependency and
every request URL can be recorded alongside the data it produced. The official
`chembl_webresource_client` is deliberately not used: it hides the underlying
URLs and its implicit paging makes a pull hard to reproduce exactly. Here the
full request log is available on the client and gets archived with each pull.
"""

from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from typing import Any, Iterator

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

# Retried on the assumption that they are transient. Anything else (400, 404)
# indicates a malformed query and is raised immediately.
RETRY_STATUS = frozenset({429, 500, 502, 503, 504})


class ChemblError(RuntimeError):
    """A ChEMBL request failed in a way that retrying will not fix."""


@dataclass
class ChemblClient:
    """Paging HTTP client for the ChEMBL data API.

    Attributes:
        timeout: Per-request socket timeout in seconds.
        max_retries: Attempts per request before giving up.
        backoff: Base for exponential retry delay, in seconds.
        pause: Delay between successful requests, to stay polite to the API.
        request_log: Every URL requested, in order, for provenance.
    """

    base_url: str = BASE_URL
    timeout: float = 60.0
    max_retries: int = 4
    backoff: float = 2.0
    pause: float = 0.2
    request_log: list[str] = field(default_factory=list)

    def _request(self, url: str) -> dict[str, Any]:
        """GET a JSON document, retrying transient failures."""
        last_error: Exception | None = None

        for attempt in range(self.max_retries):
            if attempt:
                time.sleep(self.backoff**attempt)
            try:
                request = urllib.request.Request(
                    url, headers={"Accept": "application/json"}
                )
                with urllib.request.urlopen(request, timeout=self.timeout) as response:
                    payload: dict[str, Any] = json.load(response)
                self.request_log.append(url)
                time.sleep(self.pause)
                return payload
            except urllib.error.HTTPError as error:
                if error.code not in RETRY_STATUS:
                    raise ChemblError(f"HTTP {error.code} for {url}") from error
                last_error = error
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
                last_error = error

        raise ChemblError(
            f"{self.max_retries} attempts failed for {url}"
        ) from last_error

    def _get(self, path: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
        """GET an endpoint by path, which may be relative or server-absolute.

        `page_meta.next` values returned by ChEMBL are absolute paths rather
        than full URLs, so both forms are accepted here.
        """
        if path.startswith("http"):
            url = path
        elif path.startswith("/"):
            scheme, netloc = urllib.parse.urlsplit(self.base_url)[:2]
            url = urllib.parse.urlunsplit((scheme, netloc, path, "", ""))
        else:
            url = f"{self.base_url}/{path.lstrip('/')}"

        if params:
            separator = "&" if "?" in url else "?"
            url = f"{url}{separator}{urllib.parse.urlencode(params)}"

        return self._request(url)

    def status(self) -> dict[str, Any]:
        """Service status, including the current database release."""
        return self._get("status.json")

    def release(self) -> str:
        """Current ChEMBL release identifier, for example `ChEMBL_37`."""
        version = self.status().get("chembl_db_version")
        if not version:
            raise ChemblError("status endpoint did not report chembl_db_version")
        return str(version)

    def search_targets(self, query: str, limit: int = 100) -> list[dict[str, Any]]:
        """Free-text target search, ordered by the relevance score ChEMBL assigns."""
        payload = self._get("target/search.json", {"q": query, "limit": limit})
        return list(payload.get("targets", []))

    def count(self, endpoint: str, **filters: Any) -> int:
        """Total records matching `filters`, without retrieving them.

        Asks for a single record and reads `page_meta.total_count`, which is the
        cheapest way to size a query before committing to a full pull.
        """
        payload = self._get(f"{endpoint}.json", {**filters, "limit": 1})
        return int(payload["page_meta"]["total_count"])

    @staticmethod
    def _records(page: dict[str, Any]) -> list[dict[str, Any]]:
        """Pull the record list out of a page.

        ChEMBL keys the list by a pluralised form of the endpoint (`activity`
        returns `activities`, `target` returns `targets`), so the key is found
        by elimination rather than derived from the endpoint name.
        """
        keys = [key for key in page if key != "page_meta"]
        if len(keys) != 1:
            raise ChemblError(f"cannot identify record list in page keys: {keys}")
        return list(page[keys[0]])

    def iter_records(
        self, endpoint: str, page_size: int = 1000, **filters: Any
    ) -> Iterator[dict[str, Any]]:
        """Yield every record matching `filters`, following paging to the end."""
        page = self._get(f"{endpoint}.json", {**filters, "limit": page_size})

        while True:
            yield from self._records(page)
            next_path = page.get("page_meta", {}).get("next")
            if not next_path:
                return
            page = self._get(next_path)

    def iter_activities(
        self, page_size: int = 1000, **filters: Any
    ) -> Iterator[dict[str, Any]]:
        """Yield activity records matching `filters`."""
        return self.iter_records("activity", page_size=page_size, **filters)
