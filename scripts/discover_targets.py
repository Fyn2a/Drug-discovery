"""Phase 1: find and size candidate flaviviral protease data in ChEMBL.

Produces two review tables and retrieves no data beyond the metadata needed to
build them, so the target and assay selection can be confirmed before the
archived bulk pull in Phase 2.

Targets are found by filtering the `target` endpoint on organism rather than by
free-text search. ChEMBL's relevance search is not stable across releases: as of
ChEMBL_37 a search for "dengue" returns only ORGANISM-level entries and does
not surface CHEMBL5980, the single-protein target this project is built on.

Flaviviral targets are also mostly registered as "Genome polyprotein", so one
target mixes protease, helicase and NS5 measurements. Enzyme assignment
therefore has to come from the assay description, not from the target name, and
that is what `assay_candidates.csv` is for.

Usage:
    python scripts/discover_targets.py
"""

from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from flavi_qsar.chembl import ChemblClient  # noqa: E402

# Organism substrings matched case-insensitively against the target endpoint.
ORGANISM_QUERIES = (
    "dengue",
    "Zika",
    "West Nile",
    "encephalitis",
    "yellow fever",
    "Flavivirus",
    "Usutu",
    "Powassan",
)

OUTPUT_DIR = PROJECT_ROOT / "data" / "interim"

# Assay ceiling values are the top of the dose series rather than measurements.
# Anything at or above this is treated as non-informative for a regression.
CEILING_NM = 99_000.0

# Checked in order: the NS5 and helicase patterns are tested before the protease
# pattern because "NS2B-NS3 protease" and "NS3 helicase" both contain "NS3".
ENZYME_PATTERNS = (
    ("NS5/capping", re.compile(
        r"NS5|methyltransfer|polymerase|RdRp|capping|guanylyltransfer|SAM site",
        re.IGNORECASE)),
    ("helicase", re.compile(r"helicase|NTPase|unwinding|ATPase", re.IGNORECASE)),
    ("protease", re.compile(r"NS2B|NS3\s*pro|protease|peptidase", re.IGNORECASE)),
)


def classify_enzyme(description: str) -> str:
    """Assign an assay to an enzyme activity from its description."""
    for label, pattern in ENZYME_PATTERNS:
        if pattern.search(description):
            return label
    return "unclassified"


def summarise_assay(records: list[dict[str, Any]]) -> dict[str, Any]:
    """Size and characterise one assay's IC50 records."""
    first = records[0]
    description = first.get("assay_description") or ""

    values = [
        float(record["standard_value"])
        for record in records
        if record.get("standard_value") and record.get("standard_units") == "nM"
    ]
    exact = [
        record
        for record in records
        if record.get("standard_relation") in ("=", None, "")
        and record.get("standard_units") == "nM"
        and record.get("standard_value")
        and float(record["standard_value"]) < CEILING_NM
    ]
    censored = [
        record for record in records if record.get("standard_relation") in (">", ">=")
    ]
    at_ceiling = [value for value in values if value >= CEILING_NM]

    return {
        "assay_chembl_id": first["assay_chembl_id"],
        "target_chembl_id": first["target_chembl_id"],
        "document_chembl_id": first.get("document_chembl_id") or "",
        "assay_type": first.get("assay_type") or "",
        "enzyme": classify_enzyme(description),
        "n_records": len(records),
        "n_usable_exact": len(exact),
        "n_censored": len(censored),
        "n_at_ceiling": len(at_ceiling),
        "n_molecules": len({record.get("molecule_chembl_id") for record in records}),
        "median_nm": round(statistics.median(values), 1) if values else "",
        "min_nm": round(min(values), 1) if values else "",
        "description": description,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    client = ChemblClient()
    release = client.release()
    print(f"ChEMBL release: {release}\n")

    targets: dict[str, dict[str, Any]] = {}
    for query in ORGANISM_QUERIES:
        found = list(client.iter_records("target", page_size=200, organism__icontains=query))
        single = [t for t in found if t.get("target_type") == "SINGLE PROTEIN"]
        for target in single:
            targets.setdefault(target["target_chembl_id"], target)
        print(f"  organism~{query:24s} {len(found):3d} targets, {len(single):3d} single-protein")

    print(f"\n{len(targets)} unique single-protein targets. Sizing IC50 data.\n")

    target_rows: list[dict[str, Any]] = []
    with_data: list[str] = []
    for chembl_id, target in sorted(targets.items()):
        n_ic50 = client.count("activity", target_chembl_id=chembl_id, standard_type="IC50")
        if n_ic50 == 0:
            continue
        with_data.append(chembl_id)
        target_rows.append(
            {
                "target_chembl_id": chembl_id,
                "pref_name": target.get("pref_name") or "",
                "organism": target.get("organism") or "",
                "tax_id": target.get("tax_id") or "",
                "n_ic50": n_ic50,
                "n_ic50_biochemical": client.count(
                    "activity",
                    target_chembl_id=chembl_id,
                    standard_type="IC50",
                    assay_type="B",
                ),
            }
        )

    target_rows.sort(key=lambda row: row["n_ic50"], reverse=True)

    print("Retrieving assay-level metadata.\n")
    by_assay: dict[str, list[dict[str, Any]]] = defaultdict(list)
    cache: list[dict[str, Any]] = []
    for chembl_id in with_data:
        records = list(
            client.iter_activities(
                page_size=1000, target_chembl_id=chembl_id, standard_type="IC50"
            )
        )
        cache.extend(records)
        for record in records:
            by_assay[record["assay_chembl_id"]].append(record)

    assay_rows = [summarise_assay(records) for records in by_assay.values()]
    assay_rows.sort(key=lambda row: row["n_usable_exact"], reverse=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUTPUT_DIR / "target_candidates.csv", target_rows)
    write_csv(OUTPUT_DIR / "assay_candidates.csv", assay_rows)
    (OUTPUT_DIR / "activities_ic50_cache.json").write_text(
        json.dumps(cache), encoding="utf-8"
    )
    (OUTPUT_DIR / "target_discovery_requests.json").write_text(
        json.dumps({"chembl_release": release, "requests": client.request_log}, indent=2),
        encoding="utf-8",
    )

    print(f"{'target':14s} {'IC50':>6s} {'bioch':>6s}  name | organism")
    print("-" * 96)
    for row in target_rows:
        print(
            f"{row['target_chembl_id']:14s} {row['n_ic50']:6d} {row['n_ic50_biochemical']:6d}  "
            f"{row['pref_name']} | {row['organism']}"
        )

    print(f"\nUsable exact IC50 by enzyme and assay type:")
    totals: dict[tuple[str, str], int] = defaultdict(int)
    for row in assay_rows:
        totals[(row["enzyme"], row["assay_type"])] += row["n_usable_exact"]
    for (enzyme, assay_type), total in sorted(totals.items(), key=lambda kv: -kv[1]):
        print(f"  {enzyme:14s} type {assay_type:2s}  {total:6d}")

    print(f"\nWrote target_candidates.csv  ({len(target_rows)} targets)")
    print(f"Wrote assay_candidates.csv   ({len(assay_rows)} assays)")
    print(f"Cached {len(cache)} IC50 records for Phase 2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
