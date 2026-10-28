# Legacy analysis (October 2024)

The original version of this project, preserved unchanged for reference. It is
superseded by the pipeline in `src/flavi_qsar/` and is not part of the current
analysis. Nothing here should be imported or reused.

| File | What it was |
| --- | --- |
| `DataGen.ipynb` | ChEMBL pull for target CHEMBL5980 (DENV2 NS3), IC50 filtering, activity binning |
| `Exploratory_data_analysis.ipynb` | RDKit Lipinski descriptors, pIC50 conversion, plots, Mann-Whitney test |
| `Bioactivity.csv` | Raw pull, 1,284 IC50 rows (ChEMBL release not recorded) |
| `Bioactivity_preprocessed.csv` | Derived dataset. **Contains incorrect labels, see below.** |
| `mannwhitneyu_pIC50.csv` | Saved test result. **Invalid, see below.** |

## Why it was replaced

**The activity labels are misaligned.** In `DataGen.ipynb`, `bioactivity_class`
was built as a fresh `RangeIndex` Series and then joined with
`pd.concat(axis=1)` onto a DataFrame that had kept its original gapped index
after 47 rows were dropped. Pandas aligned on index, so from the first dropped
row onward every label is attached to the wrong molecule. In
`Bioactivity_preprocessed.csv`, 104 of 1,192 rows contradict the code's own
thresholds: 28 of the 31 rows labelled "active" have IC50 above the 1,000 nM
cutoff, with a median of 6.7 uM.

**The Mann-Whitney test was circular.** Groups were defined by thresholding
pIC50, then compared on pIC50, so a small p-value was guaranteed by
construction regardless of the data. The saved cell also cannot run as written
(it returns an undefined `results`), so the CSV came from an earlier revision.

**The pooled data was confounded.** 1,013 of the 1,284 rows came from a single
PubChem HTS screen (CHEMBL1794550) whose lowest value is 1,610 nM and whose
median is 46 uM: it contains no actives at all. Every active in the dataset came
from the medicinal-chemistry papers, so "active" was perfectly correlated with
data source rather than with protease inhibition.

The current pipeline keeps only biochemical (assay type B) NS2B-NS3 protease
measurements, deduplicates on standardized InChIKey, treats censored values
separately, and models pIC50 as a continuous target.
