# An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics

**Liu J, Lichtenberg T, Hoadley KA, et al. Cell. 2018 Apr 5;173(2):400-416.e11**

- DOI: 10.1016/j.cell.2018.02.052
- Cell URL: <https://www.cell.com/cell/fulltext/S0092-8674(18)30229-0>
- PubMed: <https://pubmed.ncbi.nlm.nih.gov/29625055/>

## Why this paper is in our reference list

It's the canonical curated source of TCGA clinical outcome endpoints — the
authors went back through TCGA's heterogeneous clinical data and produced
recommended definitions for the four major survival endpoints across all 33
TCGA tumor types:

- **OS** (overall survival)
- **DSS** (disease-specific survival)
- **DFI** (disease-free interval)
- **PFI** (progression-free interval)

They also flag which endpoints they consider *usable* for each tumor type
(e.g. PFI is more reliable than OS for low-mortality cancers like THCA).

## How we'll use it

When we add survival-analysis pipeline examples to this project (downstream
of the per-patient parquet rows), we should:

1. Compare survival curves derived from our parquet rows against the
   recommended endpoint definitions in this paper.
2. Use the per-tumor-type usability table (Table 3 / supplementary Table S1)
   when picking which TCGA projects to demo on for each survival endpoint —
   e.g. don't show OS curves for cancers where the paper says event counts
   are too low.
3. Cross-check that the days_to_* fields in our patient rows produce the
   same OS/DSS/DFI/PFI numbers the paper reports per project.

## Relevant supplementary data

The paper ships a curated CDR (Clinical Data Resource) table — it's hosted
on the [PanCanAtlas Publications page](https://gdc.cancer.gov/about-data/publications/pancanatlas)
as `TCGA-CDR-SupplementalTableS1.xlsx`. That's the operational reference
we'd diff against; the PDF is for the methodology rationale.
