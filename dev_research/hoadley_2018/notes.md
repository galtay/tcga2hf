# Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer

**Hoadley KA, Yau C, Hinoue T, ... Stuart JM, Benz CC, Laird PW. Cell. 2018 Apr 5;173(2):291-304**

- DOI: 10.1016/j.cell.2018.03.022
- Cell URL: <https://www.sciencedirect.com/science/article/pii/S0092867418303027>
- Local PDF: [`1-s2.0-S0092867418303027-main.pdf`](1-s2.0-S0092867418303027-main.pdf)

## Why this paper is in our reference list

Companion to the Liu 2018 CDR paper — same TCGA Pan-Cancer Atlas issue, same
~10,000-tumor cohort. Where Liu curated the *clinical* outcome variables,
Hoadley produced the canonical *molecular* classification: 28 integrative
clusters ("iClusters") derived by jointly clustering across six data
modalities — mRNA-Seq, miRNA-Seq, DNA methylation, copy number, somatic
mutation, RPPA protein. Their grouping defines pan-cancer organ-system
cohorts (Pan-Squamous, Pan-GYN, Pan-Kidney, Pan-GI) and serves as the
reference label set for downstream multi-modal benchmarks.

## How we'll use it

This is the next paper in the *reproduce-and-validate* program after Liu
2018 CDR. The plan parallels what `survival.py` did for clinical
endpoints:

1. Pull each patient's published iCluster assignment (supplementary
   Table S2 of the paper) and attach as a `hoadley_iCluster` column on
   the patient rows — same pattern as `cdr_*` (curated stream alongside
   raw data).
2. Re-derive the iCluster assignment from current GDC data using the
   paper's documented clustering pipeline, where feasible. Validate the
   re-derivation against the published assignments.
3. The re-derivation requires modalities we don't yet ship: DNA
   methylation, copy number, miRNA-Seq, RPPA. Adding those is a
   prerequisite — currently we have somatic mutations + gene expression
   only. That's a separate workstream, but the validation rig is the
   same shape as `scripts/validate_survival.py`.

## Relevant supplementary data

- Table S2: iCluster assignments per sample (the validation target).
- Tables S3–S6: per-modality cluster assignments. Useful for partial
  validation when only some modalities are available.
- The 28 iClusters are not labels we should generate from scratch — even
  with all modalities present, the integrative clustering depends on
  joint factor decomposition that can drift across runs. The audit pattern
  is "same data, same pipeline, same labels" rather than re-clustering
  from scratch.

## Related papers

- Liu et al. 2018 CDR ([note](../liu_2018/notes.md)) — clinical endpoints;
  first reproduction target, now wired into the pipeline as `cdr_*` curated
  columns + `survival.py` re-derivations.
