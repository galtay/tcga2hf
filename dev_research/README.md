# dev_research

Notes on external papers and resources we want to consult or cross-check
against as the pipeline grows. Not user-facing — meant for developers of this
project to track relevant prior art.

## Layout

One directory per paper, named `<first_author>_<year>/`. Each holds the
original PDF, our markdown notes (`notes.md`), and any analysis
notebooks or scratch files we generate while reproducing the paper's
results.

## Contents

- [`liu_2018/`](liu_2018/) — Liu et al., Cell 2018. Curated TCGA clinical
  data resource for survival outcome analytics. First reproduce-and-validate
  target — now wired into the pipeline as `cdr_*` curated columns +
  `survival.py` re-derivations.
- [`hoadley_2018/`](hoadley_2018/) — Hoadley et al., Cell 2018. 28
  integrative clusters across 10,000 tumors / 33 cancer types via joint
  mRNA / miRNA / methylation / CNV / mutation / RPPA clustering. Next
  reproduce-and-validate target.
