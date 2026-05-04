# tcga2hf

Download public, **open-access** TCGA data from the NCI Genomic Data Commons
(GDC) and stage it as a HuggingFace Hub dataset.

**Status:** prototype. Schema, CLI surface, and dataset layout are still
evolving — the bits documented below are the stable ones.

## Install

```bash
uv sync
```

## Auth

Only `tcga2hf upload` needs auth. Put a write-scoped HF token (from
https://huggingface.co/settings/tokens) into a local `.env` file in the project
root — gitignored, auto-loaded by the CLI:

```
HF_TOKEN=hf_xxxxxxxxxxxxxxxxxxxx
```

The CLI calls `load_dotenv(override=True)` so this project's `.env` always
wins over any inherited shell `HF_TOKEN`.

## Data location

All data lives outside the repo, under `$TCGA2HF_DATA_DIR` (default
`$HOME/data/tcga2hf`), overridable per command via `--data-dir`.

## Commands

```
tcga2hf fetch-clinical     # GDC -> raw cases.json per project
tcga2hf fetch-mutations    # GDC -> raw MAF files per project (DNA, WXS)
tcga2hf fetch-expression   # GDC -> raw STAR gene-count TSVs per project (RNA-Seq)
tcga2hf build              # raw -> per-project Parquet + dataset card
tcga2hf upload             # push processed/ to a HF dataset repo
```

`tcga2hf <cmd> --help` for arguments.

## GDC data model primer

Quick orientation for anyone extending the schema. The GDC organizes
everything around a **case** (one patient) with a hierarchy of physical
biospecimens beneath it, plus assay-derived data files keyed to those
biospecimens.

```
case          one patient (TCGA-XX-1234)
└── sample    physical specimen taken at one timepoint
              (Primary Tumor, Solid Tissue Normal, Blood Derived Normal, ...)
    └── portion    a piece of that sample for a specific lab process
        └── analyte    extracted material of one type (DNA or RNA)
            └── aliquot    a vial of that analyte handed off for sequencing
```

Our patient-row schema **preserves this hierarchy verbatim** — no flattening,
no field hoisting. Portion/analyte-level QC fields like `is_ffpe`,
`a260_a280_ratio`, and `normal_tumor_genotype_snp_match` survive intact.
Convenience flat-views are provided by `TcgaHfPatient` in `tcga2hf.models`
(see below).

**Tumor vs normal.** Somatic variant calling pairs a tumor sample with a
matched normal (blood or adjacent tissue) from the same patient to identify
acquired mutations. Each MAF variant carries both aliquot UUIDs and we
resolve them back to `tumor_sample_id` / `matched_normal_sample_id`.

**GDC file taxonomy.** Each fetched data file has three GDC labels:
- `data_category` — broad bucket (Simple Nucleotide Variation, Transcriptome Profiling)
- `data_type` — specific data product (Masked Somatic Mutation, Gene Expression Quantification)
- `experimental_strategy` — underlying assay (WXS, RNA-Seq)

We name top-level molecular columns `samples_<data_type_snake_case>` so the
source product is unambiguous (e.g. `samples_masked_somatic_mutation`).

Authoritative references:
[GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/) ·
[Biospecimen Encyclopedia](https://docs.gdc.cancer.gov/Encyclopedia/pages/Biospecimen/) ·
[MAF format spec](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) ·
[Expression pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/) ·
[Sample type codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) ·
[TCGA barcode reference](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/).

## Pydantic reference implementation

`src/tcga2hf/models.py` ships a fully-typed `TcgaHfPatient` pydantic model that
mirrors the parquet schema and adds convenience joins for the common
operations: `aliquot_to_sample`, `tumor_normal_pairs`, `mutations_by_gene`,
`expression_for_gene`, `timeline`, etc. Slow but explicit — intended as a
reference for what the data means and a spec other implementations
(polars / pyarrow / SQL) can be benchmarked against.

```python
import pyarrow.parquet as pq
from tcga2hf.models import TcgaHfPatient

t = pq.read_table("~/data/tcga2hf/processed/TCGA-CHOL/train.parquet")
for row in t.to_pylist():
    patient = TcgaHfPatient.model_validate(row)
    for tumor, normal in patient.tumor_normal_pairs():
        print(tumor.submitter_id, "vs", normal.submitter_id)
    for variant in patient.mutations_by_gene().get("TP53", []):
        print(variant.HGVSp_Short, variant.t_alt_count, "/", variant.t_depth)
    for event in patient.timeline():
        print(f"day {event.day:+.0f}: {event.category} — {event.label}")
```

## Tests

```bash
uv run pytest                   # all, incl. live API smoke
uv run pytest -m "not network"  # offline only
```

## Data redistribution

Only the GDC's **open-access** tier is fetched (`/cases` with `access=open`).

Per the [NCI GDC Data Analysis Policy][gdc-policy]:

> The GDC itself places no restrictions (other than attempts at reidentification)
> on analysis or publication of open access data provided through the GDC Data Portal.

Per the [NCI TCGA citation page][nci-cite]:

> Moratoria on all cancer types are now lifted and all TCGA data are available
> without restrictions on their use in publications or presentations.

Per the [GDC Data Access Processes and Tools page][gdc-access]:

> Open access data generally includes high level genomic data that is not
> individually identifiable, as well as most clinical and all biospecimen data
> elements.

Two obligations carry over to anyone using or redistributing this data:

1. **No re-identification.** From the GDC policy:
   > Users of any data provided by GDC, whether open or controlled access, agree
   > not to attempt to reidentify any individual participant in any study
   > represented by GDC data, for any purpose whatever.

2. **Acknowledgement on publication** (per NCI):
   > The results <published or shown> here are in whole or part based upon data
   > generated by the TCGA Research Network: https://www.cancer.gov/tcga.

The dataset card emitted by `tcga2hf build` carries the same quotes and links
so they travel with the dataset onto the HF Hub.

[gdc-policy]: https://gdc.cancer.gov/analyze-data/data-analysis-policies
[nci-cite]: https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing
[gdc-access]: https://gdc.cancer.gov/access-data/data-access-processes-and-tools

References:
[GDC Policies](https://gdc.cancer.gov/about-gdc/gdc-policies),
[GDC Encyclopedia — Controlled Access](https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access/) (defines what is *not* fetched),
[NIH Genomic Data Sharing Policy](https://sharing.nih.gov/genomic-data-sharing).
