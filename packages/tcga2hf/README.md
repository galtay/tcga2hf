# tcga2hf

Typed pydantic models + pyarrow schemas for the public
[`gabrielaltay/tcga-patients-open`][patients] and
[`gabrielaltay/tcga-tabular-open`][tabular] HuggingFace datasets.

This is the **read-side** companion package — install it to validate
patient rows or browse them with full type information. To produce or
publish the datasets yourself, see
[`tcga2hf-pipeline`](../tcga2hf-pipeline/).

## Install

```bash
pip install tcga2hf
```

Dependencies are limited to `pydantic` and `pyarrow`.

## Usage

```python
import pyarrow.parquet as pq
from tcga2hf import TcgaHfPatient

# Read one project's parquet from the consolidated dataset.
t = pq.read_table("TCGA-LUAD/data.parquet")
for row in t.to_pylist():
    patient = TcgaHfPatient.model_validate(row)
    for tumor, normal in patient.tumor_normal_pairs():
        print(tumor.submitter_id, "vs", normal.submitter_id)
    for variant in patient.mutations_by_gene().get("TP53", []):
        print(variant.HGVSp_Short, variant.t_alt_count, "/", variant.t_depth)
    for event in patient.timeline():
        print(f"day {event.day:+.0f}: {event.category} — {event.label}")
```

`TcgaHfPatient` mirrors the consolidated parquet schema field-for-field
(`extra="forbid"` strict mode) and adds convenience joins:

- `aliquot_to_sample` / `aliquot_lookup` — flatten the GDC biospecimen tree
- `tumor_normal_pairs` — distinct tumor/normal sample pairs from MAF rows
- `mutations_by_gene`, `mutations_by_consequence`
- `expression_for_gene` — per-aliquot expression lookup by gene symbol
- `timeline` — every dated event (consent, diagnosis, treatments,
  follow-ups, sample procurement, BCR receipt, death) sorted on the
  case's `index_date` anchor
- `consistency_check` — count of GDC quirks worth surfacing for
  downstream consumers (e.g. samples with `days_to_collection >
  days_to_death`)

The pyarrow `*_FIELDS` lists in `tcga2hf.schema` are the single source of
truth for the dataset shape (regenerated from gdcdictionary YAMLs); the
pydantic models are derived from those same lists, so the two stay in
sync by construction.

## Source

This package ships from the [`galtay/tcga2hf`][repo] monorepo as one of
two workspace members; see the repo for the build pipeline, dataset
cards, and full documentation.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
[repo]: https://github.com/galtay/tcga2hf
