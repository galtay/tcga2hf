from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from tcga2hf.genomic import MODALITY_FILTERS as _MODALITY_FILTERS

# ---------------------------------------------------------------------------
# Static body sections. Single source of truth for content that's also
# referenced from the repo README — the README links here rather than
# duplicating the text.
# ---------------------------------------------------------------------------


_DATA_MODEL = """\
## Data model

Closely follows the GDC data model — see the [GDC Data Dictionary][gdc-dict]
for the canonical entity-by-entity definitions. Each row is one `case` (one
patient) with the full biospecimen subtree:

```
case          one patient (TCGA-XX-1234)
└── sample    physical specimen taken from the patient at one timepoint
              (Primary Tumor, Solid Tissue Normal, Blood Derived Normal, ...)
    └── portion    a piece of that sample for a specific lab process
        └── analyte    extracted material of one type (DNA or RNA)
            └── aliquot    a vial of that analyte handed off for sequencing
```

Every `days_to_*` field anchors to the case's `index_date` (TCGA: almost
always `"Diagnosis"`), per the dictionary, so clinical and biospecimen events
share a single timeline.

### Where this dataset deviates from the GDC

The few places this row layout differs from a direct mapping of the GDC `case`
tree:

- **Top-level convenience columns.** Each row carries `gdc_portal_url`
  (templated link to the patient's GDC Data Portal page) and
  `samples_<gdc_data_type_snake_case>` molecular vectors (e.g.
  `samples_masked_somatic_mutation`,
  `samples_gene_expression_quantification`). These let consumers
  column-project just the modalities they need; each entry carries foreign
  keys (FKs) back to `samples[].portions[].analytes[].aliquots[]`.
- **Resolved sample FKs on Mutation Annotation Format (MAF) rows.** The GDC
  ships MAF variants with aliquot UUIDs in `Tumor_Sample_UUID` /
  `Matched_Norm_Sample_UUID`; we additionally resolve those to
  `tumor_sample_id` / `matched_normal_sample_id` so consumers can join
  straight to `samples[]`.
- **Lifted expression QC counts.** Each Gene Expression Quantification
  record has the STAR per-feature quality-control counts `N_unmapped`,
  `N_multimapping`, `N_noFeature`, `N_ambiguous` lifted from the source
  Tab-Separated Values (TSV) file onto the row as scalar fields. The
  `stranded_first` / `stranded_second` columns are dropped — the GDC
  pipeline harmonizes by [treating all RNA-Seq reads as
  unstranded](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/),
  so `unstranded` is the canonical column.
"""


_LOADING = """\
## Loading

The [`tcga2hf` package][repo] ships a typed `TcgaHfPatient` pydantic model
that mirrors this schema and adds convenience joins (tumor/normal pairs,
mutations-by-gene, expression-by-gene, longitudinal timeline).

```python
import pyarrow.parquet as pq
from tcga2hf.models import TcgaHfPatient

t = pq.read_table("TCGA-CHOL/data.parquet")
patients = [TcgaHfPatient.model_validate(r) for r in t.to_pylist()]
```
"""


def _render_gdc_requests(projects: list[str]) -> str:
    """Render a concise summary of the GDC requests this build issues.

    Templated from the live constants in `tcga2hf.clinical` and
    `tcga2hf.genomic` so the docs can't drift from the code. Full
    request payloads (the JSON `filters` + `fields` + `expand` blocks)
    live in those source files; this section only shows the headline
    filters that distinguish each table.
    """
    project_md = ", ".join(f"`{p}`" for p in sorted(projects))

    # Per-modality clause table. One row per (table, modality filter set);
    # populated from `genomic.MODALITY_FILTERS` so adding a modality
    # automatically documents itself here.
    modality_rows: list[str] = []
    table_for_data_type = {
        "Masked Somatic Mutation": "masked_somatic_mutation",
        "Gene Expression Quantification": "gene_expression_quantification",
    }
    for dtype, extras in _MODALITY_FILTERS.items():
        table = table_for_data_type.get(dtype, dtype)
        cols = [f"`{table}`", f"`{dtype}`"]
        for col_field in (
            "data_format",
            "data_category",
            "experimental_strategy",
            "analysis.workflow_type",
        ):
            cols.append(f"`{extras.get(col_field, '')}`")
        modality_rows.append("| " + " | ".join(cols) + " |")
    modality_table = "\n".join(modality_rows)

    return f"""\
## How this dataset is built

Three GDC REST endpoints feed every row, with filters and field lists
constructed in [`src/tcga2hf/clinical.py`][src-clinical] (the `cases`
table) and [`src/tcga2hf/genomic.py`][src-genomic] (molecular and
provenance tables) of the [`tcga2hf` repo][repo].

**Projects fetched in this build:** {project_md}.

### `POST /cases` → `cases` table

Filter: `project.project_id IN [<projects above>]`. The request `expand`s
the full nested `case` structure (demographic + diagnoses → treatments +
follow_ups + exposures + family_histories + samples → portions →
analytes → aliquots), and the response is captured verbatim into
`<data-dir>/raw/<project>/cases.json` per project. See
`tcga2hf.clinical.TOP_LEVEL_FIELDS` and `tcga2hf.clinical.EXPANSIONS`
for the exact field/expand lists.

### `POST /files` → `masked_somatic_mutation`, `gene_expression_quantification`, `files` tables

One POST per (project, modality). All requests share
`cases.project.project_id = <project>` AND `access = open`. The
remaining clauses lock the format / experimental strategy / workflow
type so future GDC additions can't silently ship a different pipeline
under the same `data_type`:

| Table | data_type | data_format | data_category | experimental_strategy | analysis.workflow_type |
|---|---|---|---|---|---|
{modality_table}

The combined `/files` responses (one per modality per project) feed the
`files` table. See `tcga2hf.genomic.MODALITY_FILTERS` and
`tcga2hf.genomic.FILE_FIELDS` for the full request payload.

### `POST /data` → file bytes

UUIDs returned by `/files` are batched (≤50 per request) into `POST
/data`; the response is a tar.gz of those files. Mutations files are
parsed row-by-row into the `masked_somatic_mutation` config; expression
files are parsed gene-row-by-gene-row into the
`gene_expression_quantification` config. See
`tcga2hf.gdc.bulk_download` for the batching and retry logic.

### Provenance pinned per build

- `GET /status` → `data_release` / `tag` / `commit` saved in each
  project's `gdc_status.json`.
- `GET /v0/submission/_dictionary/_all` → schema dictionary snapshot
  saved alongside the raw data; its SHA-256 is recorded in
  `gdc_status.json`.

[src-clinical]: https://github.com/galtay/tcga2hf/blob/main/src/tcga2hf/clinical.py
[src-genomic]: https://github.com/galtay/tcga2hf/blob/main/src/tcga2hf/genomic.py
"""


_GDC_REFERENCES = """\
## GDC references

- [Data dictionary][gdc-dict] (every entity + field definition)
- [Biospecimen Encyclopedia](https://docs.gdc.cancer.gov/Encyclopedia/pages/Biospecimen/)
- [MAF format spec](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)
- [Gene Expression Quantification spec](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
- [Sample Type codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
- [TCGA Barcode reference](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
"""


_LICENSE_AND_REDISTRIBUTION = """\
## License & redistribution

Per the [NCI GDC Data Analysis Policy](https://gdc.cancer.gov/analyze-data/data-analysis-policies):

> The GDC itself places no restrictions (other than attempts at reidentification)
> on analysis or publication of open access data provided through the GDC Data Portal.

Per the [NCI TCGA citation page](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing):

> Moratoria on all cancer types are now lifted and all TCGA data are available
> without restrictions on their use in publications or presentations.

Per the [GDC Data Access Processes and Tools page](https://gdc.cancer.gov/access-data/data-access-processes-and-tools):

> Open access data generally includes high level genomic data that is not
> individually identifiable, as well as most clinical and all biospecimen data
> elements.

## Restrictions on use

> Users of any data provided by GDC, whether open or controlled access, agree
> not to attempt to reidentify any individual participant in any study
> represented by GDC data, for any purpose whatever.
> ([source](https://gdc.cancer.gov/analyze-data/data-analysis-policies))

## Required acknowledgement

If you publish or present results derived from this dataset, include the
[NCI-required TCGA acknowledgement](https://www.cancer.gov/ccg/research/genome-sequencing/tcga/using-tcga-data/citing):

> The results <published or shown> here are in whole or part based upon data
> generated by the TCGA Research Network: https://www.cancer.gov/tcga.

Suggested citations:

- Grossman, R. L., et al. (2016). Toward a Shared Vision for Cancer Genomic Data.
  *NEJM*, 375(12), 1109-1112.
- The Cancer Genome Atlas Research Network. https://www.cancer.gov/tcga
- NCI Genomic Data Commons. https://gdc.cancer.gov

Policy references:
[GDC Policies](https://gdc.cancer.gov/about-gdc/gdc-policies),
[GDC Encyclopedia — Controlled Access][controlled] (defines what is *not* in this dataset),
[NIH Genomic Data Sharing Policy](https://sharing.nih.gov/genomic-data-sharing).

[controlled]: https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access/

## Disclaimer

Prototype dataset. Schema, included projects, and column coverage are still
evolving. Re-derive from the GDC for any analysis where freshness matters.
"""


# Markdown reference-style link definitions referenced from multiple sections.
# Defined once at the bottom of the document so they resolve everywhere.
_LINK_REFS = """\
[gdc-dict]: https://docs.gdc.cancer.gov/Data_Dictionary/
[repo]: https://github.com/galtay/tcga2hf
"""


# ---------------------------------------------------------------------------
# Dynamic helpers — produce per-build header sections from project list /
# release info passed in at write time.
# ---------------------------------------------------------------------------


def _configs_yaml(projects: list[str]) -> str:
    lines = ["configs:"]
    for project in sorted(projects):
        lines.append(f"  - config_name: {project}")
        lines.append("    data_files:")
        lines.append("      - split: train")
        lines.append(f"        path: {project}/data.parquet")
    return "\n".join(lines)


def _release_md(gdc_releases: dict[str, str] | None) -> str:
    if not gdc_releases:
        return "- **GDC data release:** unknown (status file missing)"
    unique = set(gdc_releases.values())
    if len(unique) == 1:
        return f"- **GDC data release:** {next(iter(unique))}"
    lines = ["- **GDC data releases (per project):**"]
    for proj, release in sorted(gdc_releases.items()):
        lines.append(f"    - `{proj}`: {release}")
    return "\n".join(lines)


def write_card(
    processed_dir: Path,
    projects: list[str],
    gdc_releases: dict[str, str] | None = None,
) -> Path:
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    projects_md = ", ".join(f"`{p}`" for p in sorted(projects))
    configs_block = _configs_yaml(projects)
    release_md = _release_md(gdc_releases)

    frontmatter = f"""---
license: other
license_name: nih-genomic-data-sharing
license_link: https://gdc.cancer.gov/analyze-data/data-analysis-policies
pretty_name: TCGA Patients (Open Access)
tags:
  - cancer
  - tcga
  - clinical
  - genomics
{configs_block}
---
"""

    header = f"""
# TCGA Patients (Open Access)

Open-access patient data from The Cancer Genome Atlas (TCGA), pulled from the
National Cancer Institute (NCI) Genomic Data Commons (GDC). One HuggingFace
(HF) subset per TCGA project; one row per patient.

- **Projects included:** {projects_md}
- **Generated:** {timestamp}
- **Source:** GDC `/cases` endpoint, open-access tier only.
- **Schema:** derived from the [GDC Data Dictionary][gdc-dict]; the live
  dictionary the GDC was serving when the data was fetched is hashed into
  each project's `gdc_status.json` for provenance.
{release_md}
"""

    body = "\n".join(
        [
            header,
            _DATA_MODEL,
            _render_gdc_requests(projects),
            _LOADING,
            _GDC_REFERENCES,
            _LICENSE_AND_REDISTRIBUTION,
            _LINK_REFS,
        ]
    )

    out_path = processed_dir / "README.md"
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(frontmatter + body)
    return out_path


# ---------------------------------------------------------------------------
# Tabular dataset card (`gabrielaltay/tcga-tabular-open`)
# Companion to the consolidated patients card above. Reuses the same
# license / GDC references / disclaimer blocks; replaces the data-model
# section with a tabular-specific overview + SQL example.
# ---------------------------------------------------------------------------


_TABULAR_DATA_MODEL = """\
## Data model

Four tables per project. Same source data as the consolidated companion
dataset [`gabrielaltay/tcga-patients-open`][patients], different shape:
the molecular rows here are flat (one row per variant, one row per
aliquot-gene) and easy to query in SQL, while the `cases` table keeps
the GDC `cases.json` structure intact (one row per patient, with
demographic / diagnoses / follow_ups / exposures / family_histories /
samples nested under their original GDC keys).

The biospecimen hierarchy lives inside `cases.samples[]`:

```
case          one patient (TCGA-XX-1234)        ← cases row
└── sample    physical specimen at one timepoint
    └── portion    a piece of that sample for a lab process
        └── analyte    DNA or RNA extracted
            └── aliquot    a vial handed off for sequencing
```

Each (project, table) pair is one HF **config** named
`<project>_<table>` (e.g. `TCGA_LUAD_cases`,
`TCGA_LUAD_masked_somatic_mutation`); each config has the canonical
`train` split. We use a config per pair rather than splits-per-project
within one config because HF Data Studio expects splits to be
train / validation / test, and we have no such semantics.

| Table | Grain | Source |
|---|---|---|
| `cases` | 1 / patient | GDC `cases.json` verbatim (JSON keys → columns; nested fields kept) |
| `masked_somatic_mutation` | 1 / variant | open-access MAF (Mutation Annotation Format) files |
| `gene_expression_quantification` | 1 / (aliquot, gene) | RNA-Seq STAR counts TSV gene rows |
| `files` | 1 / (file, case) | per-modality `manifest.json` (file provenance / GDC metadata) |

The `cases` row mirrors the GDC `case` JSON: scalar fields stay scalar,
`demographic` stays a struct, `diagnoses` / `follow_ups` / `exposures` /
`family_histories` / `samples` stay lists of structs (with the
`samples[].portions[].analytes[].aliquots[]` biospecimen tree intact). To
walk it from SQL, use DuckDB's `UNNEST` on the array columns; for ergonomic
Python access, the consolidated [`tcga-patients-open`][patients] dataset
ships a typed `TcgaHfPatient` pydantic model with helper methods.

### Foreign keys

Every flat table (`masked_somatic_mutation`,
`gene_expression_quantification`, `files`) carries `case_id` and
`case_submitter_id` columns so users can filter and join without
roundtripping through `cases`. UUID FKs are paired with their
human-readable `*_submitter_id` form (e.g. `case_submitter_id` =
`TCGA-3X-AAVB`):

- `masked_somatic_mutation`: `case_id`, `case_submitter_id`,
  `tumor_sample_id`, `matched_normal_sample_id`, `source_file_id`
- `gene_expression_quantification`: `case_id`, `case_submitter_id`,
  `aliquot_id`, `aliquot_submitter_id`, `source_file_id`
- `files`: `case_id`, `case_submitter_id`, `file_id`

### Example: TP53 mutations in TCGA-LUAD via SQL

```sql
-- HF Data Studio's SQL console exposes each (config, split) pair as a
-- single table named `<config>_<split>` (lowercased). With config
-- `TCGA_LUAD_masked_somatic_mutation` × split `train` →
-- `tcga_luad_masked_somatic_mutation_train`.
SELECT case_submitter_id, Hugo_Symbol, HGVSp_Short, t_alt_count, t_depth
FROM tcga_luad_masked_somatic_mutation_train
WHERE Hugo_Symbol = 'TP53' AND Variant_Classification != 'Silent'
ORDER BY t_alt_count DESC
LIMIT 50;
```

Want all the data nested per patient instead? See the consolidated
[`tcga-patients-open`][patients] dataset.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
"""


_TABULAR_LOADING = """\
## Loading

One config per (project, table). Config names use underscores
throughout (e.g. `TCGA_LUAD_masked_somatic_mutation`); each config has
a single conventional `train` split. File paths on disk keep the
canonical dash-separated GDC project_id (`TCGA-LUAD/...`).

### `datasets`

```python
from datasets import load_dataset

# LUAD cases (one row per patient, nested entities preserved):
luad_cases = load_dataset(
    "gabrielaltay/tcga-tabular-open", "TCGA_LUAD_cases", split="train"
)
luad_muts = load_dataset(
    "gabrielaltay/tcga-tabular-open",
    "TCGA_LUAD_masked_somatic_mutation",
    split="train",
)
```

### DuckDB (remote, via `hf://`)

DuckDB can query the parquets in this repo directly without `datasets`,
using the [`hf://` URL scheme][duckdb-hf]. Row groups + page indexes
mean column-projection + filter pushdown skip most of the data on
selective queries.

```python
import duckdb, os
con = duckdb.connect()
con.execute(f"CREATE SECRET hf (TYPE huggingface, TOKEN '{os.environ['HF_TOKEN']}')")
con.sql('''
    SELECT case_submitter_id, Hugo_Symbol, HGVSp_Short, t_alt_count, t_depth
    FROM "hf://datasets/gabrielaltay/tcga-tabular-open/TCGA-LUAD/masked_somatic_mutation/data.parquet"
    WHERE Hugo_Symbol = 'TP53' AND Variant_Classification != 'Silent'
    ORDER BY t_alt_count DESC LIMIT 10
''').show()
```

The `CREATE SECRET` step is only needed while the dataset is private —
DuckDB picks up `HF_TOKEN` from the env automatically once it's public.

[duckdb-hf]: https://huggingface.co/docs/hub/datasets-duckdb-sql
"""


def _tabular_configs_yaml(projects: list[str], tables: list[str]) -> str:
    """One config per (project, table); one canonical `train` split each.

    HF Data Studio expects splits to be train / validation / test (its
    dataset-server warns when a config has more than three). We don't have
    train/val/test slices, so we model each (project, table) pair as its
    own config and leave the split slot at the conventional `train`.

    Config names are `<project>_<table>` with dashes in the project_id
    normalized to underscores (e.g. `TCGA_CHOL_cases`) so the config name
    is a valid identifier in HF Data Studio's SQL console without further
    sanitization. File paths on disk keep the canonical dash form
    (`TCGA-CHOL/cases/data.parquet`).
    """
    lines = ["configs:"]
    for project in sorted(projects):
        config_proj = project.replace("-", "_")
        for table in tables:
            lines.append(f"  - config_name: {config_proj}_{table}")
            lines.append("    data_files:")
            lines.append("      - split: train")
            lines.append(f"        path: {project}/{table}/data.parquet")
    return "\n".join(lines)


def write_tabular_card(
    processed_dir: Path,
    projects: list[str],
    tables: list[str],
    gdc_releases: dict[str, str] | None = None,
) -> Path:
    """Write the dataset card for the tabular HF dataset.

    `tables` is the ordered list of table names emitted per project (the
    keys of `tcga2hf.schema.TABULAR_TABLES`). Same license / GDC refs /
    disclaimer blocks as the patients card.
    """
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    projects_md = ", ".join(f"`{p}`" for p in sorted(projects))
    configs_block = _tabular_configs_yaml(projects, tables)
    release_md = _release_md(gdc_releases)

    frontmatter = f"""---
license: other
license_name: nih-genomic-data-sharing
license_link: https://gdc.cancer.gov/analyze-data/data-analysis-policies
pretty_name: TCGA Tabular (Open Access)
tags:
  - cancer
  - tcga
  - clinical
  - genomics
{configs_block}
---
"""

    header = f"""
# TCGA Tabular (Open Access)

Open-access TCGA data from the National Cancer Institute (NCI) Genomic
Data Commons (GDC), reshaped as one HuggingFace (HF) subset per (project,
table). Companion to the consolidated nested-patient dataset
[`gabrielaltay/tcga-patients-open`][patients].

- **Projects included:** {projects_md}
- **Tables per project:** {", ".join(f"`{t}`" for t in tables)}
- **Generated:** {timestamp}
- **Source:** GDC `/cases` endpoint + per-modality MAFs and TSVs,
  open-access tier only.
- **Schema:** derived from the [GDC Data Dictionary][gdc-dict].
{release_md}
"""

    body = "\n".join(
        [
            header,
            _TABULAR_DATA_MODEL,
            _render_gdc_requests(projects),
            _TABULAR_LOADING,
            _GDC_REFERENCES,
            _LICENSE_AND_REDISTRIBUTION,
            _LINK_REFS,
        ]
    )

    out_path = processed_dir / "README.md"
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(frontmatter + body)
    return out_path
