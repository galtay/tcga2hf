from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from tcga2hf.clinical import EXPANSIONS as _CLINICAL_EXPAND
from tcga2hf.clinical import TOP_LEVEL_FIELDS as _CASE_FIELDS_REQUESTED
from tcga2hf.genomic import FILE_FIELDS as _FILE_FIELDS_REQUESTED
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
    """Render the exact GDC requests this build issues, populated from the
    real constants in `tcga2hf.clinical` and `tcga2hf.genomic` so the docs
    can't drift from the code. Templated per dataset build with the project
    list it was actually fetched for.
    """
    project_list = sorted(projects)

    # /cases payload — the one the clinical fetcher posts. Keys mirror the
    # GDC REST shape: filters + fields + expand + size + from + format.
    cases_payload = {
        "filters": {
            "op": "in",
            "content": {
                "field": "project.project_id",
                "value": project_list,
            },
        },
        "fields": ",".join(_CASE_FIELDS_REQUESTED),
        "expand": ",".join(_CLINICAL_EXPAND),
        "format": "JSON",
        "size": 200,
        "from": 0,
    }
    cases_json = json.dumps(cases_payload, indent=2)

    # /files payload — same shape but per (project, data_type) pair. Show
    # one example with the LUAD-style filter; the actual builds run this
    # once per project × per data_type. The clauses are templated from
    # `MODALITY_FILTERS` so the example reflects whatever modality-locking
    # filters the build is currently using.
    example_project = project_list[-1] if project_list else "<PROJECT>"
    example_data_type = "Masked Somatic Mutation"
    files_clauses: list[dict[str, Any]] = [
        {
            "op": "=",
            "content": {
                "field": "cases.project.project_id",
                "value": example_project,
            },
        },
        {"op": "=", "content": {"field": "access", "value": "open"}},
        {"op": "=", "content": {"field": "data_type", "value": example_data_type}},
    ]
    for field, value in _MODALITY_FILTERS.get(example_data_type, {}).items():
        files_clauses.append({"op": "=", "content": {"field": field, "value": value}})
    files_payload = {
        "filters": {"op": "and", "content": files_clauses},
        "fields": ",".join(_FILE_FIELDS_REQUESTED),
        "format": "JSON",
        "size": 500,
        "from": 0,
    }
    files_json = json.dumps(files_payload, indent=2)

    # Render the per-modality filter table as markdown rows so consumers
    # see the full set of (data_type, additional clauses) without having
    # to read the source.
    modality_rows: list[str] = []
    for dtype, extras in _MODALITY_FILTERS.items():
        clauses = [f"`data_type={dtype!r}`"] + [f"`{k}={v!r}`" for k, v in extras.items()]
        modality_rows.append(f"- {' AND '.join(clauses)}")
    modalities_md = "\n".join(modality_rows)

    return f"""\
## How this dataset is built

Three GDC REST endpoints feed every row in this dataset. Both endpoints
and payloads below are reproduced verbatim from the build code; the
on-disk `gdc_status.json` for each project additionally pins the GDC
data release and dictionary SHA-256 the data was fetched against.

### `POST /cases` — clinical entities

One paginated POST per build, returning the nested case JSON for every
patient in the requested projects.

```json
{cases_json}
```

The full nested response (top-level scalars + `demographic` + `diagnoses[]`
with `treatments[]` inside + `follow_ups[]` + `exposures[]` +
`family_histories[]` + `samples[]` with `portions[].analytes[].aliquots[]`
inside) is written to `<data-dir>/raw/<project>/cases.json` and feeds the
`cases` config (and the consolidated patient row).

### `POST /files` — molecular file discovery

One POST per (project, modality) pair, listing every open-access file of
that modality. Each modality is identified by a `data_type` plus a small
set of clauses that lock the format / experimental strategy / workflow
type so future GDC additions can't silently ship under the same
`data_type`:

{modalities_md}

Below is the full request shape for `Masked Somatic Mutation` against
`{example_project}`; the build runs the same shape with the
`Gene Expression Quantification` clauses too. Responses are written to
`<data-dir>/raw/<project>/<modality>/manifest.json` and feed the `files`
config.

```json
{files_json}
```

### `POST /data` — file bytes

For every file UUID returned by `/files`, the build batches up to 50
UUIDs per request:

```json
{{"ids": ["<file_uuid_1>", "<file_uuid_2>", "..."]}}
```

The response is a tar.gz of the files. Mutations files are parsed
row-by-row into the `masked_somatic_mutation` config; expression files
are parsed gene-row-by-gene-row into the
`gene_expression_quantification` config. (See `tcga2hf.gdc.bulk_download`
for the batching + retry logic.)

### Other GDC calls captured for provenance

- `GET /status` — once at the top of every build; `data_release` /
  `tag` / `commit` are saved in each project's `gdc_status.json`.
- `GET /v0/submission/_dictionary/_all` — once at the top of every build;
  written to `<data-dir>/raw/gdc_dictionary.<major>.<minor>.json` and its
  SHA-256 is recorded in `gdc_status.json`. The static schema in
  `tcga2hf.schema` is regenerated from this dictionary; see the package
  `scripts/regenerate_clinical_fields.py`.
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

Each table is one HF
**config** named after the table; each TCGA project is one **split**
inside that config. (HF's `split` slot is arbitrary string-typed; we have
no train/val/test semantics so we repurpose it to identify the project.)

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
-- single table named `<config>_<split>` (lowercased). So
-- `masked_somatic_mutation` × `TCGA_LUAD` →
-- `masked_somatic_mutation_tcga_luad`.
SELECT case_submitter_id, Hugo_Symbol, HGVSp_Short, t_alt_count, t_depth
FROM masked_somatic_mutation_tcga_luad
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

One config per table; the TCGA project is the split. Splits use
underscores (`TCGA_LUAD`) since HF's dataset-server validator forbids
dashes in split names; file paths on disk keep the canonical
dash-separated GDC project_id (`TCGA-LUAD/`).

### `datasets`

```python
from datasets import load_dataset

# Just LUAD's cases table (one row per patient, nested entities preserved):
luad_cases = load_dataset("gabrielaltay/tcga-tabular-open", "cases", split="TCGA_LUAD")

# All projects' mutations (concatenated across splits):
all_muts = load_dataset(
    "gabrielaltay/tcga-tabular-open", "masked_somatic_mutation", split="all"
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
    """One config per table; one split per project inside each config.

    HF `split` is just an arbitrary string label — we have no train/val/test
    semantics, so the slot is repurposed to identify the TCGA project. This
    collapses the cross-product (#projects × #tables configs) down to one
    config per table with the project list folded into splits, and matches
    the actual data shape: each project/table combo carries the same schema.

    HF's dataset-server validator rejects dashes in split names, so the
    GDC `project_id` (e.g. `TCGA-CHOL`) is normalized to `TCGA_CHOL` for
    the split slot only — file paths on disk keep the canonical
    dash-separated GDC name.
    """
    lines = ["configs:"]
    for table in tables:
        lines.append(f"  - config_name: {table}")
        lines.append("    data_files:")
        for project in sorted(projects):
            split = project.replace("-", "_")
            lines.append(f"      - split: {split}")
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
