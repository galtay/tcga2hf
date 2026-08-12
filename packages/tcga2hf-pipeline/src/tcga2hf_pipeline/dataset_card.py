from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from tcga2hf_pipeline.genomic import MODALITY_FILTERS as _MODALITY_FILTERS

# ---------------------------------------------------------------------------
# Both dataset cards (`gabrielaltay/tcga-patients-open` and
# `gabrielaltay/tcga-tabular-open`) are presented as two views of the same
# underlying data, not a primary + secondary pair. They share the same
# section order; per-view variation lives only in the sentences/tables
# that genuinely depend on the shape (one-row-per-patient vs
# one-table-per-source).
#
# Section order (both cards):
#   1. Header (mostly shared)
#   2. Data model (combined with how-it-was-built; mostly shared)
#   3. Survival endpoints supplement (mostly shared)
#   4. Loading (per-view, both use load_dataset)
#   5. GDC references (shared)
#   6. License & redistribution (shared)
#   7. Restrictions on use (shared)
#   8. Required acknowledgement (shared)
#   9. Disclaimer (shared)
# ---------------------------------------------------------------------------


# ===========================================================================
# Shared content blocks (byte-identical text in both cards)
# ===========================================================================


_SHARED_PROVENANCE = """\
### Provenance pinned per build

- `GET /status` → `data_release` / `tag` / `commit` saved in each
  project's `gdc_status.json`.
- `GET /v0/submission/_dictionary/_all` → schema dictionary snapshot
  saved alongside the raw data; its SHA-256 is recorded in
  `gdc_status.json`.

See the [repository][repo] for full request payloads, filter clauses,
and the build pipeline source.
"""


def _shared_filters_table() -> str:
    """Canonical `/files` filter clauses per source. Same in both cards.

    Templated from `genomic.MODALITY_FILTERS` so molecular filter docs
    can't drift from code; the Clinical Supplement row is appended
    manually since `clinical_supplement.py` uses its own filter set.
    """
    header = (
        "| data_type | data_format | data_category | "
        "experimental_strategy | analysis.workflow_type |"
    )
    sep = "|---|---|---|---|---|"
    rows: list[str] = []
    for dtype, extras in _MODALITY_FILTERS.items():
        cols = [f"`{dtype}`"]
        for col_field in (
            "data_format",
            "data_category",
            "experimental_strategy",
            "analysis.workflow_type",
        ):
            cols.append(f"`{extras.get(col_field, '')}`")
        rows.append("| " + " | ".join(cols) + " |")
    rows.append(
        "| `Clinical Supplement` | `bcr biotab` | `Clinical` | | |"
    )
    return f"{header}\n{sep}\n" + "\n".join(rows)


# ===========================================================================
# Data model section (combined with how-it-was-built)
# ===========================================================================


def _data_model(*, consolidated: bool) -> str:
    """Section 2: Data model & how this dataset was built.

    Same shape on both cards. Per-view content slots in via the
    `consolidated` branch (intro sentence + per-view source-mapping
    table + per-view "specific to this view" prose).
    """
    if consolidated:
        per_view_table = _PATIENT_VIEW_MAPPING
        per_view_specifics = _PATIENT_VIEW_SPECIFICS
    else:
        per_view_table = _TABULAR_VIEW_MAPPING
        per_view_specifics = _TABULAR_VIEW_SPECIFICS

    return f"""\
## Data model

### Where the data comes from

Three sources feed each project's data, all open-access:

- **Case-level clinical structure** — fetched from the GDC `/cases`
  endpoint, returning the full nested case JSON (demographic + diagnoses
  → treatments + follow_ups + exposures + family_histories + samples →
  portions → analytes → aliquots). The biospecimen subtree on each case:

  ```
  case          one patient (TCGA-XX-1234)
  └── sample    physical specimen taken from the patient at one timepoint
                (Primary Tumor, Solid Tissue Normal, Blood Derived Normal, ...)
      └── portion    a piece of that sample for a specific lab process
          └── analyte    extracted material of one type (DNA or RNA)
              └── aliquot    a vial of that analyte handed off for sequencing
  ```

- **Per-modality files** — discovered via `/files` (filtered by the
  clauses in the table below) and downloaded via `/data`. Each
  combination locks one `data_type` to a specific GDC pipeline so a
  future GDC addition can't quietly substitute a different pipeline
  under the same `data_type`. This covers both the molecular modalities
  and the scanned Pathology Report PDFs, which are carried verbatim —
  no text extraction is applied, so consumers can run whichever parser
  they trust against the original document.
- **BCR Clinical Supplement biotabs** — original Biospecimen Core
  Resource (BCR) clinical forms shipped as per-project TSVs (one per
  form: patient, follow_up, nte, drug, radiation, etc.). The harmonized
  `/cases` endpoint drops or under-populates a number of clinical fields
  the BCR-original biotabs preserve. The schema varies by cancer type
  (e.g. BLCA's BCG-response columns don't exist in CHOL's hepatic-marker
  forms), so each project's biotabs ship only the columns they actually
  carry. Discovered the same way (`/files` then `/data`) — see the
  filter table below.

### Source data filters (canonical)

Same in both views of the dataset; each row locks the `/files` query
for one source:

{_shared_filters_table()}

### How each source appears in this view

{per_view_table}

### Specific to this view

{per_view_specifics}

{_SHARED_PROVENANCE}
"""


_PATIENT_VIEW_MAPPING = """\
| Source | Where it lands |
|---|---|
| GDC `/cases` | nested fields on each patient row (`demographic`, `diagnoses`, `follow_ups`, `exposures`, `family_histories`, `samples`); `gdc_portal_url` link added |
| Masked Somatic Mutation MAFs | `samples_masked_somatic_mutation` array on each patient row (sample FKs resolved alongside GDC's aliquot UUIDs) |
| Gene Expression Quantification | `samples_gene_expression_quantification` array on each patient row (`stranded_first` / `stranded_second` dropped — GDC harmonizes as unstranded) |
| BCR Clinical Supplements | `clinical_supplement` struct on each patient row, with sub-fields `patient` (1 dict) and `follow_ups` / `ntes` / `drugs` / `radiations` / `ablations` / `omfs` (lists of dicts). Sub-fields with no data for the project are omitted. |
| Pathology Reports | `samples_pathology_report` array on each patient row; `pdf_bytes` holds the scanned PDF verbatim, joined to its sample via `pathology_report_uuid` |
"""


_TABULAR_VIEW_MAPPING = """\
| Source | Where it lands |
|---|---|
| GDC `/cases` | `cases` table — one row per patient with the GDC `case` JSON structure preserved as nested struct columns (`demographic`, `diagnoses[]`, `follow_ups[]`, `samples[]`, ...); `gdc_portal_url` link added |
| Masked Somatic Mutation MAFs | `masked_somatic_mutation` table — one row per variant (sample FKs resolved alongside GDC's aliquot UUIDs) |
| Gene Expression Quantification | `gene_expression_quantification` table — one row per (aliquot, gene); `stranded_first` / `stranded_second` dropped (GDC harmonizes as unstranded) |
| BCR Clinical Supplements | `clinical_supplement_*` tables — one per BCR form (patient, follow_up, nte, drug, radiation, ablation, omf) |
| Pathology Reports | `pathology_report` table — one row per report, with the scanned PDF verbatim in `pdf_bytes` |
| Per-modality manifests | `files` table — one row per (file, case) with `file_id` / `md5sum` / `data_type` / size; useful for joining a row back to its source file or replaying a specific modality fetch |
"""


_PATIENT_VIEW_SPECIFICS = """\
- Convenience: each row carries `samples_<modality>` array columns so
  you can column-project just the molecular data you need without
  walking the nested GDC entities.
- Loading: the [`tcga2hf` package][repo] ships a typed `TcgaHfPatient`
  pydantic model that mirrors this schema and adds convenience joins
  (tumor/normal pairs, mutations-by-gene, expression-by-gene,
  longitudinal timeline).
"""


_TABULAR_VIEW_SPECIFICS = """\
- Joining: every flat row carries `case_submitter_id` (and
  `aliquot_submitter_id` where applicable) for direct joins back to
  `cases` without re-resolving UUIDs.
- Each (project, table) pair is one HuggingFace **config** named
  `<project>_<table>` (e.g. `TCGA_LUAD_cases`).
- BCR fields in the `clinical_supplement_*` tables are all typed as
  strings, with sentinel values like `[Not Available]` preserved
  verbatim.
"""


# ===========================================================================
# Survival endpoints supplement (mostly shared; one location sentence varies)
# ===========================================================================


def _shared_survival_endpoints(*, consolidated: bool) -> str:
    """Section 3: re-derived survival endpoints supplement.

    Forward-compatible phrasing — additional supplements in the future
    can be added without breaking the framing here ("a supplement..."
    not "the only supplement...").
    """
    if consolidated:
        location = (
            "Each patient row carries a top-level **`survival_derived` struct** "
            "with eight sub-fields: `os_event` / `os_time`, `dss_event` / "
            "`dss_time`, `pfi_event` / `pfi_time`, `dfi_event` / `dfi_time`."
        )
    else:
        location = (
            "Surfaced as a standalone **`survival_derived` table** (one row per "
            "patient, joined to `cases` on `case_submitter_id`) with eight "
            "columns: `os_event` / `os_time`, `dss_event` / `dss_time`, "
            "`pfi_event` / `pfi_time`, `dfi_event` / `dfi_time`."
        )
    return f"""\
## Survival endpoints (`survival_derived`)

We have provided a supplement to the GDC source data: re-derived
survival endpoints — Overall Survival (OS), Disease-Specific Survival
(DSS), Progression-Free Interval (PFI), Disease-Free Interval (DFI) —
following the algorithm published by **Liu et al. 2018**
([DOI 10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052)).

{location} `*_event` is 0/1 (event observed vs censored); `*_time` is
days from `index_date` (TCGA: diagnosis date). DFI is null for SKCM /
THYM / UVM / LAML — Liu specifies no DFI for those tumor types.

We've reimplemented Liu's method against the current TCGA data and find
broad agreement with the original curated CDR. Differences exist and are
expected: this is a newer release of the underlying GDC data, so
re-curated clinical values, post-2018 patient additions, and schema
migrations all contribute to the gap. This work is evolving; see the
[repository][repo] for the full reproduction report and per-endpoint
methodology.

**Why we don't ship Liu's curated 2018 values directly:** the CDR is a
frozen 2018 snapshot derived from a since-modified GDC release.
Including those values would lock in irreproducible source-data drift.
We re-derive on every build, so the values reflect the current GDC and
are reproducible from this dataset's other tables alone.
"""


# ===========================================================================
# Loading examples (per-view; both use load_dataset)
# ===========================================================================


_PATIENT_LOADING = """\
## Loading

```python
from datasets import load_dataset

# One config per TCGA project.
luad = load_dataset("gabrielaltay/tcga-patients-open", "TCGA-LUAD")
```

Each row is one patient with the full GDC `case` structure nested
in-place plus the `survival_derived` struct.
"""


_TABULAR_LOADING = """\
## Loading

```python
from datasets import load_dataset

# One config per (project, table).
# Config names are <project>_<table> with dashes replaced by underscores.
luad_cases = load_dataset("gabrielaltay/tcga-tabular-open", "TCGA_LUAD_cases")
luad_muts = load_dataset(
    "gabrielaltay/tcga-tabular-open", "TCGA_LUAD_masked_somatic_mutation"
)
```
"""


# ===========================================================================
# GDC references / License / Disclaimer (all shared)
# ===========================================================================


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

**This project is not affiliated with the NCI, GDC, or the TCGA Research
Network.** It is an experimental open-source pipeline that may change
significantly between versions. Pipeline source: [`galtay/tcga2hf`][repo].
"""


# Markdown reference-style link definitions referenced from multiple sections.
# Defined once at the bottom of the document so they resolve everywhere.
_LINK_REFS = """\
[gdc-dict]: https://docs.gdc.cancer.gov/Data_Dictionary/
[repo]: https://github.com/galtay/tcga2hf
[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
"""


# ===========================================================================
# Header builder (mostly shared; one sentence varies per view)
# ===========================================================================


def _header(*, consolidated: bool, timestamp: str, release_md: str) -> str:
    if consolidated:
        title = "TCGA Patients (Open Access)"
        view_sentence = (
            "**This view presents one HuggingFace subset per TCGA project**, "
            "with one row per patient. See the [`tcga-tabular-open`][tabular] "
            "companion for a per-table view of the same underlying data."
        )
    else:
        title = "TCGA Tabular (Open Access)"
        view_sentence = (
            "**This view presents one HuggingFace subset per (project, table)**. "
            "See the [`tcga-patients-open`][patients] companion for a per-row "
            "view of the same underlying data."
        )
    return f"""
# {title}

Open-access TCGA data from the NCI Genomic Data Commons (GDC). Covers
all 33 TCGA projects.

{view_sentence}

- **Generated:** {timestamp}
- **Schema:** derived from the [GDC Data Dictionary][gdc-dict].
{release_md}
"""


# ===========================================================================
# Dynamic helpers (per-build header inputs + configs YAML)
# ===========================================================================


def _configs_yaml(processed_dir: Path, projects: list[str]) -> str:
    """One config per project whose parquet actually exists on disk.

    Filesystem-aware for the same reason `_tabular_configs_yaml` is: a
    config pointing at a missing path makes HF Data Studio emit "file not
    found" rather than simply omitting the subset.
    """
    lines = ["configs:"]
    for project in sorted(projects):
        if not (processed_dir / project / "data.parquet").exists():
            continue
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


def _tabular_configs_yaml(
    processed_dir: Path,
    projects: list[str],
    tables: list[str],
) -> str:
    """One config per (project, table) where the parquet actually exists.

    HF Data Studio expects splits to be train / validation / test (its
    dataset-server warns when a config has more than three). We don't have
    train/val/test slices, so we model each (project, table) pair as its
    own config and leave the split slot at the conventional `train`.

    Config names are `<project>_<table>` with dashes in the project_id
    normalized to underscores (e.g. `TCGA_CHOL_cases`) so the config name
    is a valid identifier in HF Data Studio's SQL console without further
    sanitization. File paths on disk keep the canonical dash form
    (`TCGA-CHOL/cases/data.parquet`).

    Filesystem-aware: `clinical_supplement_*` tables are emitted per
    project only when the corresponding biotab form has any rows (e.g.
    LAML has `clinical_supplement_patient` but no follow_up; only LIHC has
    `clinical_supplement_ablation`). Skipping non-existent paths keeps HF
    Data Studio from emitting "file not found" warnings.
    """
    lines = ["configs:"]
    for project in sorted(projects):
        config_proj = project.replace("-", "_")
        for table in tables:
            parquet_path = processed_dir / project / table / "data.parquet"
            if not parquet_path.exists():
                continue
            lines.append(f"  - config_name: {config_proj}_{table}")
            lines.append("    data_files:")
            lines.append("      - split: train")
            lines.append(f"        path: {project}/{table}/data.parquet")
    return "\n".join(lines)


# ===========================================================================
# Per-view writers — choose configs YAML, header, loading; body identical order
# ===========================================================================


def write_card(
    processed_dir: Path,
    projects: list[str],
    gdc_releases: dict[str, str] | None = None,
) -> Path:
    """Write the patients dataset card."""
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    configs_block = _configs_yaml(processed_dir, projects)
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

    body = "\n".join(
        [
            _header(consolidated=True, timestamp=timestamp, release_md=release_md),
            _data_model(consolidated=True),
            _shared_survival_endpoints(consolidated=True),
            _PATIENT_LOADING,
            _GDC_REFERENCES,
            _LICENSE_AND_REDISTRIBUTION,
            _LINK_REFS,
        ]
    )

    out_path = processed_dir / "README.md"
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(frontmatter + body)
    return out_path


def write_tabular_card(
    processed_dir: Path,
    projects: list[str],
    tables: list[str],
    gdc_releases: dict[str, str] | None = None,
) -> Path:
    """Write the tabular dataset card.

    `tables` is the ordered list of table names emitted per project (the
    keys of `tcga2hf.schema.TABULAR_TABLES` plus the `clinical_supplement_*`
    flex-schema tables). Same body order as the patients card; only
    configs YAML + header view-sentence + loading example differ.
    """
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    configs_block = _tabular_configs_yaml(processed_dir, projects, tables)
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

    body = "\n".join(
        [
            _header(consolidated=False, timestamp=timestamp, release_md=release_md),
            _data_model(consolidated=False),
            _shared_survival_endpoints(consolidated=False),
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
