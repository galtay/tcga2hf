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
| Allele-specific Copy Number Segment | `samples_allele_specific_copy_number_segment` array — one record per (aliquot, caller) with segments as index-aligned arrays; filter on `workflow_type` |
| Masked Copy Number Segment | `samples_masked_copy_number_segment` array — one record per aliquot, log2 ratios in `segment_mean` |
| miRNA Expression Quantification | `samples_mirna_expression_quantification` array — one record per aliquot, ~1,881 miRNAs as index-aligned arrays |
| Protein Expression Quantification | `samples_protein_expression_quantification` array — one record per **portion** (not aliquot), ~487 antibodies |
| BCR Biospecimen Supplements | `biospecimen_supplement` struct on each patient row, with list-valued sub-fields (`sample`, `portion`, `analyte`, `aliquot`, `slide`, `protocol`, `ssf_*`, ...). Sub-fields with no data for the project are omitted. |
| MSigDB gene sets + RNA-Seq | `samples_ssgsea_<collection>` array columns — pathway activity per aliquot; see the ssGSEA section below |
"""


_TABULAR_VIEW_MAPPING = """\
| Source | Where it lands |
|---|---|
| GDC `/cases` | `cases` table — one row per patient with the GDC `case` JSON structure preserved as nested struct columns (`demographic`, `diagnoses[]`, `follow_ups[]`, `samples[]`, ...); `gdc_portal_url` link added |
| Masked Somatic Mutation MAFs | `masked_somatic_mutation` table — one row per variant (sample FKs resolved alongside GDC's aliquot UUIDs) |
| Gene Expression Quantification | `gene_expression_quantification` table — one row per (aliquot, gene); `stranded_first` / `stranded_second` dropped (GDC harmonizes as unstranded) |
| BCR Clinical Supplements | `clinical_supplement_*` tables — one per BCR form (patient, follow_up, nte, drug, radiation, ablation, omf) |
| Pathology Reports | `pathology_report` table — one row per report, with the scanned PDF verbatim in `pdf_bytes` |
| Allele-specific Copy Number Segment | `allele_specific_copy_number_segment` table — one row per segment, integer total / major / minor copy number from ASCAT2, ASCAT3 and AscatNGS (filter on `workflow_type`) |
| Masked Copy Number Segment | `masked_copy_number_segment` table — one row per DNAcopy segment, log2 ratio in `segment_mean`, germline CNVs masked out |
| miRNA Expression Quantification | `mirna_expression_quantification` table — one row per (aliquot, miRBase v21 mature miRNA) |
| Protein Expression Quantification | `protein_expression_quantification` table — one row per (portion, antibody) RPPA measurement |
| BCR Biospecimen Supplements | `biospecimen_supplement_*` tables — one per BCR form (sample, portion, analyte, aliquot, slide, protocol, ssf_*, ...) |
| MSigDB gene sets + RNA-Seq | `ssgsea_scores_*` / `ssgsea_stats_*` tables — pathway activity per aliquot; see the ssGSEA section below |
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


_MSIGDB_COLLECTION_URL = "https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp"
_MSIGDB_LICENCE_URL = "https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp"

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
# ssGSEA pathway activity (tabular view only)
# ===========================================================================


def _ssgsea_section(*, consolidated: bool) -> str:
    """Section describing the ssGSEA data, with links to every collection.

    Shared by both cards; only the shape sentence and the normalization
    pointer differ, since the reference distributions live in the tabular
    view's `ssgsea_stats_*` tables and have no consolidated counterpart.

    Templated from `msigdb.COLLECTIONS` so the documented collections and
    md5s cannot drift from what the pipeline actually fetches.
    """
    from tcga2hf.schema import SSGSEA_COLLECTIONS

    from tcga2hf_pipeline import msigdb as _msigdb

    rows = []
    for key in SSGSEA_COLLECTIONS:
        c = _msigdb.COLLECTIONS[key]
        rows.append(
            f"| [`{key}`]({_MSIGDB_COLLECTION_URL}) | {c.description} | "
            f"`{c.file_name}` | `{c.md5}` |"
        )
    table = "\n".join(rows)
    if consolidated:
        shape = (
            "Single-sample gene set enrichment scores for every RNA-Seq\n"
            "aliquot, in one `samples_ssgsea_<collection>` column per MSigDB\n"
            "collection. Each entry is one scored aliquot, with `pathway`,\n"
            "`pathway_url`, the gene counts and `score_raw` as index-aligned\n"
            "arrays — the same struct-of-arrays shape\n"
            "`samples_gene_expression_quantification` uses."
        )
        normalization = (
            "The reference distributions needed to normalize live in the\n"
            "[tabular view][tabular]'s `ssgsea_stats_<collection>` tables, which\n"
            "have no counterpart here: they are cohort-level aggregates, and\n"
            "this view is one row per patient. Load them from there to\n"
            "z-score against a population, or to recover GSVA's divisor as\n"
            "`MAX(max) - MIN(min)` over their `pan_cancer` rows. The scores\n"
            "themselves are identical in both views."
        )
    else:
        shape = (
            "Single-sample gene set enrichment scores for every RNA-Seq\n"
            "aliquot, one `ssgsea_scores_<collection>` table per MSigDB\n"
            "collection plus a matching `ssgsea_stats_<collection>` table of\n"
            "reference distributions."
        )
        normalization = (
            "The `ssgsea_stats_*` tables carry that composition-dependent\n"
            "information instead. Each project ships its own reference\n"
            "distributions **and** the pan-cancer ones, so you can normalize\n"
            "without scanning every config:\n\n"
            "- **GSVA-equivalent normalization**: divide by\n"
            "  `MAX(max) - MIN(min)` over the `pan_cancer` rows.\n"
            "- **z-score against a reference population**:\n"
            "  `(score_raw - mean) / sd` for the `population` and\n"
            "  `sample_type` you care about.\n\n"
            "Everything in `ssgsea_stats_*` is derivable from\n"
            "`ssgsea_scores_*` by aggregation; it is a materialized\n"
            "convenience view, not independent evidence."
        )
    return f"""\
## Pathway activity (ssGSEA)

{shape}

`pathway_url` links to the authoritative MSigDB definition of each gene
set — so what a score means is one click away from the score itself.

### Collections

Pinned to MSigDB **{_msigdb.MSIGDB_VERSION}** and verified by md5, because
gene-set membership changes between MSigDB releases and feeds directly
into every score.

| collection | contents | file | md5 |
|---|---|---|---|
{table}

MSigDB is released under CC BY 4.0; some constituent collections carry
extra restrictions, so we ship only collections we can redistribute
scores from. See the [MSigDB licence terms]({_MSIGDB_LICENCE_URL}).

### Method

Barbie et al. (2009) ssGSEA as implemented by Bioconductor GSVA,
transcribed to Python and validated against GSVA 2.6.6 to floating-point
noise (Pearson/Spearman 1.0000000000, max relative difference 4.8e-13).
`alpha=0.25`, gene sets filtered to a minimum of 10 genes **after**
mapping onto the expression matrix; no maximum size.

Scored on `tpm_unstranded` over a gene universe of protein-coding genes
plus the functional immunoglobulin / T-cell-receptor segments. That last
inclusion matters for tumour-immune biology: GENCODE gives Ig/TCR
segments their own biotypes, and without them MSigDB's B-cell-receptor
and complement pathways match as little as 8% of their genes. Note that
V/D/J segments are somatically rearranged, so their expression reports
lymphocyte infiltration rather than regulation of a fixed locus.

Because ssGSEA weights **ranks** rather than expression values, any
strictly monotonic transform of the input leaves scores unchanged — there
is no reason to log-transform before scoring.

### Why `score_raw`, and how to normalize

`score_raw` is the only score column, and it is a property of its own
sample: it does not depend on which other samples or gene sets were
scored alongside it. GSVA's optional normalization divides by the range
of the entire score matrix, which would make every value depend on cohort
and collection composition — adding Reactome to a Hallmark run widens
that divisor by ~49% on this data, silently restating previously
published scores.

{normalization}
"""


# ===========================================================================
# Pathology reports (shared; shape sentence differs per view)
# ===========================================================================


def _pathology_section(*, consolidated: bool) -> str:
    """Section describing the scanned pathology reports.

    Worth its own section rather than a row in the source table: it is the
    only modality here that ships a source document rather than parsed
    values, and the decision not to extract text is one consumers need to
    know about before they reach for a `text` column that does not exist.
    """
    if consolidated:
        shape = (
            "Each patient row carries a `samples_pathology_report` array; "
            "`pdf_bytes` holds the document."
        )
    else:
        shape = (
            "The `pathology_report` table has one row per report, with the "
            "document in `pdf_bytes`."
        )
    return f"""\
## Pathology reports

Scanned surgical pathology reports as GDC serves them — **11,208 reports
covering 11,121 cases across 32 projects**. {shape}

TCGA-LAML has none, which is expected rather than missing: acute myeloid
leukaemia has no surgical resection specimen to report on.

### The bytes, not a text extraction

The PDFs are carried **verbatim, with no text extraction applied**. Any
parse is specific to the tool and version that produced it, so extracting
at publication time would freeze one tool's output into the dataset and
lose the original. Shipping the source document means a better parser can
be run later without re-downloading from GDC, and a canonical parse — if
one is added — becomes an additional clearly-labelled column rather than
a replacement.

Practical notes for anyone parsing them:

- These are page scans. Most carry an OCR text layer added upstream of
  GDC, so a pure-Python extractor returns several hundred to a few
  thousand characters for nearly every report — but that layer transcribes
  the barcode strip and handwritten margin notes as noise, and its
  fidelity varies by submitting institution.
- Patient identifiers are redacted out of the page image by GDC before
  distribution.

### Joining to a sample

Every report links to the sample it describes. The GDC file name is
`<case_submitter_id>.<REPORT_UUID>.PDF`, and that UUID is the same value
GDC reports on `sample.pathology_report_uuid` — a key this dataset has
always carried, so reports join to samples without anything new being
invented. Where GDC names the sample directly in the file's
`associated_entities`, that is preferred, with the file-name UUID as
fallback.
"""


def _copy_number_section(*, consolidated: bool) -> str:
    """Section on the two copy-number segment modalities.

    Needs its own section because of one trap: three callers ship for
    overlapping aliquots and a query that ignores `workflow_type` silently
    pools them.
    """
    if consolidated:
        where = (
            "in two array columns on each patient row that answer different "
            "questions. Each record covers one assay run, with the per-segment "
            "values as index-aligned arrays inside it."
        )
        names = (
            "| Column | Caller | Measurement | Files |\n"
            "|---|---|---|---|\n"
            "| `samples_allele_specific_copy_number_segment` | ASCAT2, ASCAT3, AscatNGS | **Integer** total copy number plus its split into `major_copy_number` / `minor_copy_number` | 23,225 |\n"  # noqa: E501
            "| `samples_masked_copy_number_segment` | DNAcopy | **Relative** log2(sample / reference) in `segment_mean`, germline CNVs masked out | 22,629 |"  # noqa: E501
        )
    else:
        where = "in two tables that answer different questions."
        names = (
            "| Table | Caller | Measurement | Files |\n"
            "|---|---|---|---|\n"
            "| `allele_specific_copy_number_segment` | ASCAT2, ASCAT3, AscatNGS | **Integer** total copy number plus its split into `major_copy_number` / `minor_copy_number` | 23,225 |\n"  # noqa: E501
            "| `masked_copy_number_segment` | DNAcopy | **Relative** log2(sample / reference) in `segment_mean`, germline CNVs masked out | 22,629 |"  # noqa: E501
        )
    return f"""\
## Copy number

Copy number ships at **segment level, exactly as GDC serves it**, {where}

{names}

`copy_number = major_copy_number + minor_copy_number` holds everywhere.
`minor_copy_number = 0` with `major_copy_number > 0` is loss of
heterozygosity.

### Filter on `workflow_type`

All three allele-specific callers ship for overlapping aliquots, and each
fits tumour purity and ploidy independently, **so they can disagree**. On
TCGA-CHOL, ASCAT2 and ASCAT3 give the same length-weighted modal copy
number for 33 of 36 shared aliquots — but where they differ they differ
substantially (one aliquot is modal 2 under ASCAT2 and modal 4 under
ASCAT3), and ASCAT3 segments far more coarsely (2,469 segments against
ASCAT2's 6,580 over the same aliquots). AscatNGS is the WGS-based caller;
the other two run on genotyping arrays, recorded in
`experimental_strategy`.

A query that does not filter on `workflow_type` is pooling three different
answers to the same question. ASCAT3 is GDC's current standard.

### The two views are not interchangeable

Nesting each masked segment inside its containing ASCAT3 segment on
TCGA-CHOL (2,590 pairs) gives Spearman **+0.56**, with median
`segment_mean` rising monotonically across integer copy number:

| `copy_number` | 0 | 1 | 2 | 4 | 8 |
|---|---|---|---|---|---|
| median `segment_mean` | −1.93 | −0.44 | +0.07 | +0.22 | +1.25 |

The correlation is only moderate, and that is a property of the data
rather than a defect: ASCAT corrects for purity and ploidy while DNAcopy's
ratio is against a diploid reference, so in a hyperdiploid tumour integer
copy number 3 is copy-neutral relative to its own baseline yet still reads
near log2 0. **Use the allele-specific calls for absolute copy number, the
masked segments for reference-relative ratio.**

One formatting difference is carried through from the source rather than
normalized: the allele-specific segments write `chr1`, and the masked
segments write bare `1`.

### A small tail of over-fragmented masked segments

Most masked files hold 60-100 segments (median 77 in TCGA-LAML, 91 in
TCGA-BRCA, 67 in TCGA-CHOL). A handful hold tens of thousands: **32 of
22,629 files (0.14%) exceed 200 KB**, the largest carrying 50,780 segments
against TCGA-BRCA's per-file maximum of 1,029. They cluster in TCGA-LAML
(12), TCGA-BLCA (9) and TCGA-BRCA (7), and the most extreme are all `-11A-`
matched normals.

This is the signature of a noisy genotyping array, where circular binary
segmentation fails to merge and emits many tiny spurious calls. It is
genuine GDC content and is shipped unmodified, but it is a real trap: an
unfiltered query over this table gets a few samples contributing tens of
thousands of junk rows each, enough to skew any per-segment aggregate.
**`num_probes` is the filter** — the spurious segments are supported by
very few probes.

### Gene-level copy number is deliberately absent

GDC also serves `Gene Level Copy Number` — the same calls projected onto
GENCODE v36 — at roughly 34 GB per workflow. It is not shipped here
because it is exactly reproducible from the allele-specific segments
rather than being independent evidence. (Verified against GDC's own files
on TCGA-CHOL: projecting segments onto the gene model reproduced every
gene call with zero mismatches across three aliquots, and for genes
straddling a segment boundary GDC's `min_copy_number` / `max_copy_number`
are the min and max over the overlapping segments.) It may be added later
as a clearly-labelled derived table.
"""


def _mirna_and_protein_section(*, consolidated: bool) -> str:
    """Section on the two remaining per-specimen molecular assays."""
    mirna_name = (
        "`samples_mirna_expression_quantification`"
        if consolidated
        else "`mirna_expression_quantification`"
    )
    rppa_name = (
        "`samples_protein_expression_quantification`"
        if consolidated
        else "`protein_expression_quantification`"
    )
    mirna_shape = (
        "One record per aliquot, holding ~1,881 miRBase v21 mature miRNAs as "
        "index-aligned arrays"
        if consolidated
        else "One row per (aliquot, miRBase v21 mature miRNA) — ~1,881 miRNAs per aliquot"
    )
    rppa_shape = (
        "One record per portion, holding ~487 antibodies as index-aligned arrays"
        if consolidated
        else "One row per (portion, antibody)"
    )
    return f"""\
## miRNA-Seq and protein expression (RPPA)

### {mirna_name}

{mirna_shape}, from 11,441 files across TCGA. `read_count` is raw;
`reads_per_million_mirna_mapped` is normalized within the aliquot and sums
to exactly 1,000,000 per aliquot.

`cross_mapped` is `Y` when reads for that miRNA also aligned elsewhere in
the genome, so its count is not uniquely attributable. GDC ships the flag
rather than dropping the row and so do we; filter it out if you need clean
attribution. The source column is spelled `cross-mapped` — renamed here
only because the hyphen is not a legal bare SQL identifier.

Isoform-level quantification (`Isoform Expression Quantification`, ~4 GB)
is not shipped.

### {rppa_name}

Reverse Phase Protein Array. {rppa_shape}. 7,906 files
covering **7,827 of 11,428 TCGA cases — the narrowest coverage of any
modality here**, because RPPA was only run on a subset.

Three things to know before using it:

- It is the only modality that attaches to a **portion**, not an aliquot,
  so it carries `portion_id` and resolves `sample_id` through the portion.
- The antibody panel grew over the project's life, and `set_id` records
  which version a measurement came from. A `peptide_target` absent for a
  sample may mean "not on that panel" rather than "measured as zero".
- `protein_expression` is **null where the source says `NA`** — a failed
  or missing measurement, not a zero. On TCGA-CHOL that is 930 of 14,370
  cells (6.5%).

Values are replicate-based normalized log2 signal, centred near 0, and the
sign is meaningful. Agreement with matched RNA is modest and positive, as
expected for protein-vs-transcript: median Spearman +0.26 across shared
targets on TCGA-CHOL.
"""


def _biospecimen_section(*, consolidated: bool) -> str:
    """Section on the BCR biospecimen biotabs."""
    if consolidated:
        intro = (
            "The counterpart to the `clinical_supplement` struct: where that "
            "describes the patient, the `biospecimen_supplement` struct describes "
            "the **specimen chain**"
        )
        shape = (
            "one list-valued sub-field per form. Every sub-field is a list — "
            "unlike `clinical_supplement`, which has a single-dict `patient` slot "
            "— because each of these forms is one-row-per-entity."
        )
        table_col = "Sub-field"
        prefix = ""
    else:
        intro = (
            "The counterpart to `clinical_supplement_*`: where those describe the "
            "patient, `biospecimen_supplement_*` describes the **specimen chain**"
        )
        shape = "one table per form."
        table_col = "Table"
        prefix = "biospecimen_supplement_"
    return f"""\
## Biospecimen supplements

{intro} — how
a tumour got from the operating room to a sequencer, and the pathologist's
read on each slide along the way. 340 BCR biotab files across TCGA
(~76 MB) covering 11,315 cases, {shape}

Some of it restates what the case structure already nests (sample /
portion / analyte / aliquot ids and types). The forms worth reaching for
are the ones with no `/cases` equivalent:

| {table_col} | What is in it |
|---|---|
| `{prefix}slide` | Per-slide `percent_tumor_nuclei`, `percent_necrosis`, `percent_stromal_cells`, `percent_lymphocyte_infiltration`, `section_location` — the QC layer behind "is this sample actually tumour?", and the standard covariate for purity and deconvolution work |
| `{prefix}analyte` | `a260_a280_ratio`, `concentration`, extraction method — nucleic-acid quality, which drives batch effects |
| `{prefix}protocol`, `{prefix}shipment_portion` | Plate, shipment and centre each specimen moved through — the raw material for batch-effect analysis |
| `{prefix}ssf_tumor_samples`, `{prefix}ssf_normal_controls` | Site-specific factors: the disease-specific pathology fields the pan-cancer clinical schema has no column for |
| `{prefix}cqcf` | The submitting centre's clinical quality control form (TCGA-LUAD only) |

Like the clinical supplements these are **flex-schema**: the column set
differs by project and by submitting centre, so the shape is inferred per
project rather than padded into a pan-cancer union. Forms with no data for
a project are omitted entirely (only TCGA-LUAD has `cqcf`; only 9 projects
have `auxiliary`). All fields are typed as strings with BCR sentinels like
`[Not Available]` preserved verbatim.

Records are keyed to the patient by BCR barcode. The specimen-level forms
are keyed on their own entity and several omit the patient barcode column
entirely, in which case it is recovered as the first three groups of the
entity barcode (`TCGA-3X-AAV9-01A-11D-A42S-01` → `TCGA-3X-AAV9`) — a
property of the TCGA barcode grammar, not a heuristic.

Two submitters ship these files: `nationwidechildrens.org` for 334 of the
340, and `genome.wustl.edu` for 6 (all TCGA-LUAD). Where both ship the
same form for one project their records are concatenated, and the parquet
schema is the union of their columns.
"""


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
            _pathology_section(consolidated=True),
            _copy_number_section(consolidated=True),
            _mirna_and_protein_section(consolidated=True),
            _biospecimen_section(consolidated=True),
            _ssgsea_section(consolidated=True),
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
            _pathology_section(consolidated=False),
            _copy_number_section(consolidated=False),
            _mirna_and_protein_section(consolidated=False),
            _biospecimen_section(consolidated=False),
            _ssgsea_section(consolidated=False),
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
