from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

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

t = pq.read_table("TCGA-CHOL/train.parquet")
patients = [TcgaHfPatient.model_validate(r) for r in t.to_pylist()]
```
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
        lines.append(f"        path: {project}/*.parquet")
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
