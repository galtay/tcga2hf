from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path


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

    body = f"""
# TCGA Patients (Open Access)

**Open-access** TCGA patient data from the NCI Genomic Data Commons (GDC).
One HuggingFace subset per TCGA project; one row per patient.

- **Projects included:** {projects_md}
- **Generated:** {timestamp}
- **Source:** NCI GDC `/cases` endpoint, open-access tier only.
- **Schema:** derived from the [GDC Data Dictionary][gdc-dict]; the live
  dictionary snapshot the GDC was serving when the data was fetched is
  hashed into each project's `gdc_status.json` (and stored alongside the raw
  data on the producer side, not shipped with the parquet).
{release_md}

## GDC data model primer

A short orientation for downstream users. The GDC organizes everything around a
**case** (a single patient) with a hierarchy of physical biospecimens beneath
it, plus assay-derived data files attached to those biospecimens.

**Biospecimen hierarchy** (per GDC's [data dictionary][gdc-dict]):

```
case          one patient (TCGA-XX-1234)
└── sample    physical specimen taken from the patient at one timepoint
              (Primary Tumor, Solid Tissue Normal, Blood Derived Normal, ...)
    └── portion    a piece of that sample for a specific lab process
        └── analyte    extracted material of one type (DNA or RNA)
            └── aliquot    a vial of that analyte handed off for sequencing
```

We preserve this hierarchy verbatim — each `samples[i]` carries
`portions[j].analytes[k].aliquots[m]` exactly as the GDC returns it. No
flattening, no field hoisting; portion- and analyte-level fields like
`is_ffpe`, `a260_a280_ratio`, and `normal_tumor_genotype_snp_match` are
preserved. The `TcgaHfPatient` reference implementation (below) provides
flat-aliquot convenience views over this tree.

**Top-level molecular columns** are named `samples_<gdc_data_type_snake_case>`
(e.g. `samples_masked_somatic_mutation`, `samples_gene_expression_quantification`).
Each entry carries FK fields back to the patient's `samples[].portions[].analytes[].aliquots[]`
so cross-modality joins are local to the row.

**Timeline anchor is uniform per the GDC dictionary.** Every `days_to_*` field
is documented as days from the case's `index_date` (a top-level field on each
row; for TCGA usually `"Diagnosis"`). `TcgaHfPatient.timeline()` returns every
dated event for the patient on this single anchor — clinical (consent →
diagnosis → treatments → follow-ups → lost-to-follow-up → death) plus
biospecimen (`sample_procurement` from `days_to_sample_procurement`,
`bcr_receipt` from `days_to_collection`).

For some TCGA cases `days_to_collection` exceeds `days_to_death`. We don't
attempt to reinterpret these — we surface the count via
`TcgaHfPatient.consistency_check()` as `bcr_receipts_after_death` and leave
interpretation to the consumer.

**Reference Python implementation:** the `tcga2hf` package on GitHub ships a
fully-typed pydantic `TcgaHfPatient` model that mirrors this schema and adds
convenience joins (tumor/normal pairs, mutations-by-gene, expression-by-gene,
longitudinal timeline). Useful both as a loader and as documentation for what
the data means.

```python
import pyarrow.parquet as pq
from tcga2hf.models import TcgaHfPatient

t = pq.read_table("TCGA-CHOL/train.parquet")
patients = [TcgaHfPatient.model_validate(r) for r in t.to_pylist()]
```

**Authoritative GDC references:**
- [Data dictionary][gdc-dict] (every entity + field definition)
- [Biospecimen Encyclopedia](https://docs.gdc.cancer.gov/Encyclopedia/pages/Biospecimen/)
- [MAF format spec](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/)
- [Gene Expression Quantification spec](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
- [Sample Type codes](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes)
- [TCGA Barcode reference](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)

[gdc-dict]: https://docs.gdc.cancer.gov/Data_Dictionary/

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

References:
[GDC Policies](https://gdc.cancer.gov/about-gdc/gdc-policies),
[GDC Encyclopedia — Controlled Access][controlled] (defines what is *not* here),
[NIH Genomic Data Sharing Policy](https://sharing.nih.gov/genomic-data-sharing).

[controlled]: https://docs.gdc.cancer.gov/Encyclopedia/pages/Controlled_Access/

## Disclaimer

Prototype dataset. Schema, included projects, and column coverage are still
evolving. Re-derive from the GDC for any analysis where freshness matters.
"""

    out_path = processed_dir / "README.md"
    processed_dir.mkdir(parents=True, exist_ok=True)
    out_path.write_text(frontmatter + body)
    return out_path
