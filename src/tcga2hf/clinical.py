from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

from tcga2hf.gdc import GDCClient, in_
from tcga2hf.schema import (
    ALIQUOT_FIELDS,
    ANALYTE_FIELDS,
    DEMOGRAPHIC_FIELDS,
    DIAGNOSIS_FIELDS,
    EXPOSURE_FIELDS,
    FAMILY_HISTORY_FIELDS,
    FOLLOW_UP_FIELDS,
    PATIENTS,
    PORTION_FIELDS,
    SAMPLE_FIELDS,
    TREATMENT_FIELDS,
)

EXPANSIONS = [
    "demographic",
    "diagnoses",
    "diagnoses.treatments",
    "follow_ups",
    "exposures",
    "family_histories",
    "project",
    # GDC expand is per-level: every intermediate must be listed for its own fields
    # to come back. Without this, e.g. analyte.analyte_type would be missing.
    "samples",
    "samples.portions",
    "samples.portions.analytes",
    "samples.portions.analytes.aliquots",
]

# `fields` is required by the API but does not restrict expanded entities.
# Case-level scalar fields (per GDC dictionary): index_date names the timeline
# anchor (TCGA: "Diagnosis"); the days_to_* and consent/follow-up status all
# share that anchor. Have to request each explicitly — they don't come back
# with `expand`.
TOP_LEVEL_FIELDS = [
    "case_id",
    "submitter_id",
    "primary_site",
    "disease_type",
    "index_date",
    "consent_type",
    "days_to_consent",
    "days_to_lost_to_followup",
    "lost_to_followup",
]


def fetch_clinical(projects: list[str], client: GDCClient) -> list[dict[str, Any]]:
    return client.cases(
        filters=in_("project.project_id", projects),
        fields=TOP_LEVEL_FIELDS,
        expand=EXPANSIONS,
        page_size=200,
    )


def _pick(obj: dict[str, Any] | None, fields: list[pa.Field]) -> dict[str, Any] | None:
    """Return a dict containing exactly the named fields (None for missing).

    Returns None if the source object itself is missing — pyarrow renders that as a
    null struct, distinct from a struct of all-null fields.
    """
    if obj is None:
        return None
    return {f.name: obj.get(f.name) for f in fields}


def _sort_temporal(
    items: list[dict[str, Any]], days_field: str, id_field: str
) -> list[dict[str, Any]]:
    """Sort by (days_field asc, id_field asc) with null `days_field` last.

    The `is None` sentinel keeps None values from being compared to ints — Python
    short-circuits tuple comparison on the first differing element, so a True/False
    bool prefix safely segregates the nullable group.
    """
    return sorted(
        items,
        key=lambda x: (x.get(days_field) is None, x.get(days_field), x.get(id_field) or ""),
    )


def _sort_by_id(items: list[dict[str, Any]], id_field: str) -> list[dict[str, Any]]:
    return sorted(items, key=lambda x: x.get(id_field) or "")


def _diagnosis_dict(dx: dict[str, Any]) -> dict[str, Any]:
    out = {f.name: dx.get(f.name) for f in DIAGNOSIS_FIELDS if f.name != "treatments"}
    treatments = [_pick(tx, TREATMENT_FIELDS) for tx in (dx.get("treatments") or [])]
    out["treatments"] = _sort_temporal(treatments, "days_to_treatment_start", "treatment_id")
    return out


def _aliquot_dict(aliquot: dict[str, Any]) -> dict[str, Any]:
    return {f.name: aliquot.get(f.name) for f in ALIQUOT_FIELDS}


def _analyte_dict(analyte: dict[str, Any]) -> dict[str, Any]:
    out = {f.name: analyte.get(f.name) for f in ANALYTE_FIELDS if f.name != "aliquots"}
    aliquots = [_aliquot_dict(a) for a in (analyte.get("aliquots") or [])]
    out["aliquots"] = _sort_by_id(aliquots, "aliquot_id")
    return out


def _portion_dict(portion: dict[str, Any]) -> dict[str, Any]:
    out = {f.name: portion.get(f.name) for f in PORTION_FIELDS if f.name != "analytes"}
    analytes = [_analyte_dict(a) for a in (portion.get("analytes") or [])]
    out["analytes"] = _sort_by_id(analytes, "analyte_id")
    return out


def _sample_dict(sample: dict[str, Any]) -> dict[str, Any]:
    out = {f.name: sample.get(f.name) for f in SAMPLE_FIELDS if f.name != "portions"}
    portions = [_portion_dict(p) for p in (sample.get("portions") or [])]
    out["portions"] = _sort_by_id(portions, "portion_id")
    return out


def _patient_row(case: dict[str, Any]) -> dict[str, Any]:
    project = case.get("project") or {}
    diagnoses = [_diagnosis_dict(dx) for dx in (case.get("diagnoses") or [])]
    follow_ups = [_pick(fu, FOLLOW_UP_FIELDS) for fu in (case.get("follow_ups") or [])]
    exposures = [_pick(ex, EXPOSURE_FIELDS) for ex in (case.get("exposures") or [])]
    family_histories = [
        _pick(fh, FAMILY_HISTORY_FIELDS) for fh in (case.get("family_histories") or [])
    ]
    samples = [_sample_dict(s) for s in (case.get("samples") or [])]
    case_id = case.get("case_id")
    return {
        "case_id": case_id,
        "case_submitter_id": case.get("submitter_id"),
        "project_id": project.get("project_id"),
        "gdc_portal_url": (
            f"https://portal.gdc.cancer.gov/cases/{case_id}" if case_id else None
        ),
        "primary_site": case.get("primary_site"),
        "disease_type": case.get("disease_type"),
        "index_date": case.get("index_date"),
        "consent_type": case.get("consent_type"),
        "days_to_consent": case.get("days_to_consent"),
        "days_to_lost_to_followup": case.get("days_to_lost_to_followup"),
        "lost_to_followup": case.get("lost_to_followup"),
        "demographic": _pick(case.get("demographic"), DEMOGRAPHIC_FIELDS),
        "diagnoses": _sort_temporal(diagnoses, "days_to_diagnosis", "diagnosis_id"),
        "follow_ups": _sort_temporal(follow_ups, "days_to_follow_up", "follow_up_id"),
        "exposures": _sort_by_id(exposures, "exposure_id"),
        "family_histories": _sort_by_id(family_histories, "family_history_id"),
        "samples": _sort_temporal(samples, "days_to_collection", "sample_id"),
        # Molecular columns default to []; build step fills them in if data is fetched.
        "samples_masked_somatic_mutation": [],
        "samples_gene_expression_quantification": [],
    }


def to_patient_rows(cases: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return [_patient_row(case) for case in cases]


# Row group size targeting the HuggingFace Dataset Viewer's 300 MB scan limit
# (see https://huggingface.co/docs/hub/datasets-data-files-configuration). Each
# patient row is ~1-1.5 MB compressed, dominated by the 60k-gene expression
# arrays, so 50 rows per group ≈ 50-75 MB — comfortably below the limit and
# small enough that the Viewer's per-page random access is snappy.
_ROW_GROUP_SIZE = 50


def write_patients(rows: list[dict[str, Any]], processed_dir: Path, project_id: str) -> Path:
    """Write a single project's patient rows to <project_id>/train.parquet.

    The per-project directory layout matches the HuggingFace `configs:` convention
    so each TCGA project surfaces as its own loadable subset.

    Row groups are sized for the HF Dataset Viewer (see `_ROW_GROUP_SIZE`) and
    a page index is written so the Viewer can read only the bytes it needs to
    render a page of rows rather than scanning the whole row group.
    """
    out_path = processed_dir / project_id / "train.parquet"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist(rows, schema=PATIENTS)
    pq.write_table(
        table,
        out_path,
        row_group_size=_ROW_GROUP_SIZE,
        write_page_index=True,
    )
    return out_path
