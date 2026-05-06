"""Row emitters for the tabular HF dataset (`tcga2hf-tabular-open`).

Companion to `tcga2hf.clinical` (which produces the consolidated patient
row with both clinical and molecular vectors nested). The tabular layout
splits each project into four tables, named after stable GDC concepts:

  - `cases`                            GDC `cases.json` row, nested
                                       structure preserved (one row per
                                       patient)
  - `masked_somatic_mutation`          one row per MAF variant — snake-
                                       cased GDC data_type for somatic
                                       calls (today: ensemble-masked
                                       WXS-derived MAFs)
  - `gene_expression_quantification`   one row per (aliquot, gene) —
                                       snake-cased GDC data_type for
                                       RNA-Seq expression (today: STAR
                                       counts TSVs)
  - `files`                            one row per (file, case) from the
                                       per-modality manifests

The `cases` table reuses `clinical.to_patient_rows` and drops the two
molecular-vector columns. The molecular and provenance tables read the
same raw inputs the consolidated build does.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from tcga2hf.schema import TABULAR_CASES_FIELDS, TABULAR_TABLES
from tcga2hf_pipeline import clinical as _clinical_mod
from tcga2hf_pipeline import expression as _expression_mod
from tcga2hf_pipeline import mutations as _mutations_mod

# Names of `cases` table columns, used to re-shape patient rows. Excludes
# the two molecular-vector columns that the consolidated row carries (they
# have their own flat tables in the tabular layout).
_CASES_COLS: list[str] = [f.name for f in TABULAR_CASES_FIELDS]


# ---------------------------------------------------------------------------
# Cases (one nested row per patient — mirrors `cases.json` shape)
# ---------------------------------------------------------------------------


def _cases_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Build one nested row per case for the `cases` tabular table.

    Reuses `clinical.to_patient_rows` (which already mirrors the GDC
    `cases.json` structure with demographic / diagnoses / follow_ups /
    exposures / family_histories / samples nested), then drops the two
    molecular-vector columns those rows carry by default.
    """
    rows = _clinical_mod.to_patient_rows(cases)
    return [{k: row[k] for k in _CASES_COLS} for row in rows]


# ---------------------------------------------------------------------------
# Molecular emitters (mutations + expression)
# ---------------------------------------------------------------------------


def _aliquot_to_sample(case: dict[str, Any]) -> dict[str, str]:
    """Map aliquot_id → sample_id for one case (walks samples → … → aliquots)."""
    out: dict[str, str] = {}
    for s in case.get("samples") or []:
        sid = s.get("sample_id")
        if not sid:
            continue
        for p in s.get("portions") or []:
            for an in p.get("analytes") or []:
                for aq in an.get("aliquots") or []:
                    aid = aq.get("aliquot_id")
                    if aid:
                        out[aid] = sid
    return out


def _aliquot_to_submitter(case: dict[str, Any]) -> dict[str, str]:
    """Map aliquot_id → aliquot.submitter_id (for browseability on the
    expression table, where we don't otherwise carry the submitter id)."""
    out: dict[str, str] = {}
    for s in case.get("samples") or []:
        for p in s.get("portions") or []:
            for an in p.get("analytes") or []:
                for aq in an.get("aliquots") or []:
                    aid = aq.get("aliquot_id")
                    if aid:
                        out[aid] = aq.get("submitter_id")
    return out


def _mutations_rows(
    cases: list[dict[str, Any]],
    by_case: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """One row per MAF variant, with FKs resolved against the case's samples."""
    rows: list[dict[str, Any]] = []
    case_index = {c.get("case_id"): c for c in cases}
    for case_id, variants in by_case.items():
        case = case_index.get(case_id)
        if case is None:
            continue
        case_submitter_id = case.get("submitter_id")
        a2s = _aliquot_to_sample(case)
        for v in variants:
            row = dict(v)
            row["case_submitter_id"] = case_submitter_id
            row["tumor_sample_id"] = a2s.get(v.get("Tumor_Sample_UUID"))
            row["matched_normal_sample_id"] = a2s.get(v.get("Matched_Norm_Sample_UUID"))
            rows.append(row)
    return rows


def _expression_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per gene in each aliquot's expression TSV, with FKs prepended.

    Concatenates the per-gene rows from every aliquot's TSV under the
    project. The 4 N_* alignment-summary rows that sit at the top of each
    TSV (gene_id values `N_unmapped`, `N_multimapping`, `N_noFeature`,
    `N_ambiguous` — they have empty gene_name / gene_type and aren't real
    genes) are filtered out: they occupy the same column shape but have
    different semantics, and including them complicates SQL filtering.
    """
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return []

    case_index = {c.get("case_id"): c for c in cases}
    manifest = json.loads(manifest_path.read_text())
    rows: list[dict[str, Any]] = []
    for entry in manifest:
        file_path = expr_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id, aliquot_id = _expression_mod._file_aliquot_and_case(entry)
        if not case_id or not aliquot_id:
            continue
        case = case_index.get(case_id)
        if case is None:
            continue
        df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)
        # Drop the 4 alignment-summary rows; keep only true gene records.
        df = df[~df["gene_id"].str.startswith("N_", na=False)].reset_index(drop=True)
        a2sub = _aliquot_to_submitter(case)
        case_submitter_id = case.get("submitter_id")
        aliquot_submitter_id = a2sub.get(aliquot_id)
        source_file_id = entry["file_id"]
        for r in df.itertuples(index=False):
            row = {
                "case_id": case_id,
                "case_submitter_id": case_submitter_id,
                "aliquot_id": aliquot_id,
                "aliquot_submitter_id": aliquot_submitter_id,
                "source_file_id": source_file_id,
                "gene_id": _na_to_none(r.gene_id),
                "gene_name": _na_to_none(r.gene_name),
                "gene_type": _na_to_none(r.gene_type),
                "unstranded": _to_int_or_none(r.unstranded),
                "stranded_first": _to_int_or_none(r.stranded_first),
                "stranded_second": _to_int_or_none(r.stranded_second),
                "tpm_unstranded": _to_float_or_none(r.tpm_unstranded),
                "fpkm_unstranded": _to_float_or_none(r.fpkm_unstranded),
                "fpkm_uq_unstranded": _to_float_or_none(r.fpkm_uq_unstranded),
            }
            rows.append(row)
    return rows


def _na_to_none(v: Any) -> Any:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return v


def _to_int_or_none(v: Any) -> int | None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return int(v)


def _to_float_or_none(v: Any) -> float | None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return None
    return float(v)


# ---------------------------------------------------------------------------
# Files (provenance from per-modality manifests)
# ---------------------------------------------------------------------------


def _files_rows(project_raw_dir: Path) -> list[dict[str, Any]]:
    """One row per (file, case) in the modality manifests under `<project>/`.

    Concatenates every modality directory's `manifest.json` and explodes
    each entry's `cases` list — most files reference exactly one case, but
    we explode for safety (a future modality might attach multi-case files).
    """
    rows: list[dict[str, Any]] = []
    if not project_raw_dir.exists():
        return rows
    for modality_dir in sorted(p for p in project_raw_dir.iterdir() if p.is_dir()):
        manifest_path = modality_dir / "manifest.json"
        if not manifest_path.exists():
            continue
        manifest = json.loads(manifest_path.read_text())
        modality = modality_dir.name
        for entry in manifest:
            common = {
                "file_id": entry.get("file_id"),
                "file_name": entry.get("file_name"),
                "file_size": entry.get("file_size"),
                "md5sum": entry.get("md5sum"),
                "data_category": entry.get("data_category"),
                "data_type": entry.get("data_type"),
                "data_format": entry.get("data_format"),
                "experimental_strategy": entry.get("experimental_strategy"),
                "workflow_type": entry.get("workflow_type"),
                "access": entry.get("access"),
                "modality": modality,
            }
            cases = entry.get("cases") or []
            if not cases:
                rows.append({"case_id": None, "case_submitter_id": None, **common})
                continue
            for case in cases:
                rows.append(
                    {
                        "case_id": case.get("case_id"),
                        "case_submitter_id": case.get("submitter_id"),
                        **common,
                    }
                )
    return rows


# ---------------------------------------------------------------------------
# Top-level builder
# ---------------------------------------------------------------------------


def build_tables(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> dict[str, list[dict[str, Any]]]:
    """Build every tabular table for one project. Keys match `TABULAR_TABLES`."""
    mut_by_case = _mutations_mod.load_for_project(project_raw_dir)

    return {
        "cases": _cases_rows(cases),
        "masked_somatic_mutation": _mutations_rows(cases, mut_by_case),
        "gene_expression_quantification": _expression_rows(cases, project_raw_dir),
        "files": _files_rows(project_raw_dir),
    }


# ---------------------------------------------------------------------------
# Parquet writer
# ---------------------------------------------------------------------------

# Smaller row groups for the "expression" table since each row is narrow
# (~12 cols of small scalars) but the count is huge — 60660 rows per
# aliquot × 600 aliquots ≈ 36M rows for LUAD. 100k rows per group keeps
# each group well under the HF Data Studio 300 MB scan limit and gives
# the SQL viewer fast random access. The `cases` table's rows are wide
# (full nested case JSON); 50 rows per group there matches the
# consolidated build's sizing.
_ROW_GROUP_SIZE_DEFAULT = 50
_ROW_GROUP_SIZE_BY_TABLE: dict[str, int] = {
    "gene_expression_quantification": 100_000,
}


def write_tables(
    tables: dict[str, list[dict[str, Any]]],
    processed_dir: Path,
    project_id: str,
) -> dict[str, Path]:
    """Write each table to `<processed_dir>/<project_id>/<table>/data.parquet`.

    Uses the same row-group + page-index settings as the consolidated build
    so HF Data Studio can scan large tables (notably `expression`) without
    hitting its scan-size limit.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    out_paths: dict[str, Path] = {}
    project_dir = processed_dir / project_id
    project_dir.mkdir(parents=True, exist_ok=True)
    for table_name, rows in tables.items():
        schema = TABULAR_TABLES[table_name]
        out_path = project_dir / table_name / "data.parquet"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        table = pa.Table.from_pylist(rows, schema=schema)
        pq.write_table(
            table,
            out_path,
            row_group_size=_ROW_GROUP_SIZE_BY_TABLE.get(table_name, _ROW_GROUP_SIZE_DEFAULT),
            write_page_index=True,
        )
        out_paths[table_name] = out_path
    return out_paths
