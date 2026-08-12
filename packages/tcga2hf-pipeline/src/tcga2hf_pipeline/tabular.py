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
  - `pathology_report`                 one row per scanned report PDF —
                                       the document bytes verbatim, no
                                       text extraction (see
                                       `tcga2hf_pipeline.pathology`)
  - `ssgsea_scores_<collection>`       one row per (aliquot, pathway) of
                                       raw ssGSEA pathway activity, one
                                       table per MSigDB collection
  - `ssgsea_stats_<collection>`        reference distributions for those
                                       scores (this project + pan-cancer),
                                       written after every project scores
  - `files`                            one row per (file, case) from the
                                       per-modality manifests

Two further tables are built elsewhere and land in the same tree:
`survival_derived` (projected off the cases rows by `derived_survival_rows`
after the caller attaches survival) and the flex-schema
`clinical_supplement_*` set (one per BCR biotab form).

The `cases` table reuses `clinical.to_patient_rows` and drops the
per-file modality columns. The molecular, document, and provenance tables
read the same raw inputs the consolidated build does.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from tcga2hf import schema as _schema_mod
from tcga2hf.schema import SSGSEA_COLLECTIONS, TABULAR_CASES_FIELDS, TABULAR_TABLES

from tcga2hf_pipeline import clinical as _clinical_mod
from tcga2hf_pipeline import clinical_supplement as _clinical_supplement_mod
from tcga2hf_pipeline import expression as _expression_mod
from tcga2hf_pipeline import msigdb as _msigdb_mod
from tcga2hf_pipeline import mutations as _mutations_mod
from tcga2hf_pipeline import pathology as _pathology_mod
from tcga2hf_pipeline import ssgsea as _ssgsea_mod

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


def _pathology_report_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
) -> list[dict[str, Any]]:
    """One row per pathology report PDF, with case FKs prepended.

    Reuses the consolidated loader (which reads the bytes and resolves the
    sample FK from `associated_entities`), then applies the same filename-UUID
    fallback `pathology.attach` does — here against the raw GDC case dict
    rather than a built patient row.
    """
    by_case = _pathology_mod.load_for_project(project_raw_dir)
    if not by_case:
        return []

    rows: list[dict[str, Any]] = []
    case_index = {c.get("case_id"): c for c in cases}
    for case_id, reports in by_case.items():
        case = case_index.get(case_id)
        if case is None:
            continue
        uuid_to_sample = {
            s["pathology_report_uuid"]: s
            for s in (case.get("samples") or [])
            if s.get("pathology_report_uuid")
        }
        for r in reports:
            row = dict(r)
            if row["sample_id"] is None:
                fallback = uuid_to_sample.get(row["pathology_report_uuid"] or "")
                if fallback:
                    row["sample_id"] = fallback.get("sample_id")
                    row["sample_submitter_id"] = fallback.get("submitter_id")
            row["case_id"] = case_id
            row["case_submitter_id"] = case.get("submitter_id")
            rows.append(row)
    rows.sort(key=lambda r: (r.get("case_submitter_id") or "", r.get("file_name") or ""))
    return rows


# ---------------------------------------------------------------------------
# ssGSEA pathway activity
# ---------------------------------------------------------------------------


def _ssgsea_scores_rows(
    cases: list[dict[str, Any]],
    project_raw_dir: Path,
    collection: str,
    msigdb_dir: Path,
    load_matrix: Any = None,
) -> list[dict[str, Any]]:
    """One row per (aliquot, pathway) of raw ssGSEA scores for one project.

    Depends on nothing outside this project: raw scores are a function of a
    sample's own expression, the gene universe and the gene set. That is
    what lets pathway activity be built inside the ordinary per-project
    loop; the cohort-level reference distributions are a separate pass (see
    `ssgsea_stats_rows`).
    """
    gmt_path = msigdb_dir / _msigdb_mod.COLLECTIONS[collection].file_name
    if not gmt_path.exists():
        return []
    # Reading every STAR TSV dominates the cost, and the matrix is the same
    # for every collection, so the caller passes a memoized loader rather
    # than paying for it once per collection.
    loader = load_matrix or (lambda: _expression_mod.tpm_matrix_for_project(project_raw_dir))
    genes, matrix, records = loader()
    if not records:
        return []

    kept, stats = _ssgsea_mod.map_gene_sets(_ssgsea_mod.load_gmt(gmt_path), genes)
    if not kept:
        return []
    names = list(kept)
    by_name = {s["pathway"]: s for s in stats}
    scores = _ssgsea_mod.ssgsea_scores(matrix.astype(np.float64), [kept[n] for n in names])

    case_index = {c.get("case_id"): c for c in cases}
    rows: list[dict[str, Any]] = []
    for j, rec in enumerate(records):
        case = case_index.get(rec["case_id"])
        if case is None:
            continue
        a2s = _aliquot_to_sample(case)
        a2sub = _aliquot_to_submitter(case)
        sample_id = a2s.get(rec["aliquot_id"])
        sample = next(
            (s for s in (case.get("samples") or []) if s.get("sample_id") == sample_id), {}
        )
        common = {
            "case_id": rec["case_id"],
            "case_submitter_id": case.get("submitter_id"),
            "sample_id": sample_id,
            "sample_submitter_id": sample.get("submitter_id"),
            "sample_type": sample.get("sample_type"),
            "aliquot_id": rec["aliquot_id"],
            "aliquot_submitter_id": a2sub.get(rec["aliquot_id"]),
            "source_file_id": rec["source_file_id"],
        }
        for i, name in enumerate(names):
            rows.append(
                {
                    **common,
                    "pathway": name,
                    "pathway_url": _msigdb_mod.geneset_url(name),
                    "matched_gene_count": int(by_name[name]["matched_gene_count"]),
                    "original_gene_count": int(by_name[name]["original_gene_count"]),
                    "score_raw": float(scores[i, j]),
                }
            )
    return rows


def _describe(series: pd.Series) -> dict[str, Any]:
    """Reference distribution for one (population, pathway) cell."""
    n = int(series.size)
    return {
        "n_aliquots": n,
        "mean": float(series.mean()),
        # ddof=1 is undefined for a single observation; emit null rather
        # than a misleading 0.0.
        "sd": float(series.std(ddof=1)) if n > 1 else None,
        "min": float(series.min()),
        "q25": float(series.quantile(0.25)),
        "median": float(series.median()),
        "q75": float(series.quantile(0.75)),
        "max": float(series.max()),
    }


def ssgsea_stats_rows(
    processed_dir: Path,
    collection: str,
    sample_types: tuple[str, ...] = ("Primary Tumor",),
) -> dict[str, list[dict[str, Any]]]:
    """Reference distributions per project, computed from the written scores.

    Returns `{project_id: rows}`. Every project's rows carry that project's
    own distributions **and** the pan-cancer ones, so a consumer who loads a
    single project can still normalize against all of TCGA without scanning
    every config.

    This reads the `ssgsea_scores_*` parquets rather than rescoring: the
    stats are pure aggregates of published values, and computing them from
    the published bytes is the strongest guarantee that they agree. It also
    means `--table ssgsea_stats_<collection>` can refresh them alone.

    Everything here is composition-dependent — adding a project changes the
    pan-cancer rows — which is precisely why it lives in this small table
    instead of as extra columns on the immutable scores.
    """
    import pyarrow.parquet as pq

    table_name = _schema_mod.ssgsea_scores_table(collection)
    frames: dict[str, pd.DataFrame] = {}
    for path in sorted(processed_dir.glob(f"*/{table_name}/data.parquet")):
        df = pq.read_table(
            path, columns=["pathway", "sample_type", "score_raw"]
        ).to_pandas()
        if len(df):
            frames[path.parent.parent.name] = df
    if not frames:
        return {}

    def _rows_for(df: pd.DataFrame, population: str, project_id: str | None) -> list[dict]:
        out = []
        for sample_type in (None, *sample_types):
            sub = df if sample_type is None else df[df["sample_type"] == sample_type]
            if not len(sub):
                continue
            for pathway, series in sub.groupby("pathway")["score_raw"]:
                out.append(
                    {
                        "population": population,
                        "project_id": project_id,
                        "sample_type": sample_type,
                        "pathway": pathway,
                        **_describe(series),
                    }
                )
        return out

    pan = _rows_for(pd.concat(frames.values(), ignore_index=True), "pan_cancer", None)
    return {proj: _rows_for(df, "project", proj) + pan for proj, df in frames.items()}


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
                # Absent from manifests written before version provenance
                # was recorded; those rows report null rather than a guess.
                "gdc_version": entry.get("gdc_version"),
                "gdc_first_release": entry.get("gdc_first_release"),
                "gdc_superseded": entry.get("gdc_superseded"),
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
    only: set[str] | None = None,
    msigdb_dir: Path | None = None,
) -> dict[str, list[dict[str, Any]]]:
    """Build tabular tables for one project (all of them unless `only` is given).

    Returns one entry per table keyed by table name. The fixed-schema
    keys (`cases`, `masked_somatic_mutation`, `gene_expression_quantification`,
    `pathology_report`, `files`, `survival_derived`) match `TABULAR_TABLES`. The
    `clinical_supplement_*` keys are flex-schema tables — their parquet
    column set is inferred per project, since BCR biotab forms vary by
    cancer type (e.g. BLCA's clinical_patient form has bladder-specific
    fields that don't exist in CHOL's). Cross-project queries union via
    HF `concatenate_datasets` with NULL padding — the same pattern
    cBioPortal and the GDC use for their per-study clinical exports.

    `only` restricts the build to the named tables. Each table's rows are
    produced by a thunk that is called only when its table is wanted, so
    the expensive emitters — `gene_expression_quantification` re-reads
    every STAR TSV in the project, `masked_somatic_mutation` every MAF —
    cost nothing when a caller is appending one new table to an existing
    tree. `cases` is still built whenever `survival_derived` is requested,
    since the derived endpoints are projected off those rows.

    Note: the `survival_derived` table starts empty here. The caller must
    `survival.attach_survival(tables["cases"])` and then call
    `derived_survival_rows(tables["cases"])` to populate it. Splitting
    that off keeps `build_tables` purely deterministic on raw inputs;
    survival re-derivation depends on Clinical Supplement data the caller
    may attach beforehand.
    """

    def _mutation_rows() -> list[dict[str, Any]]:
        return _mutations_rows(cases, _mutations_mod.load_for_project(project_raw_dir))

    thunks: dict[str, Any] = {
        "cases": lambda: _cases_rows(cases),
        "masked_somatic_mutation": _mutation_rows,
        "gene_expression_quantification": lambda: _expression_rows(cases, project_raw_dir),
        "files": lambda: _files_rows(project_raw_dir),
        "survival_derived": list,
        "pathology_report": lambda: _pathology_report_rows(cases, project_raw_dir),
    }
    # ssGSEA scores are per-project pure; the matching stats tables are a
    # cohort-level aggregate written by the caller after every project's
    # scores exist (see `ssgsea_stats_rows`), so they start empty here for
    # the same reason `survival_derived` does.
    msigdb = msigdb_dir if msigdb_dir is not None else project_raw_dir.parent / "msigdb"
    _matrix_cache: list[Any] = []

    def _matrix():
        if not _matrix_cache:
            _matrix_cache.append(_expression_mod.tpm_matrix_for_project(project_raw_dir))
        return _matrix_cache[0]

    for coll in SSGSEA_COLLECTIONS:
        thunks[_schema_mod.ssgsea_scores_table(coll)] = (
            lambda c=coll: _ssgsea_scores_rows(cases, project_raw_dir, c, msigdb, _matrix)
        )
        thunks[_schema_mod.ssgsea_stats_table(coll)] = list

    wanted = set(thunks) if only is None else set(only)
    # survival_derived is projected off the cases rows, so it drags `cases`
    # in even when the caller didn't ask for it directly.
    if "survival_derived" in wanted:
        wanted.add("cases")

    tables: dict[str, list[dict[str, Any]]] = {
        name: thunk() for name, thunk in thunks.items() if name in wanted
    }

    if only is None or any(t.startswith("clinical_supplement_") for t in wanted):
        supp_rows = _clinical_supplement_mod.build_tabular_rows(
            project_raw_dir / "clinical_supplement"
        )
        for suffix, rows in supp_rows.items():
            name = f"clinical_supplement_{suffix}"
            if only is None or name in wanted:
                tables[name] = rows
    return tables


def derived_survival_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Project the `survival_derived` struct off each case row into a flat
    table. Caller must run `survival.attach_survival(cases)` first.

    Output: one row per case with `case_submitter_id` plus the 8 derived
    survival fields (os_event, os_time, dss_event, dss_time, pfi_event,
    pfi_time, dfi_event, dfi_time). Rows whose `survival_derived` is None
    or all-null are skipped — they have nothing to contribute.
    """
    out: list[dict[str, Any]] = []
    for r in cases:
        sd = r.get("survival_derived") or {}
        if all(sd.get(k) is None for k in (
            "os_event", "os_time", "dss_event", "dss_time",
            "pfi_event", "pfi_time", "dfi_event", "dfi_time",
        )):
            continue
        out.append({"case_submitter_id": r["case_submitter_id"], **sd})
    return out


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
    # ssGSEA scores are the same shape as expression -- narrow rows, very
    # many of them (575k per project for Hallmark, ~17M for a Reactome-sized
    # collection) -- so they get the same row-group treatment.
    **{f"ssgsea_scores_{c}": 100_000 for c in SSGSEA_COLLECTIONS},
}


def write_tables(
    tables: dict[str, list[dict[str, Any]]],
    processed_dir: Path,
    project_id: str,
) -> dict[str, Path]:
    """Write each table to `<processed_dir>/<project_id>/<table>/data.parquet`.

    Tables in `TABULAR_TABLES` use their fixed pan-cancer schema. Tables
    not in that map (today: `clinical_supplement_*`) get schema inferred
    from rows — each project ships only the columns its biotab forms
    contain.

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
        out_path = project_dir / table_name / "data.parquet"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if table_name in TABULAR_TABLES:
            table = pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table_name])
        elif rows:
            # Flex-schema table — let pyarrow infer column types per project.
            # Empty-table case skipped: zero-row biotab forms (e.g. LAML's
            # missing follow_up) write nothing rather than a 0×0 parquet.
            table = pa.Table.from_pylist(rows)
        else:
            continue
        pq.write_table(
            table,
            out_path,
            row_group_size=_ROW_GROUP_SIZE_BY_TABLE.get(table_name, _ROW_GROUP_SIZE_DEFAULT),
            write_page_index=True,
        )
        out_paths[table_name] = out_path
    return out_paths
