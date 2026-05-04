"""Parse open-access GDC Gene Expression Quantification TSVs (RNA-Seq STAR counts).

Each TSV is one RNA-Seq aliquot's STAR-aligned gene quantifications: ~60,660
GENCODE v36 gene rows × 9 columns, plus 4 N_* summary rows at the top. We:

  - Lift the 4 N_* QC rows out as scalar fields (using their `unstranded` value).
  - Keep all 4 unstranded value columns (raw counts + TPM + FPKM + FPKM-UQ).
  - Drop `stranded_first` and `stranded_second` (TCGA used unstranded library
    prep, so the stranded counts are uninformative).
  - Preserve full gene metadata (gene_id with GENCODE version suffix, gene_name,
    gene_type) so the dataset is self-describing without a side vocab file.

GDC pipeline reference:
  https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

_QC_NAMES = ["N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"]


def _file_aliquot_and_case(entry: dict[str, Any]) -> tuple[str | None, str | None]:
    """Return (case_id, aliquot_id) for one manifest entry.

    Expression files are per-aliquot, so we expect exactly one case and one
    aliquot. Returns (None, None) if the manifest's nested cases data is
    ambiguous — caller should skip the file.
    """
    case_ids: set[str] = set()
    aliquot_ids: set[str] = set()
    for case in entry.get("cases", []):
        case_ids.add(case["case_id"])
        for sample in case.get("samples", []):
            for portion in sample.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for aliquot in analyte.get("aliquots") or []:
                        aliquot_ids.add(aliquot["aliquot_id"])
    if len(case_ids) != 1 or len(aliquot_ids) != 1:
        return None, None
    return next(iter(case_ids)), next(iter(aliquot_ids))


def _parse_expression_tsv(path: Path) -> tuple[dict[str, int | None], pd.DataFrame]:
    """Read the TSV; return ({QC_name: int_or_None}, gene_rows_df).

    The TSV starts with one `# gene-model: GENCODE v36` line (skipped via
    `comment="#"`), then a header line, then 4 N_* summary rows, then the
    ENSG gene rows. N_* rows have empty `gene_name`/`gene_type` cells.
    """
    df = pd.read_csv(path, sep="\t", comment="#", low_memory=False)
    qc_mask = df["gene_id"].str.startswith("N_", na=False)
    qc_rows = df[qc_mask].set_index("gene_id")
    gene_rows = df[~qc_mask].reset_index(drop=True)

    qc: dict[str, int | None] = {}
    for name in _QC_NAMES:
        if name in qc_rows.index and pd.notna(qc_rows.loc[name, "unstranded"]):
            qc[name] = int(qc_rows.loc[name, "unstranded"])
        else:
            qc[name] = None
    return qc, gene_rows


def _to_int_list(series: pd.Series) -> list[int | None]:
    """Pandas float-with-NaN -> Python list of int-or-None for pyarrow."""
    out: list[int | None] = []
    for v in series:
        if pd.isna(v):
            out.append(None)
        else:
            out.append(int(v))
    return out


def _to_float_list(series: pd.Series) -> list[float | None]:
    out: list[float | None] = []
    for v in series:
        if pd.isna(v):
            out.append(None)
        else:
            out.append(float(v))
    return out


def _build_record(
    qc: dict[str, int | None],
    gene_rows: pd.DataFrame,
    *,
    aliquot_id: str,
    source_file_id: str,
) -> dict[str, Any]:
    return {
        "sample_id": None,  # resolved at attach time from the patient's samples
        "aliquot_id": aliquot_id,
        "source_file_id": source_file_id,
        "N_unmapped": qc["N_unmapped"],
        "N_multimapping": qc["N_multimapping"],
        "N_noFeature": qc["N_noFeature"],
        "N_ambiguous": qc["N_ambiguous"],
        "gene_id": gene_rows["gene_id"].tolist(),
        "gene_name": gene_rows["gene_name"].tolist(),
        "gene_type": gene_rows["gene_type"].tolist(),
        "unstranded": _to_int_list(gene_rows["unstranded"]),
        "tpm_unstranded": _to_float_list(gene_rows["tpm_unstranded"]),
        "fpkm_unstranded": _to_float_list(gene_rows["fpkm_unstranded"]),
        "fpkm_uq_unstranded": _to_float_list(gene_rows["fpkm_uq_unstranded"]),
    }


def load_for_project(project_raw_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Return {case_id: [expression_record, ...]} for raw/<PROJECT>/expression/.

    Empty dict if the directory or its manifest is missing — letting `build`
    proceed without expression data when it hasn't been fetched.
    """
    expr_dir = project_raw_dir / "expression"
    manifest_path = expr_dir / "manifest.json"
    if not manifest_path.exists():
        return {}

    manifest = json.loads(manifest_path.read_text())
    by_case: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in manifest:
        file_path = expr_dir / entry["file_name"]
        if not file_path.exists():
            continue
        case_id, aliquot_id = _file_aliquot_and_case(entry)
        if not case_id or not aliquot_id:
            continue
        qc, gene_rows = _parse_expression_tsv(file_path)
        record = _build_record(
            qc, gene_rows, aliquot_id=aliquot_id, source_file_id=entry["file_id"]
        )
        by_case[case_id].append(record)
    return dict(by_case)


def attach(
    rows: list[dict[str, Any]], by_case: dict[str, list[dict[str, Any]]]
) -> list[dict[str, Any]]:
    """Mutate `rows` to populate samples_gene_expression_quantification.

    Resolves each record's sample_id by looking up its aliquot_id in the
    patient row's `samples` field (source of truth for every aliquot the case
    has). Sorts records by aliquot_id for deterministic output.
    """
    for row in rows:
        records = by_case.get(row["case_id"], [])
        # Walk the full GDC tree: samples → portions → analytes → aliquots
        aliquot_to_sample: dict[str, str] = {}
        for s in row.get("samples") or []:
            sid = s.get("sample_id")
            if not sid:
                continue
            for portion in s.get("portions") or []:
                for analyte in portion.get("analytes") or []:
                    for a in analyte.get("aliquots") or []:
                        aq = a.get("aliquot_id")
                        if aq:
                            aliquot_to_sample[aq] = sid
        for r in records:
            r["sample_id"] = aliquot_to_sample.get(r["aliquot_id"])
        records.sort(key=lambda r: r.get("aliquot_id") or "")
        row["samples_gene_expression_quantification"] = records
    return rows
