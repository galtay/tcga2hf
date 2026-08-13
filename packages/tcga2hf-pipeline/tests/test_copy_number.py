"""Tests for the two copy-number segment tables.

The interesting logic is FK resolution, not parsing: allele-specific files
name two aliquots and we have to pick the tumour one correctly, while masked
files name one. Both identify their subject through the file's own
`GDC_Aliquot` column rather than by reading sample-type digits out of a
barcode, and these tests pin that down.
"""

from __future__ import annotations

import json
from pathlib import Path

import pyarrow as pa
from tcga2hf.schema import TABULAR_TABLES
from tcga2hf_pipeline import tabular

_CASE_UUID = "5fb3affa-3661-48e8-9d0f-8f0f7f6b0f11"
_SAMPLE_UUID = "e235e45d-9cbe-4c11-9a3a-1cf1d1e0aa10"
_NORMAL_SAMPLE_UUID = "aa1d5a74-0687-4d12-abcc-7f21c2e9fb00"
_TUMOR_ALIQUOT = "b3c68473-78a7-44d2-95bf-7997a7e6c0e8"
_NORMAL_ALIQUOT = "451f8cc6-05fe-48f6-9f2e-1b8c0d3a77aa"

_ALLELE_SPECIFIC_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tCopy_Number\tMajor_Copy_Number\tMinor_Copy_Number\n"
    f"{_TUMOR_ALIQUOT}\tchr1\t62920\t33086177\t1\t1\t0\n"
    f"{_TUMOR_ALIQUOT}\tchr1\t33086544\t64924592\t3\t2\t1\n"
)

_MASKED_TSV = (
    "GDC_Aliquot\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n"
    f"{_TUMOR_ALIQUOT}\t1\t3301765\t6028611\t1917\t-0.2829\n"
    f"{_TUMOR_ALIQUOT}\tX\t6040795\t6041047\t2\t-2.0991\n"
)


def _case() -> dict:
    """A GDC case dict with one tumour sample and one matched normal sample."""
    return {
        "case_id": _CASE_UUID,
        "submitter_id": "TCGA-W5-AA33",
        "samples": [
            {
                "sample_id": _SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-01A",
                "sample_type": "Primary Tumor",
                "portions": [
                    {
                        "portion_id": "portion-tumor",
                        "submitter_id": "TCGA-W5-AA33-01A-11",
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _TUMOR_ALIQUOT,
                                        "submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
                                    }
                                ]
                            }
                        ],
                    }
                ],
            },
            {
                "sample_id": _NORMAL_SAMPLE_UUID,
                "submitter_id": "TCGA-W5-AA33-10A",
                "sample_type": "Blood Derived Normal",
                "portions": [
                    {
                        "portion_id": "portion-normal",
                        "submitter_id": "TCGA-W5-AA33-10A-01",
                        "analytes": [
                            {
                                "aliquots": [
                                    {
                                        "aliquot_id": _NORMAL_ALIQUOT,
                                        "submitter_id": "TCGA-W5-AA33-10A-01D-A419-01",
                                    }
                                ]
                            }
                        ],
                    }
                ],
            },
        ],
    }


def _build_project(
    tmp_path: Path,
    *,
    modality: str,
    file_name: str,
    content: str,
    entities: list[dict],
    workflow_type: str,
    experimental_strategy: str | None = None,
    status: str = "downloaded",
    write_file: bool = True,
) -> Path:
    project_dir = tmp_path / "TCGA-XYZ"
    mod_dir = project_dir / modality
    mod_dir.mkdir(parents=True)
    if write_file:
        (mod_dir / file_name).write_text(content)
    entry = {
        "file_id": f"synthetic-{modality}-file-id",
        "file_name": file_name,
        "workflow_type": workflow_type,
        "experimental_strategy": experimental_strategy,
        "cases": [{"case_id": _CASE_UUID, "submitter_id": "TCGA-W5-AA33"}],
        "associated_entities": entities,
        "_status": status,
    }
    (mod_dir / "manifest.json").write_text(json.dumps([entry]))
    return project_dir


def _paired_entities() -> list[dict]:
    """Both aliquots, normal listed first — order must not decide the tumour."""
    return [
        {
            "entity_id": _NORMAL_ALIQUOT,
            "entity_type": "aliquot",
            "entity_submitter_id": "TCGA-W5-AA33-10A-01D-A419-01",
        },
        {
            "entity_id": _TUMOR_ALIQUOT,
            "entity_type": "aliquot",
            "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
        },
    ]


def _allele_specific_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="TCGA-CHOL.aliquot.ascat3.allelic_specific.seg.txt",
        content=_ALLELE_SPECIFIC_TSV,
        entities=_paired_entities(),
        workflow_type="ASCAT3",
        experimental_strategy="Genotyping Array",
        **kwargs,
    )


def _masked_project(tmp_path: Path, **kwargs) -> Path:
    return _build_project(
        tmp_path,
        modality="copy_number_masked",
        file_name="SPIKY_p_TCGAb.nocnv_grch38.seg.v2.txt",
        content=_MASKED_TSV,
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
        workflow_type="DNAcopy",
        experimental_strategy="Genotyping Array",
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Allele-specific
# ---------------------------------------------------------------------------


def test_tumor_aliquot_comes_from_the_file_not_the_entity_order(tmp_path: Path) -> None:
    """The normal is listed first in `associated_entities`; the tumour still wins.

    This is the whole point of reading `GDC_Aliquot`: entity order is not a
    contract, and inferring tumour-vs-normal from barcode digits would be.
    """
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)

    assert len(rows) == 2
    for row in rows:
        assert row["aliquot_id"] == _TUMOR_ALIQUOT
        assert row["aliquot_submitter_id"] == "TCGA-W5-AA33-01A-11D-A416-01"
        assert row["matched_normal_aliquot_id"] == _NORMAL_ALIQUOT
        assert row["matched_normal_aliquot_submitter_id"] == "TCGA-W5-AA33-10A-01D-A419-01"
        # The sample FK must follow the tumour aliquot, not the normal one.
        assert row["sample_id"] == _SAMPLE_UUID
        assert row["sample_type"] == "Primary Tumor"


def test_allele_specific_carries_values_and_workflow_verbatim(tmp_path: Path) -> None:
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)

    first, second = rows
    assert (first["chromosome"], first["start"], first["end"]) == ("chr1", 62920, 33086177)
    assert (first["copy_number"], first["major_copy_number"], first["minor_copy_number"]) == (
        1,
        1,
        0,
    )
    assert second["copy_number"] == second["major_copy_number"] + second["minor_copy_number"]
    assert first["workflow_type"] == "ASCAT3"
    assert first["experimental_strategy"] == "Genotyping Array"


def test_allele_specific_keeps_the_chr_prefix(tmp_path: Path) -> None:
    """`chr1` here vs bare `1` in the masked table — both as the source writes them."""
    project = _allele_specific_project(tmp_path)
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)
    assert {r["chromosome"] for r in rows} == {"chr1"}


def test_ambiguous_gdc_aliquot_column_is_skipped_not_guessed(tmp_path: Path) -> None:
    """A file naming two different aliquots has no defensible tumour FK."""
    mixed = _ALLELE_SPECIFIC_TSV + f"{_NORMAL_ALIQUOT}\tchr2\t1\t100\t2\t1\t1\n"
    project = _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="ambiguous.seg.txt",
        content=mixed,
        entities=_paired_entities(),
        workflow_type="ASCAT3",
    )
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []


def test_missing_matched_normal_leaves_null_rather_than_a_guess(tmp_path: Path) -> None:
    project = _build_project(
        tmp_path,
        modality="copy_number_allele_specific",
        file_name="unpaired.seg.txt",
        content=_ALLELE_SPECIFIC_TSV,
        entities=[
            {
                "entity_id": _TUMOR_ALIQUOT,
                "entity_type": "aliquot",
                "entity_submitter_id": "TCGA-W5-AA33-01A-11D-A416-01",
            }
        ],
        workflow_type="AscatNGS",
    )
    rows = tabular._allele_specific_copy_number_segment_rows([_case()], project)
    assert rows
    assert all(r["matched_normal_aliquot_id"] is None for r in rows)
    assert all(r["aliquot_id"] == _TUMOR_ALIQUOT for r in rows)


# ---------------------------------------------------------------------------
# Masked
# ---------------------------------------------------------------------------


def test_masked_segments_keep_bare_chromosome_names_as_strings(tmp_path: Path) -> None:
    """`1` must stay the string "1" — not the integer 1, and X must survive."""
    project = _masked_project(tmp_path)
    rows = tabular._masked_copy_number_segment_rows([_case()], project)

    assert [r["chromosome"] for r in rows] == ["1", "X"]
    assert all(isinstance(r["chromosome"], str) for r in rows)


def test_masked_carries_probe_counts_and_segment_means(tmp_path: Path) -> None:
    project = _masked_project(tmp_path)
    rows = tabular._masked_copy_number_segment_rows([_case()], project)

    assert (rows[0]["num_probes"], rows[0]["segment_mean"]) == (1917, -0.2829)
    assert rows[0]["aliquot_id"] == _TUMOR_ALIQUOT
    assert rows[0]["workflow_type"] == "DNAcopy"


# ---------------------------------------------------------------------------
# Shared contracts
# ---------------------------------------------------------------------------


def test_manifest_only_entries_are_skipped(tmp_path: Path) -> None:
    """A capped fetch lists files it didn't download; they have no bytes to read."""
    project = _allele_specific_project(tmp_path, status="manifest_only", write_file=False)
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []


def test_missing_modality_dir_yields_no_rows(tmp_path: Path) -> None:
    project = tmp_path / "TCGA-EMPTY"
    project.mkdir()
    assert tabular._allele_specific_copy_number_segment_rows([_case()], project) == []
    assert tabular._masked_copy_number_segment_rows([_case()], project) == []


def test_unknown_case_is_skipped(tmp_path: Path) -> None:
    project = _masked_project(tmp_path)
    assert tabular._masked_copy_number_segment_rows([], project) == []


def test_rows_satisfy_the_tabular_schemas(tmp_path: Path) -> None:
    """Emitted rows must load into the published schemas without coercion."""
    allele = tabular._allele_specific_copy_number_segment_rows(
        [_case()], _allele_specific_project(tmp_path / "a")
    )
    masked = tabular._masked_copy_number_segment_rows([_case()], _masked_project(tmp_path / "m"))

    for rows, table in (
        (allele, "allele_specific_copy_number_segment"),
        (masked, "masked_copy_number_segment"),
    ):
        built = pa.Table.from_pylist(rows, schema=TABULAR_TABLES[table])
        assert built.num_rows == len(rows)
        assert built.schema.equals(TABULAR_TABLES[table])
