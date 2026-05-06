"""Generic open-access GDC file fetcher for genomic / molecular data.

Mutations and expression (and miRNA, copy number, etc. later) all share the same
shape: query /files with a project + data_type filter, download each hit's bytes,
and write a manifest mapping file_id -> case_submitter_id + sample info. The
patient-row build step can join on those FKs later — we keep raw downloads and
manifest separate so iteration on the schema doesn't trigger re-downloads.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from tcga2hf_pipeline.gdc import GDCClient, and_, eq

# Fields we request on each /files hit. The GDC `/files` endpoint
# *replaces* the default response with whatever you list here (it does
# NOT append to a default), so this enumerates every field our pipeline
# consumes — file metadata we save into the manifest, plus the case and
# aliquot linkages we need for downstream join logic. Trimmed to only
# what's actually read; cf. `expression._file_aliquot_and_case` for the
# reason `cases.samples.portions.analytes.aliquots.aliquot_id` is in
# here despite looking deep.
#
# `analysis.workflow_type` is the canonical GDC path — the unprefixed
# `workflow_type` comes back as an unrecognized field warning. Flattened
# to a top-level `workflow_type` column at manifest-write time.
FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "access",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "analysis.workflow_type",
    "cases.case_id",
    "cases.submitter_id",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
]


# Per-modality filter values that lock the `/files` request to the exact
# tooling we expect. Filtering on `data_type` alone is technically enough
# for TCGA today, but adding `data_format` + `experimental_strategy` +
# `analysis.workflow_type` guards against future GDC additions silently
# shipping different workflows under the same data_type. Values verified
# against live `/files` responses on 2026-05; if any of them changes
# upstream the filter will return zero hits and surface the regression
# loudly.
MODALITY_FILTERS: dict[str, dict[str, str]] = {
    "Masked Somatic Mutation": {
        "data_format": "MAF",
        "data_category": "Simple Nucleotide Variation",
        "experimental_strategy": "WXS",
        "analysis.workflow_type": "Aliquot Ensemble Somatic Variant Merging and Masking",
    },
    "Gene Expression Quantification": {
        "data_format": "TSV",
        "data_category": "Transcriptome Profiling",
        "experimental_strategy": "RNA-Seq",
        "analysis.workflow_type": "STAR - Counts",
    },
}


def list_open_files(
    client: GDCClient,
    project_id: str,
    data_type: str,
    extra_filters: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Return all open-access /files hits for the given project + data_type.

    Filter clauses:
      - cases.project.project_id = <project_id>
      - access = "open"
      - data_type = <data_type>
      - …plus every key/value pair in `MODALITY_FILTERS[data_type]` (locks
        format, experimental_strategy, workflow_type for the modalities we
        explicitly support).
    """
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
        eq("data_type", data_type),
    ]
    for field, value in MODALITY_FILTERS.get(data_type, {}).items():
        clauses.append(eq(field, value))
    if extra_filters:
        clauses.extend(extra_filters)
    return client.files(filters=and_(*clauses), fields=FILE_FIELDS, page_size=500)


def fetch_files(
    client: GDCClient,
    project_id: str,
    data_type: str,
    out_dir: Path,
    extra_filters: list[dict[str, Any]] | None = None,
    skip_existing: bool = True,
    batch_size: int = 50,
    max_files: int | None = None,
) -> list[dict[str, Any]]:
    """Download hits' bytes to <out_dir>/<file_name> and write manifest.json.

    Returns the manifest list. Existing files (matching size) are skipped
    unless skip_existing=False. Downloads are batched via POST /data
    (`batch_size` per request) so wall time scales with total bytes rather
    than per-file latency.

    `max_files` caps the *total* number of files we'll have on disk after
    this call (cached + freshly downloaded). Set to 0 to discover and
    write the full manifest without downloading any new bytes — useful for
    populating the `files` table across projects where you only want a
    sampler of the actual molecular content. The manifest always lists
    every hit returned by `/files`, regardless of how many were
    downloaded; entries trimmed off by `max_files` are tagged
    `_status="manifest_only"`.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = list_open_files(client, project_id, data_type, extra_filters)

    cached_ids: set[str] = set()
    eligible_for_download: list[str] = []
    for hit in hits:
        target = out_dir / hit["file_name"]
        if skip_existing and target.exists() and target.stat().st_size == hit.get("file_size"):
            cached_ids.add(hit["file_id"])
        else:
            eligible_for_download.append(hit["file_id"])

    if max_files is None:
        to_download_ids = eligible_for_download
    else:
        # `max_files` is a target *total* count on disk, so subtract what's
        # already cached. If we already have ≥ max_files cached, download
        # nothing new.
        budget = max(0, max_files - len(cached_ids))
        to_download_ids = eligible_for_download[:budget]
    skipped_ids = set(eligible_for_download) - set(to_download_ids)

    if len(to_download_ids) >= 2:
        client.bulk_download(to_download_ids, out_dir, batch_size=batch_size)
    elif len(to_download_ids) == 1:
        only = next(h for h in hits if h["file_id"] == to_download_ids[0])
        client.download(only["file_id"], out_dir / only["file_name"])

    manifest: list[dict[str, Any]] = []
    for hit in hits:
        # GDC nests workflow info under `analysis.*`; flatten it into the
        # manifest so downstream consumers see a single `workflow_type`
        # column matching the dataset's `files` schema.
        analysis = hit.get("analysis") or {}
        if hit["file_id"] in cached_ids:
            status = "cached"
        elif hit["file_id"] in skipped_ids:
            status = "manifest_only"
        else:
            status = "downloaded"
        manifest.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "file_size": hit.get("file_size"),
                "md5sum": hit.get("md5sum"),
                "data_category": hit.get("data_category"),
                "data_type": hit.get("data_type"),
                "data_format": hit.get("data_format"),
                "experimental_strategy": hit.get("experimental_strategy"),
                "workflow_type": analysis.get("workflow_type"),
                "access": hit.get("access"),
                "cases": hit.get("cases", []),
                "_status": status,
            }
        )

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest
