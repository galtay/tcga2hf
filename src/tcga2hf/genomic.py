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

from tcga2hf.gdc import GDCClient, and_, eq

# Fields we always pull on each /files hit so the manifest carries enough
# provenance to build patient-keyed tables later without re-querying GDC.
FILE_FIELDS: list[str] = [
    "file_id",
    "file_name",
    "file_size",
    "md5sum",
    "data_category",
    "data_type",
    "data_format",
    "experimental_strategy",
    "workflow_type",
    "access",
    "cases.case_id",
    "cases.submitter_id",
    "cases.project.project_id",
    "cases.samples.sample_id",
    "cases.samples.submitter_id",
    "cases.samples.sample_type",
    "cases.samples.tissue_type",
    "cases.samples.portions.analytes.aliquots.aliquot_id",
    "cases.samples.portions.analytes.aliquots.submitter_id",
]


def list_open_files(
    client: GDCClient,
    project_id: str,
    data_type: str,
    extra_filters: list[dict[str, Any]] | None = None,
) -> list[dict[str, Any]]:
    """Return all open-access /files hits for the given project + data_type."""
    clauses = [
        eq("cases.project.project_id", project_id),
        eq("access", "open"),
        eq("data_type", data_type),
    ]
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
) -> list[dict[str, Any]]:
    """Download every hit's bytes to <out_dir>/<file_name> and write manifest.json.

    Returns the manifest list. Existing files (matching size) are skipped unless
    skip_existing=False. Downloads are batched via POST /data (`batch_size` per
    request) so wall time scales with total bytes rather than per-file latency.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    hits = list_open_files(client, project_id, data_type, extra_filters)

    cached_ids: set[str] = set()
    to_download_ids: list[str] = []
    for hit in hits:
        target = out_dir / hit["file_name"]
        if skip_existing and target.exists() and target.stat().st_size == hit.get("file_size"):
            cached_ids.add(hit["file_id"])
        else:
            to_download_ids.append(hit["file_id"])

    if len(to_download_ids) >= 2:
        client.bulk_download(to_download_ids, out_dir, batch_size=batch_size)
    elif len(to_download_ids) == 1:
        only = next(h for h in hits if h["file_id"] == to_download_ids[0])
        client.download(only["file_id"], out_dir / only["file_name"])

    manifest: list[dict[str, Any]] = []
    for hit in hits:
        manifest.append(
            {
                "file_id": hit["file_id"],
                "file_name": hit["file_name"],
                "file_size": hit.get("file_size"),
                "md5sum": hit.get("md5sum"),
                "data_type": hit.get("data_type"),
                "data_format": hit.get("data_format"),
                "experimental_strategy": hit.get("experimental_strategy"),
                "workflow_type": hit.get("workflow_type"),
                "cases": hit.get("cases", []),
                "_status": "cached" if hit["file_id"] in cached_ids else "downloaded",
            }
        )

    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest
