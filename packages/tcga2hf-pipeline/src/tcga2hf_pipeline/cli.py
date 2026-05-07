from __future__ import annotations

import json
import os
import shutil
from pathlib import Path
from typing import Annotated

import typer
from dotenv import load_dotenv
from tcga2hf.schema import TABULAR_TABLES

from tcga2hf_pipeline import (
    cdr,
    clinical,
    dataset_card,
    expression,
    genomic,
    hf_upload,
    mutations,
    tabular,
)
from tcga2hf_pipeline.gdc import GDCClient, write_cases_json

# Load .env from cwd (or any parent) on import. override=True so the project's
# .env wins over any inherited shell variable: the HF_TOKEN here is scoped to
# this project, and we don't want a stale global token to silently take over.
load_dotenv(override=True)

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Download public TCGA data from the NCI GDC and stage it for the HF Hub.",
)

DEFAULT_DATA_DIR = Path.home() / "data" / "tcga2hf"


def _resolve_data_dir(data_dir: Path | None) -> Path:
    if data_dir is not None:
        return data_dir
    env = os.environ.get("TCGA2HF_DATA_DIR")
    if env:
        return Path(env).expanduser()
    return DEFAULT_DATA_DIR


DataDirOpt = Annotated[
    Path | None,
    typer.Option(
        "--data-dir",
        help="Root data dir. Defaults to $TCGA2HF_DATA_DIR or $HOME/data/tcga2hf.",
    ),
]


@app.command("fetch-clinical")
def fetch_clinical_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable), e.g. --project TCGA-CHOL --project TCGA-DLBC.",
        ),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """Fetch clinical case JSON for one or more TCGA projects from the GDC."""
    import hashlib

    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        # Snapshot the GDC server's data release once at the top so every project
        # fetched in this invocation pairs with the same release tag.
        status = client.status()
        release_str = status.get("data_release", "<unknown>")
        typer.echo(f"GDC: {release_str} (tag {status.get('tag', '?')})")

        # Capture the dictionary alongside, tagged by release. Reuses the cached
        # copy if a fetch with the same release already happened locally.
        version = status.get("data_release_version") or {}
        major = version.get("major", "x")
        minor = version.get("minor", "y")
        dict_path = raw_dir / f"gdc_dictionary.{major}.{minor}.json"
        if not dict_path.exists():
            typer.echo(f"capturing GDC dictionary -> {dict_path.name} ...")
            dict_path.parent.mkdir(parents=True, exist_ok=True)
            dictionary = client.dictionary()
            dict_text = json.dumps(dictionary, indent=2, sort_keys=True)
            dict_path.write_text(dict_text)
        else:
            dict_text = dict_path.read_text()
        # SHA-256 of the (sorted) dictionary JSON pins provenance precisely
        # without needing to diff 10MB blobs across runs.
        status["_gdc_dictionary_sha256"] = hashlib.sha256(dict_text.encode()).hexdigest()
        status["_gdc_dictionary_path"] = dict_path.name

        for proj in project:
            typer.echo(f"fetching {proj} ...")
            cases = clinical.fetch_clinical([proj], client)
            out_path = raw_dir / proj / "cases.json"
            write_cases_json(cases, out_path)
            (raw_dir / proj / "gdc_status.json").write_text(json.dumps(status, indent=2))
            typer.echo(f"  wrote {len(cases)} cases -> {out_path}")


def _fetch_modality(
    project: list[str],
    data_dir: Path | None,
    data_type: str,
    modality_dir: str,
    max_files: int | None = None,
) -> None:
    """Shared body for fetch-mutations / fetch-expression / future modalities."""
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')} (tag {status.get('tag', '?')})")
        for proj in project:
            out_dir = raw_dir / proj / modality_dir
            typer.echo(f"fetching {data_type!r} for {proj} -> {out_dir}")
            manifest = genomic.fetch_files(client, proj, data_type, out_dir, max_files=max_files)
            n_dl = sum(1 for m in manifest if m["_status"] == "downloaded")
            n_cache = sum(1 for m in manifest if m["_status"] == "cached")
            n_skip = sum(1 for m in manifest if m["_status"] == "manifest_only")
            total_mb = sum((m.get("file_size") or 0) for m in manifest) / 1e6
            extra = f", {n_skip} manifest-only" if n_skip else ""
            typer.echo(
                f"  {len(manifest):>4} files ({n_dl} downloaded, {n_cache} cached{extra}), "
                f"{total_mb:.1f} MB total in manifest"
            )
            (raw_dir / proj / "gdc_status.json").write_text(json.dumps(status, indent=2))


MaxFilesOpt = Annotated[
    int | None,
    typer.Option(
        "--max-files",
        help=(
            "Cap on total files on disk (cached + new) per project. Set to 1 to "
            "sample one file per project, or 0 to populate the manifest without "
            "downloading any new bytes. Manifest always lists every discovered file."
        ),
    ),
]


@app.command("fetch-mutations")
def fetch_mutations_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access Masked Somatic Mutation MAFs.",  # noqa: E501
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Masked Somatic Mutation MAF files (DNA)."""
    _fetch_modality(project, data_dir, "Masked Somatic Mutation", "mutations", max_files=max_files)


@app.command("fetch-expression")
def fetch_expression_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access RNA-Seq STAR gene-count TSVs.",  # noqa: E501
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Gene Expression Quantification TSVs (RNA-Seq, STAR counts)."""
    _fetch_modality(
        project, data_dir, "Gene Expression Quantification", "expression", max_files=max_files
    )


@app.command("fetch-cdr")
def fetch_cdr_cmd(
    data_dir: DataDirOpt = None,
) -> None:
    """Fetch the Liu et al. 2018 PanCanAtlas CDR workbook (curated TCGA survival).

    Downloads the GDC's PanCanAtlas auxiliary file (UUID
    1b5f413e-...) into <data-dir>/raw/cdr/. md5-pinned; safe to re-run.
    The file is frozen at the 2018 data freeze, so coverage stops there
    — post-2018 cases get `cdr_matched=False` at build time.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")
    out_path = cdr.fetch_cdr_workbook(raw_dir)
    typer.echo(f"  wrote -> {out_path}")
    index = cdr.load_cdr_index(raw_dir)
    typer.echo(f"  CDR rows indexed: {len(index)}")


@app.command("build")
def build_cmd(
    data_dir: DataDirOpt = None,
) -> None:
    """Flatten raw clinical JSON into per-project patient Parquets + dataset card."""
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = sorted(raw_dir.glob("*/cases.json"))
    if not project_files:
        raise typer.BadParameter(f"no cases.json files found under {raw_dir}.")

    # Wipe the whole processed/ tree so removed projects + any legacy layout
    # (e.g. the old patients/ directory) don't leave stale data behind.
    if processed_dir.exists():
        shutil.rmtree(processed_dir)
    processed_dir.mkdir(parents=True)

    # CDR (Liu et al 2018 curated survival) is shared across all projects;
    # load once at the top. If `fetch-cdr` hasn't been run, the index is
    # empty and every patient row gets cdr_matched=False — the build
    # still works.
    cdr_index = cdr.load_cdr_index(raw_dir)
    typer.echo(f"CDR rows indexed: {len(cdr_index)}")

    projects: list[str] = []
    gdc_releases: dict[str, str] = {}
    for path in project_files:
        proj = path.parent.name
        projects.append(proj)
        cases = json.loads(path.read_text())
        rows = clinical.to_patient_rows(cases)

        # Attach molecular modalities if their raw data has been fetched.
        mut_by_case = mutations.load_for_project(path.parent)
        expr_by_case = expression.load_for_project(path.parent)
        if mut_by_case:
            mutations.attach(rows, mut_by_case)
        if expr_by_case:
            expression.attach(rows, expr_by_case)
        # Liu's CDR survival columns; populate even when cdr_index is
        # empty so the parquet schema stays stable.
        cdr.attach_cdr(rows, cdr_index)

        n_variants = sum(len(r["samples_masked_somatic_mutation"]) for r in rows)
        n_expr = sum(len(r["samples_gene_expression_quantification"]) for r in rows)
        n_cdr = sum(1 for r in rows if r["cdr_matched"])
        typer.echo(
            f"  {proj:<12} {len(rows):>4} patients  "
            f"mutations={n_variants:>5} ({len(mut_by_case)} MAFs)  "
            f"expression={n_expr:>4} aliquots ({len(expr_by_case)} TSVs)  "
            f"cdr={n_cdr}/{len(rows)}"
        )

        out = clinical.write_patients(rows, processed_dir, proj)
        typer.echo(f"             -> {out}")

        status_path = path.parent / "gdc_status.json"
        if status_path.exists():
            status = json.loads(status_path.read_text())
            gdc_releases[proj] = status.get("data_release", "<unknown>")

    card = dataset_card.write_card(processed_dir, projects, gdc_releases=gdc_releases)
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed tree at: {processed_dir}")


@app.command("build-tabular")
def build_tabular_cmd(
    project: Annotated[
        list[str] | None,
        typer.Option(
            "--project",
            help=(
                "TCGA project id (repeatable). If omitted, all projects under "
                "<data-dir>/raw/ are processed."
            ),
        ),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Flatten raw TCGA data into per-(project, table) Parquets + dataset card.

    Companion to `build`: same raw inputs, different output shape — one HF
    subset per (project, table) under <data-dir>/processed_tabular/.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_tabular"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = sorted(raw_dir.glob("*/cases.json"))
    if project:
        wanted = set(project)
        project_files = [p for p in project_files if p.parent.name in wanted]
        missing = wanted - {p.parent.name for p in project_files}
        if missing:
            raise typer.BadParameter(
                f"requested projects missing from raw/: {sorted(missing)}"
            )
    if not project_files:
        raise typer.BadParameter(f"no cases.json files found under {raw_dir}.")

    # Wipe any prior tabular tree so removed projects/tables don't leave
    # stale data behind.
    if processed_dir.exists():
        shutil.rmtree(processed_dir)
    processed_dir.mkdir(parents=True)

    cdr_index = cdr.load_cdr_index(raw_dir)
    typer.echo(f"CDR rows indexed: {len(cdr_index)}")

    projects: list[str] = []
    gdc_releases: dict[str, str] = {}
    for path in project_files:
        proj = path.parent.name
        projects.append(proj)
        cases = json.loads(path.read_text())
        tables = tabular.build_tables(cases, path.parent)
        # The `cases` table inherits the cdr_* placeholders from
        # `clinical._patient_row`; populate them from CDR (no-op for
        # cases not in Liu's 2018 freeze).
        cdr.attach_cdr(tables["cases"], cdr_index)

        sizes = {name: len(rows) for name, rows in tables.items()}
        # Compact one-line summary of row counts per table — easier to spot
        # cardinality regressions than scrolling through 14 lines per project.
        summary = " ".join(f"{name}={n}" for name, n in sizes.items())
        n_cdr = sum(1 for r in tables["cases"] if r["cdr_matched"])
        typer.echo(f"  {proj:<12} {summary}  cdr={n_cdr}/{len(tables['cases'])}")

        out_paths = tabular.write_tables(tables, processed_dir, proj)
        # Use any one table's path to print the project's output dir.
        typer.echo(f"             -> {next(iter(out_paths.values())).parent.parent}")

        status_path = path.parent / "gdc_status.json"
        if status_path.exists():
            status = json.loads(status_path.read_text())
            gdc_releases[proj] = status.get("data_release", "<unknown>")

    card = dataset_card.write_tabular_card(
        processed_dir,
        projects,
        list(TABULAR_TABLES),
        gdc_releases=gdc_releases,
    )
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed_tabular tree at: {processed_dir}")


@app.command("upload")
def upload_cmd(
    repo_id: Annotated[
        str,
        typer.Option(
            "--repo-id",
            help="HF dataset repo id, e.g. galtay/tcga-patients-open.",
        ),
    ],
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private (default) or public. Recommend private for first pushes.",
        ),
    ] = True,
    commit_message: Annotated[
        str | None,
        typer.Option(
            "--commit-message",
            "-m",
            help="Commit message for this upload.",
        ),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push <data-dir>/processed/ to a HuggingFace dataset repo."""
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed"
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")
    if not private:
        typer.confirm(
            f"This will publish {repo_id} as PUBLIC. Continue?",
            abort=True,
        )

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message,
    )
    typer.echo(f"\nuploaded -> {url}")


@app.command("upload-tabular")
def upload_tabular_cmd(
    repo_id: Annotated[
        str,
        typer.Option(
            "--repo-id",
            help="HF dataset repo id, e.g. gabrielaltay/tcga-tabular-open.",
        ),
    ],
    private: Annotated[
        bool,
        typer.Option(
            "--private/--public",
            help="Upload as private (default) or public.",
        ),
    ] = True,
    commit_message: Annotated[
        str | None,
        typer.Option("--commit-message", "-m", help="Commit message for this upload."),
    ] = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Push <data-dir>/processed_tabular/ to a HuggingFace dataset repo.

    Companion to `upload`: same shape, just points at the tabular tree.
    """
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_tabular"
    visibility = "private" if private else "PUBLIC"
    typer.echo(f"processed dir: {processed_dir}")
    typer.echo(f"repo_id:       {repo_id} ({visibility})")
    if not private:
        typer.confirm(
            f"This will publish {repo_id} as PUBLIC. Continue?",
            abort=True,
        )

    url = hf_upload.upload_dataset(
        processed_dir=processed_dir,
        repo_id=repo_id,
        private=private,
        commit_message=commit_message,
    )
    typer.echo(f"\nuploaded -> {url}")


if __name__ == "__main__":
    app()
