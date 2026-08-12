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
    clinical_supplement,
    dataset_card,
    expression,
    genomic,
    hf_upload,
    mutations,
    pathology,
    survival,
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

ProjectFilterOpt = Annotated[
    list[str] | None,
    typer.Option(
        "--project",
        help=(
            "TCGA project id (repeatable). If omitted, all projects under "
            "<data-dir>/raw/ are processed and the output tree is rebuilt from "
            "scratch. If given, only those projects' output dirs are replaced — "
            "other projects' parquets are left untouched."
        ),
    ),
]


def _clear_processed(processed_dir: Path, project: list[str] | None) -> None:
    """Remove prior output for the projects about to be rebuilt.

    Without `--project` the whole tree goes, so projects dropped from raw/
    and any legacy layout don't leave stale files behind. With `--project`
    only the named projects' directories go: a full wipe there would delete
    parquets we aren't regenerating, and (because HF uploads diff against
    what's on disk) would turn an incremental append into a full re-upload.
    """
    if project is None:
        if processed_dir.exists():
            shutil.rmtree(processed_dir)
        processed_dir.mkdir(parents=True)
        return
    processed_dir.mkdir(parents=True, exist_ok=True)
    for proj in project:
        target = processed_dir / proj
        if target.exists():
            shutil.rmtree(target)


def _card_inputs(
    raw_dir: Path, processed_dir: Path, built_marker: str
) -> tuple[list[str], dict[str, str]]:
    """Return (projects, gdc_releases) for every project present in the output tree.

    The dataset card describes the tree as a whole, not just the projects
    this invocation rebuilt — otherwise a `--project` build would drop the
    other projects from the `configs:` block and orphan their parquets on
    the Hub. `built_marker` is a glob, relative to a project dir, that
    identifies a completed build for the layout being written.

    Each project's GDC release is read back from its own raw
    `gdc_status.json`, so projects fetched at different releases report
    their own — the card renders a per-project list rather than a single
    release when they diverge.
    """
    projects: list[str] = []
    gdc_releases: dict[str, str] = {}
    if not processed_dir.exists():
        return projects, gdc_releases
    for project_dir in sorted(p for p in processed_dir.iterdir() if p.is_dir()):
        if not any(project_dir.glob(built_marker)):
            continue
        proj = project_dir.name
        projects.append(proj)
        status_path = raw_dir / proj / "gdc_status.json"
        if status_path.exists():
            gdc_releases[proj] = json.loads(status_path.read_text()).get(
                "data_release", "<unknown>"
            )
    return projects, gdc_releases


def _select_projects(raw_dir: Path, project: list[str] | None) -> list[Path]:
    """Resolve `--project` to the cases.json paths to process (all if None)."""
    project_files = sorted(raw_dir.glob("*/cases.json"))
    if project:
        wanted = set(project)
        project_files = [p for p in project_files if p.parent.name in wanted]
        missing = wanted - {p.parent.name for p in project_files}
        if missing:
            raise typer.BadParameter(f"requested projects missing from raw/: {sorted(missing)}")
    if not project_files:
        raise typer.BadParameter(f"no cases.json files found under {raw_dir}.")
    return project_files


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
    """Shared body for fetch-mutations / fetch-expression / future modalities.

    The GDC release is recorded per modality, in the modality's own
    directory — not on the project. Modalities are fetched at different
    times (a new one added years later pulls whatever release GDC is
    serving then), so one status file per project would report only the
    most recent fetch and silently overwrite the release + dictionary hash
    that `fetch-clinical` recorded for `cases.json`.
    """
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
            (out_dir / "gdc_status.json").write_text(json.dumps(status, indent=2))


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


@app.command("fetch-pathology-reports")
def fetch_pathology_reports_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads open-access Pathology Report PDFs.",
        ),
    ],
    data_dir: DataDirOpt = None,
    max_files: MaxFilesOpt = None,
) -> None:
    """Download open-access Pathology Report PDFs (scanned BCR documents).

    One report per case's tumor sample, ~11,200 across TCGA (~2.6 GB; LAML
    has none). The PDFs are shipped into the dataset verbatim — see
    `tcga2hf_pipeline.pathology` for why no text extraction happens here.
    """
    _fetch_modality(
        project, data_dir, "Pathology Report", "pathology_reports", max_files=max_files
    )


@app.command("fetch-clinical-supplements")
def fetch_clinical_supplements_cmd(
    project: Annotated[
        list[str],
        typer.Option(
            "--project",
            help="TCGA project id (repeatable). Downloads BCR biotab Clinical Supplement files.",
        ),
    ],
    data_dir: DataDirOpt = None,
) -> None:
    """Download BCR biotab Clinical Supplement files (the original BCR forms).

    These TSV files carry clinical fields the harmonized `/cases?expand=...`
    API drops or under-populates — most importantly Liu et al. 2018's
    `treatment_outcome_first_course`, the disease-free signal that drives
    DFI re-derivation. Files land at `<data-dir>/raw/<project>/clinical_supplement/`.
    """
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    typer.echo(f"raw dir: {raw_dir}")

    with GDCClient() as client:
        status = client.status()
        typer.echo(f"GDC: {status.get('data_release', '<unknown>')} (tag {status.get('tag', '?')})")
        for proj in project:
            out_dir = raw_dir / proj / "clinical_supplement"
            typer.echo(f"fetching Clinical Supplement biotabs for {proj} -> {out_dir}")
            manifest = clinical_supplement.fetch_clinical_supplements(client, proj, out_dir)
            n_dl = sum(1 for m in manifest if m["_status"] == "downloaded")
            n_cache = sum(1 for m in manifest if m["_status"] == "cached")
            kinds = sorted({m["form_kind"] for m in manifest})
            typer.echo(
                f"  {len(manifest):>3} files ({n_dl} downloaded, {n_cache} cached)  "
                f"forms: {', '.join(kinds)}"
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
    project: ProjectFilterOpt = None,
    data_dir: DataDirOpt = None,
) -> None:
    """Flatten raw clinical JSON into per-project patient Parquets + dataset card."""
    root = _resolve_data_dir(data_dir)
    raw_dir = root / "raw"
    processed_dir = root / "processed_patient"
    typer.echo(f"raw dir:       {raw_dir}")
    typer.echo(f"processed dir: {processed_dir}")

    if not raw_dir.exists():
        raise typer.BadParameter(f"raw dir does not exist: {raw_dir}. Run fetch-clinical first.")

    project_files = _select_projects(raw_dir, project)
    _clear_processed(processed_dir, project)

    for path in project_files:
        proj = path.parent.name
        cases = json.loads(path.read_text())
        rows = clinical.to_patient_rows(cases)

        # Attach molecular modalities if their raw data has been fetched.
        mut_by_case = mutations.load_for_project(path.parent)
        expr_by_case = expression.load_for_project(path.parent)
        path_by_case = pathology.load_for_project(path.parent)
        if mut_by_case:
            mutations.attach(rows, mut_by_case)
        if expr_by_case:
            expression.attach(rows, expr_by_case)
        if path_by_case:
            pathology.attach(rows, path_by_case)
        # Attach BCR biotab Clinical Supplement data if it's been fetched.
        # survival.attach_survival reads `row["clinical_supplement"]` for the
        # ~2.4×-better-populated `treatment_outcome_first_course` field,
        # which drives DFI re-derivation. Supplement data is consumed
        # in-memory only; not serialized to the patients dataset.
        supp_dir = path.parent / "clinical_supplement"
        supps = clinical_supplement.load_supplements_for_project(supp_dir)
        if supps:
            clinical_supplement.attach_supplements(rows, supps)
        # Re-derive OS / DSS / PFI / DFI from current GDC data using Liu
        # et al. 2018's algorithm; results land in `survival_derived` struct.
        # See `dev_research/liu_2018/report.html` for validation against
        # Liu's curated 2018 CDR.
        survival.attach_survival(rows)

        n_variants = sum(len(r["samples_masked_somatic_mutation"]) for r in rows)
        n_expr = sum(len(r["samples_gene_expression_quantification"]) for r in rows)
        n_path = sum(len(r["samples_pathology_report"]) for r in rows)
        n_os = sum(1 for r in rows if (r.get("survival_derived") or {}).get("os_event") is not None)
        typer.echo(
            f"  {proj:<12} {len(rows):>4} patients  "
            f"mutations={n_variants:>5} ({len(mut_by_case)} MAFs)  "
            f"expression={n_expr:>4} aliquots ({len(expr_by_case)} TSVs)  "
            f"path_reports={n_path:>4}  "
            f"os={n_os}/{len(rows)}"
        )

        out = clinical.write_patients(rows, processed_dir, proj)
        typer.echo(f"             -> {out}")

    projects, gdc_releases = _card_inputs(raw_dir, processed_dir, "data.parquet")
    card = dataset_card.write_card(processed_dir, projects, gdc_releases=gdc_releases)
    typer.echo(f"wrote dataset card -> {card}")
    typer.echo(f"done. processed_patient tree at: {processed_dir}")


@app.command("build-tabular")
def build_tabular_cmd(
    project: ProjectFilterOpt = None,
    table: Annotated[
        list[str] | None,
        typer.Option(
            "--table",
            help=(
                "Table name (repeatable). If omitted, every table is built. If "
                "given, only those tables are re-derived and written; the "
                "project's other parquets are left untouched. Use this to append "
                "a newly-added modality without re-deriving the whole tree."
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

    project_files = _select_projects(raw_dir, project)

    only: set[str] | None = None
    if table:
        known = set(TABULAR_TABLES) | {
            f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS
        }
        unknown = set(table) - known
        if unknown:
            raise typer.BadParameter(
                f"unknown table(s): {sorted(unknown)}. Known: {sorted(known)}"
            )
        only = set(table)
        typer.echo(f"tables:        {', '.join(sorted(only))} (other tables left as-is)")

    # Clear prior output so removed projects/tables don't leave stale data
    # behind. Scope the wipe to what we're about to rebuild: with --project
    # only those project dirs go, leaving every other project's parquets in
    # place so a single project can be re-derived (or a new table appended)
    # without re-running — or re-uploading — the whole cohort.
    # A --table build must not wipe the tables it isn't rebuilding, so the
    # project-level wipe is skipped entirely; write_tables overwrites each
    # named table's parquet in place.
    if only is None:
        _clear_processed(processed_dir, project)

    for path in project_files:
        proj = path.parent.name
        cases = json.loads(path.read_text())
        tables = tabular.build_tables(cases, path.parent, only=only)
        # Re-derive survival endpoints (Liu et al. 2018 algorithm) — same
        # in-memory enrichment the consolidated `build` does. Lands on
        # `tables["cases"][i]["survival_derived"]`; tabular.write_tables
        # then projects it into the standalone `survival_derived` table.
        # Skipped when neither table was requested — nothing would read it.
        os_note = ""
        if "cases" in tables:
            supp_dir = path.parent / "clinical_supplement"
            supps = clinical_supplement.load_supplements_for_project(supp_dir)
            if supps:
                clinical_supplement.attach_supplements(tables["cases"], supps)
            survival.attach_survival(tables["cases"])
            n_cases = len(tables["cases"])
            n_os = sum(
                1
                for r in tables["cases"]
                if (r.get("survival_derived") or {}).get("os_event") is not None
            )
            os_note = f"  os={n_os}/{n_cases}"
            if only is None or "survival_derived" in only:
                tables["survival_derived"] = tabular.derived_survival_rows(tables["cases"])
            else:
                # `cases` was only pulled in as a dependency; don't write it.
                tables.pop("cases")

        sizes = {name: len(rows) for name, rows in tables.items()}
        # Compact one-line summary of row counts per table — easier to spot
        # cardinality regressions than scrolling through 14 lines per project.
        summary = " ".join(f"{name}={n}" for name, n in sizes.items())
        typer.echo(f"  {proj:<12} {summary}{os_note}")

        out_paths = tabular.write_tables(tables, processed_dir, proj)
        # Use any one table's path to print the project's output dir.
        # write_tables returns {} when every requested table was a
        # flex-schema one with no rows (nothing is written for those).
        if out_paths:
            typer.echo(f"             -> {next(iter(out_paths.values())).parent.parent}")

    projects, gdc_releases = _card_inputs(raw_dir, processed_dir, "*/data.parquet")

    # Tables list combines the fixed-schema TABULAR_TABLES with the
    # flex-schema clinical_supplement_* tables (one per BCR biotab form);
    # _tabular_configs_yaml filters to (project, table) pairs that
    # actually exist on disk.
    all_tables = list(TABULAR_TABLES) + [
        f"clinical_supplement_{kind}" for kind in clinical_supplement.TABULAR_FORM_KINDS
    ]
    card = dataset_card.write_tabular_card(
        processed_dir,
        projects,
        all_tables,
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
    """Push <data-dir>/processed_patient/ to a HuggingFace dataset repo."""
    root = _resolve_data_dir(data_dir)
    processed_dir = root / "processed_patient"
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
