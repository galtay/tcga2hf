# tcga2hf-pipeline

Build pipeline for the
[`gabrielaltay/tcga-patients-open`][patients] and
[`gabrielaltay/tcga-tabular-open`][tabular] HuggingFace datasets. Provides
the `tcga2hf-pipeline` CLI that fetches public, open-access TCGA data from
the National Cancer Institute (NCI) Genomic Data Commons (GDC), flattens
it into per-project Parquets, and pushes the result to HF Hub.

For the read-side typed models / pyarrow schemas alone, install the
companion [`tcga2hf`](../tcga2hf/) package instead — it has only `pydantic`
+ `pyarrow` as dependencies.

## Install

```bash
uv sync   # at the repo root, inside the workspace
```

(Or `pip install tcga2hf-pipeline` once published.)

## Auth

Only the `upload` and `upload-tabular` commands need auth. Put a
write-scoped HF token (from https://huggingface.co/settings/tokens) into a
local `.env` file in the repo root — gitignored, auto-loaded by the CLI:

```
HF_TOKEN=hf_xxxxxxxxxxxxxxxxxxxx
```

The CLI calls `load_dotenv(override=True)` so this project's `.env` always
wins over any inherited shell `HF_TOKEN`.

## Data location

All data lives outside the repo, under `$TCGA2HF_DATA_DIR` (default
`$HOME/data/tcga2hf`), overridable per command via `--data-dir`. The GDC
fetch step is open-access and needs no auth.

## Commands

These hit the GDC for open-access data. Acronyms used below: MAF =
Mutation Annotation Format; WXS = Whole Exome Sequencing.

```
tcga2hf-pipeline fetch-clinical    # GDC -> raw cases.json per project
tcga2hf-pipeline fetch-mutations   # GDC -> Masked Somatic Mutation MAFs from WXS (DNA)
tcga2hf-pipeline fetch-expression  # GDC -> Gene Expression Quantification TSVs from RNA-Seq (STAR counts)
tcga2hf-pipeline build             # raw -> consolidated per-project parquets + dataset card
tcga2hf-pipeline build-tabular     # raw -> per-(project, table) parquets + dataset card
tcga2hf-pipeline upload            # push processed/ to the consolidated HF dataset repo
tcga2hf-pipeline upload-tabular    # push processed_tabular/ to the tabular HF dataset repo
```

`tcga2hf-pipeline <cmd> --help` for arguments. Notable flags:

- `--project TCGA-LUAD` (repeatable) on every fetch command, plus on
  `build-tabular`. Omit on the build commands to process every project
  under `<data-dir>/raw/`.
- `--max-files N` on `fetch-mutations` / `fetch-expression`: cap the
  total number of files on disk per project (cached + freshly downloaded).
  Set to 1 to sample one file per project, or 0 to populate the manifest
  without downloading any new bytes. The manifest always lists every
  discovered file; entries trimmed off get `_status="manifest_only"`.

## Tests

Run from the workspace root:

```bash
uv run pytest -m "not network and not integration"  # offline
uv run pytest                                        # all (incl. network)
```

## Source

This package ships from the [`galtay/tcga2hf`][repo] monorepo. The
generated dataset cards live in `dataset_card.py`; build them via
`tcga2hf-pipeline build` / `build-tabular`.

[patients]: https://huggingface.co/datasets/gabrielaltay/tcga-patients-open
[tabular]: https://huggingface.co/datasets/gabrielaltay/tcga-tabular-open
[repo]: https://github.com/galtay/tcga2hf
