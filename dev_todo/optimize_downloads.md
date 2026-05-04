# Optimize per-file downloads

## Why this matters

`tcga2hf.genomic.fetch_files` downloads files one at a time, sequentially:

```python
for hit in hits:
    client.download(file_id, target)
```

Each download is dominated by HTTP round-trip latency, not bytes. For small
files (MAFs, ~30-100 KB each) we're losing most of the wall time to
connection setup and TLS handshakes per file.

## Observed cost

CHOL+DLBC fetch (98 + 95 files = 193 files):
- mutations (~5 MB total): ~30 sec wall, but the bytes alone would transfer
  in <1 sec at typical residential bandwidth.
- expression (~600 MB total): ~30 sec wall, mostly bytes.

Estimated cost at full TCGA scale (~33 projects, ~11k cases, ~10k MAFs +
~10k expression TSVs ≈ ~50 GB total):
- Sequential: hours to a full day, mostly HTTP overhead.
- Parallel (8-16 workers): bandwidth-limited, probably 1-2 hours.

## Investigation before implementation

Before writing our own threadpool / async downloader, check what already
exists. The GDC ships an official client:

- **gdc-client** — official CLI for bulk downloads, takes a manifest TSV
  and pulls in parallel with retries. Documentation:
  <https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/>
- It can be invoked from Python via subprocess, or we can use the manifest
  format directly.

Questions to answer:
1. Does `gdc-client` work on open-access data without an auth token? (For
   controlled-access we'd need a token; we don't fetch that, so probably yes.)
2. How does it handle retries / resume on partial failure?
3. Is its output layout something we can adapt our manifest pattern to,
   or does it impose a directory structure we'd have to work around?
4. Is there a Python library wrapper, or only the standalone binary?

If `gdc-client` is a clean fit, the simplest path is:
1. Build a GDC manifest TSV from the file list our `/files` query produces.
2. Shell out to `gdc-client download -m <manifest>` per project.
3. Generate our `manifest.json` (case/sample mapping) separately, since
   gdc-client doesn't carry that linkage.

If it isn't a fit, fall back to a small `httpx.AsyncClient` or
`concurrent.futures.ThreadPoolExecutor` (~8-16 workers) wrapping our
existing `client.download`. Keep tenacity retries; respect `Retry-After`
headers if GDC pushes back.

## Acceptance signals

- Full TCGA fetch (any single large project) completes in minutes, not hours.
- Re-run is cached / idempotent (same as today).
- Failure mid-fetch can be resumed without re-downloading completed files.

## Out of scope here

- Parallelizing the parquet build itself (parsing TSVs in parallel). That's
  a separate optimization once download is no longer the bottleneck.
