# dev_todo

Open engineering tasks that don't fit a current iteration. Plain markdown
files, no formal tracking. Items here describe *what* and *why*; pick one up
when there's room.

## Contents

- [`optimize_downloads.md`](optimize_downloads.md) — sequential per-file
  HTTP fetches in `genomic.fetch_files` are the dominant cost when scaling
  beyond a couple of TCGA projects. Parallelize, or adopt the existing
  `gdc-client` bulk-download tool if it covers our needs.
