"""MSigDB gene-set collections (GMT), fetched and md5-pinned.

Companion to `ssgsea.py`, which consumes these sets. Same handling as the
Liu CDR workbook in `cdr.py`: download once into `<data-dir>/raw/msigdb/`,
verify against the hash we tested with, and surface a mismatch rather than
silently consuming new bytes.

Pinning matters more here than usual. MSigDB re-releases regularly (the
collections are curated, so set membership genuinely changes between
versions), and gene-set membership feeds straight into every score. A
dataset that says "ssGSEA of Hallmark" without naming the release is not
reproducible.

Only the collections we can redistribute scores from are listed. MSigDB is
CC BY 4.0, but some constituent collections carry extra restrictions —
notably the KEGG-derived sets — so we stay with Hallmark, Reactome and
WikiPathways.

Reference: https://www.gsea-msigdb.org/gsea/msigdb
License terms: https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

import httpx

MSIGDB_VERSION = "2026.1.Hs"
_BASE_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release"


@dataclass(frozen=True)
class Collection:
    """One MSigDB collection we can fetch, pinned by md5."""

    key: str
    stem: str
    md5: str
    description: str

    @property
    def file_name(self) -> str:
        return f"{self.stem}.v{MSIGDB_VERSION}.symbols.gmt"

    @property
    def url(self) -> str:
        return f"{_BASE_URL}/{MSIGDB_VERSION}/{self.file_name}"


# md5s are of the 2026-08 snapshots we validated against. Hallmark is the
# v1 collection; Reactome is listed because it is the planned v2 addition
# and its size distribution is what makes the ssGSEA normalization
# divisor unstable (see `ssgsea.normalize_global`).
COLLECTIONS: dict[str, Collection] = {
    "hallmark": Collection(
        key="hallmark",
        stem="h.all",
        md5="367eec875967c2cfbf664a1a065b7b8d",
        description="MSigDB Hallmark — 50 coherent, deliberately non-redundant signatures",
    ),
    "reactome": Collection(
        key="reactome",
        stem="c2.cp.reactome",
        md5="1516b5d15611415d1996c92b7cb6d1cc",
        description="MSigDB C2:CP:REACTOME — 1,839 canonical pathways",
    ),
    "wikipathways": Collection(
        key="wikipathways",
        stem="c2.cp.wikipathways",
        md5="8e8a38972816a3997a557d6dd625138a",
        description="MSigDB C2:CP:WIKIPATHWAYS — 925 community-curated pathways",
    ),
}


def md5_of(path: Path, chunk: int = 1 << 16) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def fetch_collection(
    collection: str,
    raw_dir: Path,
    *,
    verify: bool = True,
    timeout: float = 120.0,
) -> Path:
    """Download one collection's GMT into `<raw_dir>/msigdb/`; return the path.

    Re-runs are cheap: an existing file whose md5 matches the pin is left
    alone. With `verify=False` the hash is reported rather than enforced,
    which is how a new collection gets its pin recorded the first time.
    """
    if collection not in COLLECTIONS:
        raise ValueError(f"unknown collection {collection!r}; known: {sorted(COLLECTIONS)}")
    spec = COLLECTIONS[collection]
    out_dir = raw_dir / "msigdb"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / spec.file_name

    if out_path.exists():
        digest = md5_of(out_path)
        if not verify or not spec.md5 or digest == spec.md5:
            return out_path
        raise ValueError(
            f"{out_path} md5 {digest} != pinned {spec.md5}. MSigDB may have re-released "
            f"{MSIGDB_VERSION}; inspect before overwriting the pin."
        )

    with httpx.Client(timeout=timeout, follow_redirects=True) as client:
        resp = client.get(spec.url)
        resp.raise_for_status()
        out_path.write_bytes(resp.content)

    digest = md5_of(out_path)
    if verify and spec.md5 and digest != spec.md5:
        out_path.unlink()
        raise ValueError(
            f"downloaded {spec.file_name} has md5 {digest}, expected {spec.md5}. "
            "Refusing to keep unverified gene sets."
        )
    return out_path
