from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import httpx
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

GDC_BASE_URL = "https://api.gdc.cancer.gov"


def eq(field: str, value: str) -> dict[str, Any]:
    return {"op": "=", "content": {"field": field, "value": value}}


def in_(field: str, values: Iterable[str]) -> dict[str, Any]:
    return {"op": "in", "content": {"field": field, "value": list(values)}}


def and_(*clauses: dict[str, Any]) -> dict[str, Any]:
    return {"op": "and", "content": list(clauses)}


_RETRYABLE = (httpx.TransportError, httpx.HTTPStatusError)


class GDCClient:
    def __init__(self, base_url: str = GDC_BASE_URL, timeout: float = 60.0) -> None:
        self._client = httpx.Client(base_url=base_url, timeout=timeout)

    def close(self) -> None:
        self._client.close()

    def __enter__(self) -> GDCClient:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def status(self) -> dict[str, Any]:
        """Return the GDC server status dict (includes data_release, tag, commit)."""
        resp = self._client.get("/status")
        resp.raise_for_status()
        return resp.json()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=60),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def dictionary(self) -> dict[str, Any]:
        """Return the GDC submission data dictionary (~10 MB; full snapshot).

        The dictionary has no explicit version field — its current contents
        are implicitly tied to whatever data release the GDC server is on
        (see `status()`). Capture both at the same moment to pin provenance.
        """
        resp = self._client.get("/v0/submission/_dictionary/_all", timeout=120.0)
        resp.raise_for_status()
        return resp.json()

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def _post(self, endpoint: str, payload: dict[str, Any]) -> dict[str, Any]:
        resp = self._client.post(endpoint, json=payload)
        if resp.status_code >= 500:
            resp.raise_for_status()
        resp.raise_for_status()
        return resp.json()

    def cases(
        self,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None = None,
        page_size: int = 500,
    ) -> list[dict[str, Any]]:
        """Page through /cases and return all hits (raw dicts)."""
        return self._paginate("/cases", filters, fields, expand, page_size)

    def files(
        self,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None = None,
        page_size: int = 500,
    ) -> list[dict[str, Any]]:
        """Page through /files and return all hits (raw dicts)."""
        return self._paginate("/files", filters, fields, expand, page_size)

    def _paginate(
        self,
        endpoint: str,
        filters: dict[str, Any],
        fields: list[str],
        expand: list[str] | None,
        page_size: int,
    ) -> list[dict[str, Any]]:
        hits: list[dict[str, Any]] = []
        offset = 0
        while True:
            payload: dict[str, Any] = {
                "filters": filters,
                "fields": ",".join(fields),
                "format": "JSON",
                "size": page_size,
                "from": offset,
            }
            if expand:
                payload["expand"] = ",".join(expand)
            data = self._post(endpoint, payload)
            page = data.get("data", {}).get("hits", [])
            hits.extend(page)
            pagination = data.get("data", {}).get("pagination", {})
            total = pagination.get("total", len(hits))
            if len(hits) >= total or not page:
                break
            offset += len(page)
        return hits

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        retry=retry_if_exception_type(_RETRYABLE),
    )
    def download(self, file_uuid: str, out_path: Path) -> None:
        """Stream a single open-access file by UUID. Used by future genomic modules."""
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with self._client.stream("GET", f"/data/{file_uuid}") as resp:
            resp.raise_for_status()
            with out_path.open("wb") as fh:
                for chunk in resp.iter_bytes():
                    fh.write(chunk)


def write_cases_json(cases: list[dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cases, indent=2))
