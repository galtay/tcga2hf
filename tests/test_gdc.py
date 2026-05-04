from __future__ import annotations

import pytest

from tcga2hf.gdc import GDCClient, and_, eq, in_


def test_filter_helpers() -> None:
    assert eq("a", "x") == {"op": "=", "content": {"field": "a", "value": "x"}}
    assert in_("a", ["x", "y"]) == {"op": "in", "content": {"field": "a", "value": ["x", "y"]}}
    combined = and_(eq("a", "x"), eq("b", "y"))
    assert combined["op"] == "and"
    assert len(combined["content"]) == 2


@pytest.mark.network
def test_cases_smoke() -> None:
    """Live GDC API smoke test. Skip with `pytest -m 'not network'`."""
    with GDCClient() as c:
        cases = c.cases(
            filters=eq("project.project_id", "TCGA-CHOL"),
            fields=["case_id", "submitter_id"],
            page_size=3,
        )
    assert len(cases) >= 3
    assert all("case_id" in case for case in cases)
