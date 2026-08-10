from __future__ import annotations

import json
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]


def test_release_manifest_locks_parent_scientific_state() -> None:
    manifest = json.loads(
        (ROOT / "reproducibility/release_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == 1
    assert manifest["role"] == "mechanistic_parent"
    assert manifest["scientific_commit"] == "dd8ee379d0d3518194c767d16402042525bc00dc"
    assert manifest["package"] == {
        "name": "eco-genetic-criticality",
        "version": "0.1.0",
        "python": ">=3.10",
    }
    for required in (
        "docs/final_evidence_ledger.md",
        "docs/eco_genetic_hypothesis_program.md",
        "manuscript/claim_evidence_map.md",
        "manuscript/main_text.md",
    ):
        assert required in manifest["canonical_files"]


def test_package_metadata_exposes_reproducibility_surface() -> None:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]
    assert project["name"] == "eco-genetic-criticality"
    assert "reproducibility" in project["optional-dependencies"]
    assert project["urls"]["Integrated submission"].endswith(
        "eco-genetic-warning-extensions"
    )


def test_reproducibility_guide_preserves_evidence_boundary() -> None:
    guide = (ROOT / "REPRODUCIBILITY.md").read_text(encoding="utf-8")
    assert "mechanistic parent" in guide
    assert "dd8ee379d0d3518194c767d16402042525bc00dc" in guide
    assert "must not overwrite the canonical ledger" in guide
    assert "separate provenance units" in guide
    assert "must not pool their trajectories" in guide
