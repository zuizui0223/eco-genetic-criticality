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


def test_standalone_submission_claim_firewall() -> None:
    manifest = json.loads(
        (ROOT / "reproducibility/release_manifest.json").read_text(encoding="utf-8")
    )
    submission = manifest["standalone_submission_claims"]
    headline = set(submission["headline_claim_ids"])
    historical = set(submission["historical_non_headline_claim_ids"])

    assert headline == {"T1", "T2", "C1", "T3", "T4", "H1-S", "H3-S", "H3-G"}
    assert historical == {"H2-R", "H2-A"}
    assert headline.isdisjoint(historical)
    assert (
        submission["predictive_validity_owner"]
        == "zuizui0223/eco-genetic-warning-extensions"
    )
    assert "does not establish discrimination" in submission["claim_ceiling"]


def test_standalone_manuscript_headline_excludes_warning_claim() -> None:
    manuscript = (ROOT / "manuscript/main_text.md").read_text(encoding="utf-8")
    title = manuscript.splitlines()[0]
    abstract = manuscript.split("## Abstract", 1)[1].split("## 1. Introduction", 1)[0]
    conclusion = manuscript.split("## 7. Conclusion", 1)[1].split(
        "## Data and code availability", 1
    )[0]

    for section in (title, abstract, conclusion):
        assert "warning" not in section.casefold()


def test_h2r_remains_historical_and_non_predictive() -> None:
    claim_map = (ROOT / "manuscript/claim_evidence_map.md").read_text(
        encoding="utf-8"
    )
    ledger = (ROOT / "docs/final_evidence_ledger.md").read_text(encoding="utf-8")
    claim_map_flat = " ".join(claim_map.split())
    ledger_flat = " ".join(ledger.split())

    assert "historical event-conditioned ordering benchmark" in claim_map_flat
    assert "does not establish discrimination, specificity" in claim_map_flat
    assert "historical event-conditioned" in ledger_flat
    assert "predictive warning validity" in ledger_flat
    assert "eco-genetic-warning-extensions" in ledger_flat


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
