from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_rach_projection_ledger_is_externally_owned():
    assert not (ROOT / "causal_model/theorem_projection_ledger.py").exists()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "zuizui0223/microdonta" in readme
    assert "sole owner of eco-genetic criticality code and evidence" in readme
