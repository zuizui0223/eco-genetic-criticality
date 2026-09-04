from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TITLE = "Interaction thresholds and state separation under fragmentation: a theorem-guided finite-model framework"


def _flat(text: str) -> str:
    return " ".join(text.split())


def test_egc_submission_documents_are_synchronized() -> None:
    manuscript = (ROOT / "manuscript/main_text.md").read_text(encoding="utf-8")
    cover = (ROOT / "manuscript/cover_letter_theoretical_ecology.md").read_text(encoding="utf-8")
    metadata = (ROOT / "manuscript/submission_metadata.md").read_text(encoding="utf-8")
    venue = (ROOT / "manuscript/VENUE_AUDIT_2026-09-05.md").read_text(encoding="utf-8")
    display = (ROOT / "manuscript/submission_display_allocation.md").read_text(encoding="utf-8")
    assert manuscript.startswith(f"# {TITLE}\n")
    assert TITLE in cover
    assert TITLE in metadata
    assert "Primary target: Theoretical Ecology" in venue
    assert "exactly five items" in display
    assert "### Figure 1" in display
    assert "### Figure 4" in display
    assert "### Table 1" in display


def test_parent_claim_firewall_excludes_downstream_headline_results() -> None:
    manuscript = (ROOT / "manuscript/main_text.md").read_text(encoding="utf-8")
    flat = _flat(manuscript)
    for token in ("0.2543", "+5.33", "+5.20", "35/35", "48/48", "33/33", "49/49"):
        assert token not in manuscript
    assert "Predictive warning validity is assessed in the separately versioned `eco-genetic-warning-extensions` repository." in flat
    assert "state sufficiency" not in manuscript.casefold()


def test_submission_bundle_builds_from_locked_parent_evidence(tmp_path: Path) -> None:
    out = tmp_path / "egc_submission"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/build_egc_submission_bundle.py"),
            "--repo-root",
            str(ROOT),
            "--output",
            str(out),
        ],
        cwd=ROOT,
        check=True,
    )
    expected = {
        "manuscript/main_text.md",
        "manuscript/references.md",
        "manuscript/table1_evidence_and_states.md",
        "manuscript/display_allocation.md",
        "manuscript/cover_letter.md",
        "manuscript/submission_metadata.md",
        "figures/figure1_architecture.svg",
        "figures/figure2_analytical_boundaries.svg",
        "figures/figure3_fragmentation_gradient.svg",
        "figures/figure4_viability_vs_occupancy.svg",
        "tables/h3_fragmentation_gradient_pooled_summary.csv",
        "provenance/final_evidence_ledger.md",
        "MANIFEST.sha256",
    }
    observed = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}
    assert expected <= observed

    fig2 = (out / "figures/figure2_analytical_boundaries.svg").read_text(encoding="utf-8")
    fig4 = (out / "figures/figure4_viability_vs_occupancy.svg").read_text(encoding="utf-8")
    assert "K=8" in fig2
    assert "Composition mixing ≠ demographic or trait rescue" in fig2
    assert "1037/1037" in fig4
    assert "0/1037" in fig4
    assert "99.6–100%" in fig4
