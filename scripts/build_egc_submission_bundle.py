from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

from causal_model.submission_figures import (
    validate_gradient_summary,
    write_analytical_boundaries_figure,
    write_architecture_figure,
    write_viability_occupancy_figure,
)


TITLE = "Interaction thresholds and state separation under fragmentation: a theorem-guided finite-model framework"


def _copy(root: Path, out: Path, source: str, destination: str) -> None:
    src = root / source
    if not src.is_file():
        raise FileNotFoundError(src)
    dst = out / destination
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _manifest(out: Path) -> None:
    rows: list[str] = []
    for path in sorted(p for p in out.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256"):
        rows.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(out).as_posix()}")
    (out / "MANIFEST.sha256").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _validate_bundle(out: Path) -> None:
    manuscript = (out / "manuscript/main_text.md").read_text(encoding="utf-8")
    cover = (out / "manuscript/cover_letter.md").read_text(encoding="utf-8")
    if not manuscript.startswith(f"# {TITLE}\n"):
        raise RuntimeError("EGC manuscript title drifted")
    if TITLE not in cover:
        raise RuntimeError("cover letter title is not synchronized")

    required = {
        "figures/figure1_architecture.svg",
        "figures/figure2_analytical_boundaries.svg",
        "figures/figure3_fragmentation_gradient.svg",
        "figures/figure4_viability_vs_occupancy.svg",
        "manuscript/table1_evidence_and_states.md",
    }
    observed = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file()}
    if not required <= observed:
        raise RuntimeError(f"EGC submission bundle is missing main displays: {sorted(required-observed)}")

    for forbidden in ("0.2543", "+5.33", "+5.20", "35/35", "48/48", "33/33", "49/49"):
        if forbidden in manuscript:
            raise RuntimeError(f"downstream EGWE result leaked into EGC manuscript: {forbidden}")

    figure4 = (out / "figures/figure4_viability_vs_occupancy.svg").read_text(encoding="utf-8")
    for token in ("1037/1037", "0/1037", "99.6–100%"):
        if token not in figure4:
            raise RuntimeError(f"Figure 4 missing locked state-separation token: {token}")


def build_bundle(root: Path, out: Path) -> None:
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    validate_gradient_summary(root / "artifacts/h3_fragmentation_gradient/h3_fragmentation_gradient_pooled_summary.csv")

    files = {
        "manuscript/main_text.md": "manuscript/main_text.md",
        "manuscript/references.md": "manuscript/references.md",
        "manuscript/supplementary_mathematical_results.md": "supplement/supplementary_mathematical_results.md",
        "manuscript/table1_evidence_and_states.md": "manuscript/table1_evidence_and_states.md",
        "manuscript/submission_display_allocation.md": "manuscript/display_allocation.md",
        "manuscript/cover_letter_theoretical_ecology.md": "manuscript/cover_letter.md",
        "manuscript/submission_metadata.md": "manuscript/submission_metadata.md",
        "manuscript/VENUE_AUDIT_2026-09-05.md": "provenance/VENUE_AUDIT_2026-09-05.md",
        "manuscript/LITERATURE_POSITIONING_2026-09-04.md": "provenance/LITERATURE_POSITIONING_2026-09-04.md",
        "manuscript/claim_evidence_map.md": "provenance/claim_evidence_map.md",
        "docs/final_evidence_ledger.md": "provenance/final_evidence_ledger.md",
        "docs/canonical_h1_bifurcation.md": "provenance/canonical_h1_bifurcation.md",
        "docs/network_migration_matrix_theory.md": "provenance/network_migration_matrix_theory.md",
        "docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_PROTOCOL.md": "provenance/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_PROTOCOL.md",
        "docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_RESULTS.md": "provenance/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_RESULTS.md",
        "artifacts/h3_fragmentation_gradient/h3_fragmentation_gradient_pooled_summary.csv": "tables/h3_fragmentation_gradient_pooled_summary.csv",
        "artifacts/h3_fragmentation_gradient/h3_fragmentation_gradient_metadata.json": "provenance/h3_fragmentation_gradient_metadata.json",
        "artifacts/h3_fragmentation_gradient/MANIFEST.sha256": "provenance/h3_fragmentation_gradient_source_MANIFEST.sha256",
    }
    for src, dst in files.items():
        _copy(root, out, src, dst)

    write_architecture_figure(out / "figures/figure1_architecture.svg")
    write_analytical_boundaries_figure(out / "figures/figure2_analytical_boundaries.svg")
    _copy(
        root,
        out,
        "artifacts/h3_fragmentation_gradient/figure_s_fragmentation_gradient.svg",
        "figures/figure3_fragmentation_gradient.svg",
    )
    write_viability_occupancy_figure(out / "figures/figure4_viability_vs_occupancy.svg")

    _validate_bundle(out)
    _manifest(out)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build standalone EGC submission bundle")
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    build_bundle(Path(args.repo_root), Path(args.output))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
