"""End-to-end smoke test: simulator CLI -> artifact -> figure report.

This test deliberately uses the smallest valid quick-profile run.  It checks the
actual public pipeline rather than only its internal helper functions.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from causal_model.theorem_boundary_cli import build_parser, run_from_namespace
from causal_model.theorem_boundary_report import write_theorem_boundary_report


def test_quick_cli_artifact_runs_through_report_generation(tmp_path: Path):
    output_dir = tmp_path / "phase"
    args = build_parser().parse_args(
        [
            "--profile",
            "quick",
            "--scenario",
            "one_large",
            "--replicates",
            "1",
            "--generations",
            "2",
            "--master-seed",
            "20260627",
            "--output-dir",
            str(output_dir),
            "--prefix",
            "e2e",
        ]
    )

    csv_path, json_path, manifest_path = run_from_namespace(args)
    figures = write_theorem_boundary_report(csv_path, output_dir=tmp_path / "report", dpi=72)

    artifact = pd.read_csv(csv_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    replicate_cells = json.loads(json_path.read_text(encoding="utf-8"))

    assert not artifact.empty
    assert manifest["profile"] == "quick"
    assert manifest["scenario_ids"] == ["one_large"]
    assert manifest["replicate_count_per_cell"] == 1
    assert len(replicate_cells) == len(artifact)
    assert {"scenario_id", "scope.maximum_canonical_update_residual.mean"}.issubset(artifact.columns)
    assert figures
    assert all(figure.path.exists() and figure.path.stat().st_size > 0 for figure in figures)
    report = (tmp_path / "report" / "REPORT.md")
    assert report.exists()
    assert "Cells marked `NA`" in report.read_text(encoding="utf-8")
