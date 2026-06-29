import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_h1_boundary_seed_ensemble import (
    run_finite_h1_boundary_seed_ensemble,
    write_finite_h1_boundary_seed_ensemble_artifacts,
)
from causal_model.finite_h1_boundary_seed_ensemble_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=1,
        replicates=1,
        master_seed=11,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def test_seed_ensemble_aligns_parameter_cells_and_preserves_per_seed_raw_runs():
    cells = run_finite_h1_boundary_seed_ensemble(
        _spec(),
        master_seeds=(101, 202),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert [run.master_seed for run in cell.seed_runs] == [101, 202]
    assert all(len(run.cell.replicates) == 1 for run in cell.seed_runs)
    summary = cell.summary
    assert summary["master_seed_count"] == 2
    assert summary["replicates_per_master_seed"] == 1
    assert summary["total_replicates"] == 2
    assert set(summary["by_master_seed"]) == {"101", "202"}
    assert summary["available_master_seed_probability"] == 1.0
    assert summary["pooled_available_replicate_probability"] == 1.0
    assert summary["pooled_resolution_stable_loop_supported_probability"] in {0.0, 1.0}
    assert isinstance(summary["all_master_seed_runs_fully_resolution_stable"], bool)
    assert isinstance(summary["all_master_seed_runs_fully_resolution_stable_h1_mechanism"], bool)


def test_seed_ensemble_rejects_single_duplicate_and_negative_seed_designs():
    spec = _spec()
    kwargs = dict(
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
    )
    with pytest.raises(ValueError, match="at least two"):
        run_finite_h1_boundary_seed_ensemble(spec, master_seeds=(101,), **kwargs)
    with pytest.raises(ValueError, match="distinct"):
        run_finite_h1_boundary_seed_ensemble(spec, master_seeds=(101, 101), **kwargs)
    with pytest.raises(ValueError, match="non-negative"):
        run_finite_h1_boundary_seed_ensemble(spec, master_seeds=(101, -1), **kwargs)


def test_artifacts_and_cli_freeze_seed_ensemble_without_outcome_selection(tmp_path):
    cells = run_finite_h1_boundary_seed_ensemble(
        _spec(),
        master_seeds=(101, 202),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )
    csv_path = tmp_path / "ensemble.csv"
    json_path = tmp_path / "ensemble.json"
    write_finite_h1_boundary_seed_ensemble_artifacts(cells, csv_path=csv_path, json_path=json_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "all_master_seed_runs_fully_resolution_stable" in rows[0]
    assert records[0]["master_seed_count"] == 2
    assert records[0]["seed_runs"][0]["master_seed"] == 101

    exit_code = main(
        [
            "--profile",
            "quick",
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "cli",
            "--replicates",
            "1",
            "--generations",
            "1",
            "--ensemble-master-seed",
            "101",
            "--ensemble-master-seed",
            "202",
            "--endpoint-padding-fraction",
            "0.5",
            "--stage-generations",
            "1",
            "--barrier-points",
            "5",
            "--barrier-points",
            "9",
            "--interaction-separation-threshold",
            "0.0",
            "--maximum-normalized-bracket-width",
            "1.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli.manifest.json").read_text(encoding="utf-8"))
    assert manifest["campaign"] == "finite_h1_boundary_seed_ensemble_v1"
    assert manifest["ensemble_master_seeds"] == [101, 202]
    assert "No seed-specific boundary" in manifest["selection_policy"]
