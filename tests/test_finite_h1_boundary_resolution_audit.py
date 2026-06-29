import csv
import json
from dataclasses import replace
from types import SimpleNamespace

import pytest

from causal_model.finite_h1_boundary_resolution_audit import (
    _falling_bracket,
    _rising_bracket,
    run_finite_h1_boundary_resolution_audit,
    write_finite_h1_boundary_resolution_artifacts,
)
from causal_model.finite_h1_boundary_resolution_audit_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=2,
        replicates=1,
        master_seed=127,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def _stage(barrier, interaction):
    return SimpleNamespace(barrier=barrier, terminal_interaction_mean=interaction)


def test_transition_brackets_store_adjacent_grid_cells_not_only_crossing_points():
    rising = _rising_bracket(
        (_stage(0.0, 0.9), _stage(0.5, 0.6), _stage(1.0, 0.1)),
        low_state_threshold=0.25,
    )
    falling = _falling_bracket(
        (_stage(1.0, 0.1), _stage(0.5, 0.3), _stage(0.0, 0.8)),
        high_state_threshold=0.75,
    )

    assert rising is not None
    assert rising.lower_barrier == 0.5
    assert rising.upper_barrier == 1.0
    assert rising.crossing_barrier == 1.0
    assert rising.width == 0.5
    assert falling is not None
    assert falling.lower_barrier == 0.0
    assert falling.upper_barrier == 0.5
    assert falling.crossing_barrier == 0.0
    assert falling.width == 0.5


def test_nested_grid_audit_keeps_same_seed_and_records_resolution_specific_observations():
    cells = run_finite_h1_boundary_resolution_audit(
        _spec(),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is not None
    assert cell.nested_barrier_points == (5, 9)
    assert len(cell.resolution_cells) == 2
    assert cell.summary["available_nested_grid_probability"] == 1.0
    replicate = cell.replicates[0]
    assert [item.barrier_points for item in replicate.observations] == [5, 9]
    assert all(item.finite_loop_bracket_supported is not None for item in replicate.observations)
    assert isinstance(replicate.loop_on_two_finest_grids, bool)
    assert isinstance(replicate.resolution_stable_loop_supported, bool)


def test_noncanonical_pairs_are_retained_as_unavailable():
    cells = run_finite_h1_boundary_resolution_audit(
        _spec(interaction_feedback_values=(0.5,)),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
    )

    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is None
    assert cell.summary["available_nested_grid_probability"] == 0.0
    assert cell.replicates[0].resolution_stable_loop_supported is None


def test_nested_grid_validation_rejects_non_nested_or_ambiguous_designs():
    spec = _spec()
    with pytest.raises(ValueError, match="at least two"):
        run_finite_h1_boundary_resolution_audit(spec, nested_barrier_points=(5,))
    with pytest.raises(ValueError, match="strictly increasing"):
        run_finite_h1_boundary_resolution_audit(spec, nested_barrier_points=(9, 5))
    with pytest.raises(ValueError, match="divisible"):
        run_finite_h1_boundary_resolution_audit(spec, nested_barrier_points=(5, 8))
    with pytest.raises(ValueError, match="positive"):
        run_finite_h1_boundary_resolution_audit(spec, endpoint_padding_fraction=0.0)


def test_artifacts_and_cli_freeze_common_padding_and_nested_grid_design(tmp_path):
    cells = run_finite_h1_boundary_resolution_audit(
        _spec(),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )
    csv_path = tmp_path / "resolution.csv"
    json_path = tmp_path / "resolution.json"
    write_finite_h1_boundary_resolution_artifacts(cells, csv_path=csv_path, json_path=json_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "resolution_stable_loop_supported_probability" in rows[0]
    assert records[0]["nested_barrier_points"] == [5, 9]

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
    assert manifest["campaign"] == "finite_h1_boundary_resolution_v1"
    assert manifest["endpoint_padding_fraction"] == 0.5
    assert manifest["nested_barrier_points"] == [5, 9]
    assert "No pair-specific observed boundary" in manifest["selection_policy"]
