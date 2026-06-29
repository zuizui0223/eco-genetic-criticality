import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_h1_sweep_boundary_audit import (
    run_finite_h1_sweep_boundary_audit,
    write_finite_h1_sweep_boundary_artifacts,
)
from causal_model.finite_h1_sweep_boundary_audit_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=3,
        replicates=2,
        master_seed=113,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def test_endpoint_ladder_preserves_same_seed_pairing_and_grid_limited_boundary_fields():
    cells = run_finite_h1_sweep_boundary_audit(
        _spec(),
        endpoint_padding_fractions=(0.1, 0.5),
        stage_generations=1,
        barrier_points=5,
        interaction_separation_threshold=0.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is not None
    assert len(cell.sweep_cells) == 2
    assert cell.endpoint_padding_fractions == (0.1, 0.5)
    for replicate in cell.replicates:
        assert [item.endpoint_padding_fraction for item in replicate.observations] == [0.1, 0.5]
        assert all(item.finite_loop_closed is not None for item in replicate.observations)
        assert replicate.observations[0].absolute_padding < replicate.observations[1].absolute_padding
        assert replicate.observations[0].barrier_step > 0.0
        assert isinstance(replicate.finite_loop_closed_at_any_padding, bool)
        assert isinstance(replicate.finite_h1_loop_mechanism_supported_at_any_padding, bool)

    summary = cell.summary
    assert summary["available_endpoint_sweep_probability"] == 1.0
    assert set(summary["by_endpoint_padding_fraction"]) == {"0.1", "0.5"}


def test_noncanonical_pairs_are_retained_as_unavailable():
    cells = run_finite_h1_sweep_boundary_audit(
        _spec(interaction_feedback_values=(0.5,)),
        endpoint_padding_fractions=(0.1, 0.5),
        stage_generations=1,
        barrier_points=5,
    )

    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is None
    assert cell.summary["available_endpoint_sweep_probability"] == 0.0
    assert all(replicate.finite_loop_closed_at_any_padding is None for replicate in cell.replicates)


def test_endpoint_fraction_validation_rejects_unordered_or_nonpositive_designs():
    spec = _spec()
    with pytest.raises(ValueError, match="nonempty"):
        run_finite_h1_sweep_boundary_audit(spec, endpoint_padding_fractions=())
    with pytest.raises(ValueError, match="positive"):
        run_finite_h1_sweep_boundary_audit(spec, endpoint_padding_fractions=(0.0, 0.5))
    with pytest.raises(ValueError, match="strictly increasing"):
        run_finite_h1_sweep_boundary_audit(spec, endpoint_padding_fractions=(0.5, 0.1))
    with pytest.raises(ValueError, match="at least three"):
        run_finite_h1_sweep_boundary_audit(spec, endpoint_padding_fractions=(0.1,), barrier_points=2)


def test_artifacts_and_cli_record_endpoint_design_and_no_selection_policy(tmp_path):
    cells = run_finite_h1_sweep_boundary_audit(
        _spec(replicates=1),
        endpoint_padding_fractions=(0.1, 0.5),
        stage_generations=1,
        barrier_points=5,
        interaction_separation_threshold=0.0,
    )
    csv_path = tmp_path / "sweep.csv"
    json_path = tmp_path / "sweep.json"
    write_finite_h1_sweep_boundary_artifacts(cells, csv_path=csv_path, json_path=json_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "finite_loop_closed_at_any_padding_probability" in rows[0]
    assert records[0]["endpoint_padding_fractions"] == [0.1, 0.5]

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
            "0.1",
            "--endpoint-padding-fraction",
            "0.5",
            "--stage-generations",
            "1",
            "--barrier-points",
            "5",
            "--interaction-separation-threshold",
            "0.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli.manifest.json").read_text(encoding="utf-8"))
    assert manifest["campaign"] == "finite_h1_sweep_boundary_v1"
    assert manifest["endpoint_padding_fractions"] == [0.1, 0.5]
    assert "no observed collapse" in manifest["selection_policy"]
