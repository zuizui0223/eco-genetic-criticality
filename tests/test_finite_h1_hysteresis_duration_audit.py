import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_h1_hysteresis_duration_audit import (
    run_finite_h1_hysteresis_duration_audit,
    write_finite_h1_hysteresis_duration_artifacts,
)
from causal_model.finite_h1_hysteresis_duration_audit_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=3,
        replicates=2,
        master_seed=97,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def test_duration_ladder_keeps_one_seed_per_replicate_across_stage_lengths():
    cells = run_finite_h1_hysteresis_duration_audit(
        _spec(),
        stage_generations=(1, 2),
        barrier_points=3,
        interaction_separation_threshold=0.0,
        gap_stability_tolerance=10.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is not None
    assert cell.parameters.interaction_barrier == pytest.approx(
        sum(cell.canonical_bistable_barrier_interval) / 2.0
    )
    assert cell.stage_generations == (1, 2)
    assert len(cell.duration_cells) == 2
    assert cell.summary["available_duration_ladder_probability"] == 1.0
    for replicate in cell.replicates:
        assert [observation.stage_generations for observation in replicate.observations] == [1, 2]
        assert all(observation.finite_hysteresis_supported is not None for observation in replicate.observations)
        assert replicate.longest_pair_gap_change is not None
        assert replicate.longest_pair_gap_stable is True
        assert isinstance(replicate.convergence_robust_hysteresis_supported, bool)
        assert isinstance(replicate.convergence_robust_h1_mechanism_supported, bool)


def test_noncanonical_pairs_are_retained_as_unavailable_not_reported_as_negative():
    cells = run_finite_h1_hysteresis_duration_audit(
        _spec(interaction_feedback_values=(0.5,)),
        stage_generations=(1, 2),
        barrier_points=3,
    )

    cell = cells[0]
    assert cell.canonical_bistable_barrier_interval is None
    assert cell.summary["available_duration_ladder_probability"] == 0.0
    assert all(replicate.convergence_robust_hysteresis_supported is None for replicate in cell.replicates)


def test_stage_duration_ladder_validation_rejects_ambiguous_designs():
    spec = _spec()
    with pytest.raises(ValueError, match="at least two"):
        run_finite_h1_hysteresis_duration_audit(spec, stage_generations=(5,))
    with pytest.raises(ValueError, match="strictly increasing"):
        run_finite_h1_hysteresis_duration_audit(spec, stage_generations=(5, 5))
    with pytest.raises(ValueError, match="strictly increasing"):
        run_finite_h1_hysteresis_duration_audit(spec, stage_generations=(10, 5))
    with pytest.raises(ValueError, match="non-negative"):
        run_finite_h1_hysteresis_duration_audit(spec, stage_generations=(1, 2), gap_stability_tolerance=-0.1)


def test_artifacts_and_cli_freeze_duration_design_and_no_selection_policy(tmp_path):
    cells = run_finite_h1_hysteresis_duration_audit(
        _spec(replicates=1),
        stage_generations=(1, 2),
        barrier_points=3,
        interaction_separation_threshold=0.0,
    )
    csv_path = tmp_path / "duration.csv"
    json_path = tmp_path / "duration.json"
    write_finite_h1_hysteresis_duration_artifacts(cells, csv_path=csv_path, json_path=json_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "convergence_robust_hysteresis_supported_probability" in rows[0]
    assert records[0]["stage_generations"] == [1, 2]
    assert records[0]["duration_cells"][0]["stage_generations"] == 1

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
            "--stage-generations",
            "1",
            "--stage-generations",
            "2",
            "--barrier-points",
            "3",
            "--interaction-separation-threshold",
            "0.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli.manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage_generations"] == [1, 2]
    assert manifest["campaign"] == "finite_h1_hysteresis_duration_v1"
    assert "without one-large strict canonical" in manifest["selection_policy"]
