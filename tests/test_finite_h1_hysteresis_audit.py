import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_h1_hysteresis_audit import (
    run_finite_h1_hysteresis_audit,
    write_finite_h1_hysteresis_artifacts,
)
from causal_model.multipatch_criticality_experiments import scenario_equal_isolated, standard_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=4,
        replicates=2,
        master_seed=31,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(0.5,),
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_finite_continuation_carries_rising_and_falling_paths_over_shared_barriers():
    spec = _spec()
    cells = run_finite_h1_hysteresis_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        barrier_points=5,
        stage_generations=2,
        interaction_separation_threshold=0.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_midpoint_h1.branch_dependent_high_trait_mode
    assert cell.barriers is not None
    assert cell.barriers == tuple(sorted(cell.barriers))
    assert cell.summary["finite_continuation_available_probability"] == 1.0
    for replicate in cell.replicates:
        assert replicate.rising is not None
        assert replicate.falling is not None
        assert tuple(stage.barrier for stage in replicate.rising) == cell.barriers
        assert tuple(stage.barrier for stage in replicate.falling) == tuple(reversed(cell.barriers))
        assert len({stage.stage_seed for stage in replicate.rising}) == len(cell.barriers)
        assert isinstance(replicate.finite_hysteresis_supported, bool)
        assert isinstance(replicate.finite_h1_hysteresis_mechanism_supported, bool)
        assert replicate.internal_barrier_gaps is not None
        assert all(
            cell.canonical_bistable_barrier_interval[0] < barrier < cell.canonical_bistable_barrier_interval[1]
            for barrier, _ in replicate.internal_barrier_gaps
        )

    assert 0.0 <= cell.summary["finite_hysteresis_supported_probability"] <= 1.0
    assert 0.0 <= cell.summary["finite_h1_hysteresis_mechanism_supported_probability"] <= 1.0


def test_noncanonical_cells_remain_unavailable_not_negative_hysteresis_results():
    spec = _spec(interaction_feedback_values=(2.0,))
    cell = run_finite_h1_hysteresis_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        barrier_points=5,
        stage_generations=2,
    )[0]

    assert not cell.canonical_midpoint_h1.branch_dependent_high_trait_mode
    assert cell.barriers is None
    assert cell.summary["finite_continuation_available_probability"] == 0.0
    assert cell.summary["finite_hysteresis_supported_probability"] is None
    assert all(replicate.finite_hysteresis_supported is None for replicate in cell.replicates)


def test_artifacts_keep_full_paths_and_flat_hysteresis_summary(tmp_path):
    spec = _spec(replicates=1)
    cells = run_finite_h1_hysteresis_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        barrier_points=5,
        stage_generations=2,
    )
    csv_path = tmp_path / "finite-hysteresis.csv"
    json_path = tmp_path / "finite-hysteresis.json"
    write_finite_h1_hysteresis_artifacts(cells, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(table) == len(records) == 1
    assert "finite_hysteresis_supported_probability" in table[0]
    assert "canonical_context.bistable_interval_lower" in table[0]
    record = records[0]
    assert record["barriers"]
    assert record["replicates"][0]["rising"]
    assert record["replicates"][0]["falling"]
    assert "h1_scope" in record["replicates"][0]["rising"][0]


def test_invalid_continuation_settings_are_rejected():
    spec = _spec()
    with pytest.raises(ValueError, match="at least three"):
        run_finite_h1_hysteresis_audit(spec, barrier_points=2)
    with pytest.raises(ValueError, match="positive"):
        run_finite_h1_hysteresis_audit(spec, barrier_padding=0.0)
    with pytest.raises(ValueError, match="state thresholds"):
        run_finite_h1_hysteresis_audit(spec, low_state_threshold=0.8, high_state_threshold=0.7)
