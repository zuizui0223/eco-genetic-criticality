import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_h1_branch_separation_audit import (
    run_finite_h1_branch_separation_audit,
    write_finite_h1_branch_separation_artifacts,
)
from causal_model.multipatch_criticality_experiments import scenario_equal_isolated, standard_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=4,
        replicates=2,
        master_seed=19,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(0.5,),
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_audit_pairs_canonical_low_and_high_starts_with_same_seed():
    spec = _spec()
    cells = run_finite_h1_branch_separation_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        interaction_separation_threshold=0.0,
        terminal_window=2,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_h1.branch_dependent_high_trait_mode
    assert cell.summary["finite_pair_available_probability"] == 1.0
    for pair in cell.replicates:
        assert pair.low_initial_interaction < pair.high_initial_interaction
        assert pair.low_start is not None
        assert pair.high_start is not None
        assert pair.low_start.seed == pair.seed == pair.high_start.seed
        assert pair.low_scope is not None
        assert pair.high_scope is not None
        assert isinstance(pair.interaction_branch_separation_supported, bool)
        assert isinstance(pair.finite_h1_mechanism_supported, bool)

    assert 0.0 <= cell.summary["interaction_branch_separation_probability"] <= 1.0
    assert 0.0 <= cell.summary["finite_h1_mechanism_supported_probability"] <= 1.0


def test_noncanonical_cells_are_retained_as_unavailable_not_false_support():
    spec = _spec(interaction_feedback_values=(2.0,))
    cell = run_finite_h1_branch_separation_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
    )[0]

    assert not cell.canonical_h1.branch_dependent_high_trait_mode
    assert cell.summary["finite_pair_available_probability"] == 0.0
    assert cell.summary["finite_h1_mechanism_supported_probability"] is None
    assert all(pair.finite_h1_mechanism_supported is None for pair in cell.replicates)


def test_artifacts_retain_canonical_context_scope_and_same_seed_records(tmp_path):
    spec = _spec(replicates=1)
    cells = run_finite_h1_branch_separation_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        interaction_separation_threshold=0.0,
    )
    csv_path = tmp_path / "finite-h1.csv"
    json_path = tmp_path / "finite-h1.json"
    write_finite_h1_branch_separation_artifacts(cells, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(table) == len(records) == 1
    assert "canonical_context.branch_dependent_high_trait_mode" in table[0]
    assert "finite_h1_mechanism_supported_probability" in table[0]
    pair = records[0]["replicates"][0]
    assert pair["low_start"]["seed"] == pair["seed"] == pair["high_start"]["seed"]
    assert "low_scope" in pair and "high_scope" in pair


def test_invalid_audit_settings_are_rejected():
    spec = _spec()
    with pytest.raises(ValueError, match="non-negative"):
        run_finite_h1_branch_separation_audit(spec, interaction_separation_threshold=-0.1)
    with pytest.raises(ValueError, match="positive"):
        run_finite_h1_branch_separation_audit(spec, terminal_window=0)
