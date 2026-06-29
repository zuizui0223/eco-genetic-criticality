import csv
import json
from dataclasses import replace

import pytest

from causal_model.h3_branch_aware_fragmentation_audit import (
    run_h3_branch_aware_fragmentation_audit,
    write_h3_branch_aware_fragmentation_artifacts,
)
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
    scenario_one_large,
    standard_profile,
)


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=5,
        replicates=2,
        master_seed=53,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_audit_uses_one_large_branches_for_all_matched_landscapes():
    spec = _spec()
    cells = run_h3_branch_aware_fragmentation_audit(spec, interaction_separation_threshold=0.0)

    assert len(cells) == 1
    cell = cells[0]
    assert cell.canonical_one_large_h1.branch_dependent_high_trait_mode
    for replicate in cell.replicates:
        assert replicate.low_initial_interaction < replicate.high_initial_interaction
        assert replicate.outcomes is not None
        assert set(replicate.outcomes) == {
            SCENARIO_ONE_LARGE,
            SCENARIO_EQUAL_ISOLATED,
            SCENARIO_EQUAL_MIGRATING,
        }
        for scenario_id in replicate.outcomes:
            low = replicate.outcomes[scenario_id]["low_start"]
            high = replicate.outcomes[scenario_id]["high_start"]
            assert low.summary.seed == replicate.seed == high.summary.seed
            assert low.branch_id == "low_start"
            assert high.branch_id == "high_start"
            assert {warning.warning_id for warning in low.warning_orders} == {
                "H_alpha", "H_gamma", "allele_loss"
            }
        if replicate.one_large_finite_h1_mechanism_supported is True:
            assert replicate.branch_retention is not None
            assert replicate.contrasts is not None
            assert len(replicate.contrasts) == 4
        else:
            assert replicate.branch_retention is None
            assert replicate.contrasts is None

    conditioning = cell.summary["one_large_h1_conditioning"]
    assert conditioning["total_replicates"] == 2
    assert 0.0 <= conditioning["finite_h1_precondition_probability"] <= 1.0


def test_noncanonical_one_large_cells_are_unavailable_not_h3_failures():
    spec = _spec(interaction_feedback_values=(0.5,))
    cell = run_h3_branch_aware_fragmentation_audit(spec)[0]

    assert not cell.canonical_one_large_h1.branch_dependent_high_trait_mode
    assert cell.summary["one_large_h1_conditioning"]["finite_h1_precondition_count"] == 0
    assert all(replicate.outcomes is None for replicate in cell.replicates)
    assert all(replicate.one_large_finite_h1_mechanism_supported is None for replicate in cell.replicates)
    assert cell.summary["branch_retention"][SCENARIO_ONE_LARGE] is None


def test_artifacts_preserve_branch_conditioning_warnings_and_contrasts(tmp_path):
    spec = _spec(replicates=1)
    cells = run_h3_branch_aware_fragmentation_audit(spec, interaction_separation_threshold=0.0)
    csv_path = tmp_path / "h3-branch-aware.csv"
    json_path = tmp_path / "h3-branch-aware.json"
    write_h3_branch_aware_fragmentation_artifacts(cells, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(table) == len(records) == 1
    assert "one_large_h1_conditioning.finite_h1_precondition_count" in table[0]
    assert "branch_retention.equal_isolated" in table[0]
    assert (
        "branch_conditioned_warning_order.equal_migrating.high_start.H_alpha.valid_pair_count"
        in table[0]
    )
    replicate = records[0]["replicates"][0]
    assert "one_large_finite_h1_mechanism_supported" in replicate
    if replicate["outcomes"] is not None:
        assert "warning_orders" in replicate["outcomes"][SCENARIO_ONE_LARGE]["low_start"]


def test_missing_required_h3_landscapes_are_rejected():
    spec = _spec()
    with pytest.raises(ValueError, match="missing required scenarios"):
        run_h3_branch_aware_fragmentation_audit(spec, scenarios=(scenario_one_large(spec),))
