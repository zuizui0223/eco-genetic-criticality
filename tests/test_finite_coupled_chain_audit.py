import csv
import json
from dataclasses import replace

import pytest

from causal_model.finite_coupled_chain_audit import (
    run_finite_coupled_chain_audit,
    write_finite_coupled_chain_artifacts,
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
        generations=3,
        replicates=2,
        master_seed=23,
        area_reference_values=(1.0,),
        interaction_feedback_values=(5.0,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_audit_pairs_all_three_landscapes_by_seed_and_records_joint_predicates():
    cells = run_finite_coupled_chain_audit(_spec())

    assert len(cells) == 1
    cell = cells[0]
    assert set(cell.canonical_h1_context) == {
        SCENARIO_ONE_LARGE,
        SCENARIO_EQUAL_ISOLATED,
        SCENARIO_EQUAL_MIGRATING,
    }
    assert len(cell.replicates) == 2
    for replicate in cell.replicates:
        assert replicate.one_large.seed == replicate.seed
        assert replicate.equal_isolated.seed == replicate.seed
        assert replicate.equal_migrating.seed == replicate.seed
        assert set(replicate.h1_scope) == {
            SCENARIO_ONE_LARGE,
            SCENARIO_EQUAL_ISOLATED,
            SCENARIO_EQUAL_MIGRATING,
        }
        assert isinstance(replicate.finite_chain_supported, bool)
        assert replicate.isolated_h_alpha_leads_trait in {True, False, None}

    fragmentation = cell.summary["fragmentation_relative_to_one_large"]
    assert fragmentation["valid_isolated_H_alpha_trait_pairs"] + fragmentation[
        "censored_isolated_H_alpha_or_trait_pairs"
    ] == 2
    assert 0.0 <= fragmentation["finite_chain_supported_probability"] <= 1.0
    assert cell.summary["canonical_h1_context"][SCENARIO_ONE_LARGE]["gain"] > cell.summary[
        "canonical_h1_context"
    ][SCENARIO_EQUAL_ISOLATED]["gain"]


def test_artifacts_preserve_matched_outcomes_scope_and_censoring(tmp_path):
    cells = run_finite_coupled_chain_audit(_spec(replicates=1))
    csv_path = tmp_path / "finite-chain.csv"
    json_path = tmp_path / "finite-chain.json"
    write_finite_coupled_chain_artifacts(cells, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(table) == len(records) == 1
    assert "fragmentation_relative_to_one_large.finite_chain_supported_probability" in table[0]
    assert "h1_theorem_scope.one_large.single_patch_canonical_theorem_limit_probability" in table[0]
    replicate = records[0]["replicates"][0]
    assert set(replicate["outcomes"]) == {
        SCENARIO_ONE_LARGE,
        SCENARIO_EQUAL_ISOLATED,
        SCENARIO_EQUAL_MIGRATING,
    }
    assert "isolated_h_alpha_leads_trait" in replicate["fragmentation"]
    assert "h1_scope" in replicate


def test_audit_rejects_missing_or_unmatched_scenarios():
    spec = _spec()
    with pytest.raises(ValueError, match="missing required scenarios"):
        run_finite_coupled_chain_audit(spec, scenarios=(scenario_one_large(spec),))
