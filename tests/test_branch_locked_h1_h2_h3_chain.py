import csv
import json
from dataclasses import replace

import pytest

from causal_model.branch_locked_h1_h2_h3_chain import (
    _validate_master_seeds,
    run_branch_locked_h1_h2_h3_chain_audit,
    write_branch_locked_h1_h2_h3_chain_artifacts,
)
from causal_model.branch_locked_h1_h2_h3_chain_cli import main
from causal_model.multipatch_criticality_experiments import quick_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=1,
        replicates=1,
        master_seed=17,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(quick_profile(), **values)


def test_chain_audit_retains_every_seed_replicate_and_separates_denominators():
    cells = run_branch_locked_h1_h2_h3_chain_audit(
        _spec(),
        master_seeds=(101, 202),
        h1_endpoint_padding_fraction=0.5,
        h1_stage_generations=1,
        h1_nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        h1_maximum_normalized_bracket_width=1.0,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.master_seeds == (101, 202)
    assert len(cell.replicates) == 2
    assert {record.master_seed for record in cell.replicates} == {101, 202}
    assert {record.replicate_index for record in cell.replicates} == {0}
    assert cell.summary["denominators"]["total_seed_replicates"] == 2
    assert "h2_isolated_high" in cell.summary
    assert "h3_high_branch_fragmentation" in cell.summary
    assert "migration_as_allele_frequency_mixing" in cell.summary
    assert "same_replicate_chain" in cell.summary
    assert set(cell.summary["by_master_seed"]) == {"101", "202"}


def test_master_seed_validation_rejects_single_duplicate_and_negative_designs():
    with pytest.raises(ValueError, match="at least two"):
        _validate_master_seeds((101,))
    with pytest.raises(ValueError, match="distinct"):
        _validate_master_seeds((101, 101))
    with pytest.raises(ValueError, match="non-negative"):
        _validate_master_seeds((101, -1))


def test_artifacts_and_cli_freeze_h1_h2_h3_conditions_without_outcome_selection(tmp_path):
    cells = run_branch_locked_h1_h2_h3_chain_audit(
        _spec(),
        master_seeds=(101, 202),
        h1_endpoint_padding_fraction=0.5,
        h1_stage_generations=1,
        h1_nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        h1_maximum_normalized_bracket_width=1.0,
    )
    csv_path = tmp_path / "chain.csv"
    json_path = tmp_path / "chain.json"
    write_branch_locked_h1_h2_h3_chain_artifacts(cells, csv_path=csv_path, json_path=json_path)
    with csv_path.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "same_replicate_chain.h_alpha.support_probability_across_all_seed_replicates" in rows[0]
    assert records[0]["master_seeds"] == [101, 202]
    assert len(records[0]["replicates"]) == 2

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
            "--master-seed",
            "101",
            "--master-seed",
            "202",
            "--h1-endpoint-padding-fraction",
            "0.5",
            "--h1-stage-generations",
            "1",
            "--h1-barrier-points",
            "5",
            "--h1-barrier-points",
            "9",
            "--interaction-separation-threshold",
            "0.0",
            "--h1-maximum-normalized-bracket-width",
            "1.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli.manifest.json").read_text(encoding="utf-8"))
    assert manifest["campaign"] == "branch_locked_h1_h2_h3_chain_v1"
    assert manifest["master_seeds"] == [101, 202]
    assert manifest["h2_primary_endpoint"]["landscape"] == "equal_isolated"
    assert "No H2/H3 outcome selects records" in manifest["selection_policy"]
