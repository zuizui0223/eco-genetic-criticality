import csv
import json
from dataclasses import replace

from causal_model.finite_validation_campaign import (
    campaign_output_paths,
    run_finite_validation_campaign,
    write_finite_validation_campaign_artifacts,
)
from causal_model.finite_validation_campaign_cli import main
from causal_model.multipatch_criticality_experiments import standard_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=3,
        replicates=1,
        master_seed=71,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(1.0,),
        migration_rate=0.2,
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_campaign_runs_all_finite_audits_on_one_shared_cell_and_ledger(tmp_path):
    spec = _spec()
    campaign = run_finite_validation_campaign(
        spec,
        interaction_separation_threshold=0.0,
        terminal_window=1,
        hysteresis_barrier_points=3,
        hysteresis_stage_generations=1,
    )

    assert campaign.scenario_ids == ("one_large", "equal_isolated", "equal_migrating")
    assert campaign.cell_count == 1
    assert len(campaign.h1_branch_cells) == 3
    assert len(campaign.h1_hysteresis_cells) == 3
    assert len(campaign.h2_warning_cells) == 3
    assert len(campaign.coupled_chain_cells) == 1
    assert len(campaign.h3_branch_aware_cells) == 1

    ledger = campaign.ledger[0]
    assert ledger["parameters"]["cell_index"] == 0
    assert set(ledger["h1_branch_by_scenario"]) == set(campaign.scenario_ids)
    assert set(ledger["h1_hysteresis_by_scenario"]) == set(campaign.scenario_ids)
    assert set(ledger["h2_warning_by_scenario"]) == set(campaign.scenario_ids)
    assert "one_large_h1_conditioning" in ledger["h3_branch_aware"]
    assert "fragmentation_relative_to_one_large" in ledger["coupled_chain"]

    paths = write_finite_validation_campaign_artifacts(
        campaign,
        output_dir=tmp_path,
        prefix="smoke",
    )
    assert all(path.exists() for path in paths.as_dict().values())
    with paths.ledger_csv.open(encoding="utf-8", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = json.loads(paths.ledger_json.read_text(encoding="utf-8"))
    assert len(rows) == len(records) == 1
    assert "h3_branch_aware.one_large_h1_conditioning.finite_h1_precondition_count" in rows[0]
    assert "coupled_chain.fragmentation_relative_to_one_large.finite_chain_supported_probability" in rows[0]


def test_campaign_output_paths_are_stable_and_complete(tmp_path):
    paths = campaign_output_paths(tmp_path, "validation")
    names = paths.as_dict()

    assert len(names) == 12
    assert names["h1_branch_csv"].name == "validation.h1_branch.csv"
    assert names["h3_branch_aware_json"].name == "validation.h3_branch_aware.json"
    assert names["ledger_json"].name == "validation.ledger.json"


def test_cli_writes_manifest_that_freezes_spec_arguments_and_selection_policy(tmp_path):
    exit_code = main(
        [
            "--profile",
            "quick",
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "cli_smoke",
            "--replicates",
            "1",
            "--generations",
            "2",
            "--terminal-window",
            "1",
            "--hysteresis-barrier-points",
            "3",
            "--hysteresis-stage-generations",
            "1",
            "--interaction-separation-threshold",
            "0.0",
        ]
    )

    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli_smoke.manifest.json").read_text(encoding="utf-8"))
    assert manifest["campaign"] == "finite_validation_campaign_v1"
    assert manifest["cell_count"] == 2
    assert manifest["replicate_count_per_cell"] == 1
    assert manifest["scenario_ids"] == ["one_large", "equal_isolated", "equal_migrating"]
    assert "censored" in manifest["selection_policy"]
    assert set(manifest["outputs"]) == {
        "h1_branch_csv",
        "h1_branch_json",
        "h1_hysteresis_csv",
        "h1_hysteresis_json",
        "h2_warning_csv",
        "h2_warning_json",
        "coupled_chain_csv",
        "coupled_chain_json",
        "h3_branch_aware_csv",
        "h3_branch_aware_json",
        "ledger_csv",
        "ledger_json",
    }
