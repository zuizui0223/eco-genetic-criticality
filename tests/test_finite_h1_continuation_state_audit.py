import json
from dataclasses import replace

from causal_model.finite_h1_continuation_state_audit import (
    _hold,
    run_finite_h1_continuation_state_audit,
    write_finite_h1_continuation_state_artifacts,
)
from causal_model.finite_h1_continuation_state_audit_cli import main
from causal_model.finite_h1_hysteresis_audit import _parameters_from_terminal
from causal_model.multipatch_criticality_dynamics import simulate
from causal_model.multipatch_criticality_experiments import (
    default_scenarios,
    parameter_grid,
    parameters_for_cell,
    quick_profile,
)


def _spec():
    return replace(
        quick_profile(),
        total_area=4.0,
        patch_count=4,
        generations=1,
        replicates=1,
        master_seed=17,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
    )


def test_fresh_hold_accepts_complete_carried_terminal_state():
    spec = _spec()
    one_large = next(scenario for scenario in default_scenarios(spec) if scenario.scenario_id == "one_large")
    parameters = parameters_for_cell(spec, one_large, parameter_grid(spec)[0], seed=101)
    terminal = simulate(parameters).snapshots[-1]
    carried = _parameters_from_terminal(parameters, terminal)
    hold = _hold(carried, anchor=0.5, seed=202, generations=1, route="test")
    assert 0.0 <= hold.terminal_interaction_mean <= 1.0
    assert 0.0 <= hold.terminal_high_allele_frequency_mean <= 1.0


def test_tiny_state_audit_retains_every_seed_replicate_and_writes_artifacts(tmp_path):
    cells = run_finite_h1_continuation_state_audit(
        _spec(),
        master_seeds=(101, 202),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        hold_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )
    assert len(cells) == 1
    assert len(cells[0].replicates) == 2
    assert cells[0].summary["denominators"]["total_seed_replicates"] == 2

    csv_path = tmp_path / "summary.csv"
    json_path = tmp_path / "records.json"
    write_finite_h1_continuation_state_artifacts(cells, csv_path=csv_path, json_path=json_path)
    assert csv_path.exists() and json_path.exists()


def test_cli_writes_declared_state_transfer_policy(tmp_path):
    exit_code = main(
        [
            "--profile", "quick",
            "--output-dir", str(tmp_path),
            "--prefix", "cli",
            "--replicates", "1",
            "--generations", "1",
            "--master-seed", "101",
            "--master-seed", "202",
            "--endpoint-padding-fraction", "0.5",
            "--stage-generations", "1",
            "--hold-generations", "1",
            "--barrier-points", "5",
            "--barrier-points", "9",
            "--interaction-separation-threshold", "0.0",
            "--maximum-normalized-bracket-width", "1.0",
        ]
    )
    assert exit_code == 0
    manifest = json.loads((tmp_path / "cli.manifest.json").read_text(encoding="utf-8"))
    assert manifest["campaign"] == "finite_h1_continuation_state_hold_v1"
    assert "population" in manifest["state_transfer"]["carried_state"]
    assert "realised_trait_abundance" in manifest["state_transfer"]["carried_state"]
