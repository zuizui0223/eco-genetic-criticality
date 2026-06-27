import json

import pytest

from causal_model.multipatch_criticality_experiments import (
    PROFILE_STANDARD,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
)
from causal_model.theorem_boundary_cli import (
    build_parser,
    profile_spec,
    run_from_namespace,
    select_scenarios,
)


def test_profile_spec_applies_explicit_runtime_overrides():
    spec = profile_spec(PROFILE_STANDARD, replicates=7, generations=9, master_seed=42)

    assert spec.profile == PROFILE_STANDARD
    assert spec.replicates == 7
    assert spec.generations == 9
    assert spec.master_seed == 42


def test_cli_runs_selected_scenario_and_writes_reproducibility_manifest(tmp_path):
    args = build_parser().parse_args(
        [
            "--profile",
            "quick",
            "--scenario",
            SCENARIO_ONE_LARGE,
            "--replicates",
            "1",
            "--generations",
            "2",
            "--master-seed",
            "7",
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "smoke",
        ]
    )

    csv_path, json_path, manifest_path = run_from_namespace(args)

    assert csv_path.exists()
    assert json_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["profile"] == "quick"
    assert manifest["scenario_ids"] == [SCENARIO_ONE_LARGE]
    assert manifest["replicate_count_per_cell"] == 1
    assert manifest["spec"]["master_seed"] == 7


def test_cli_refuses_to_overwrite_without_explicit_flag(tmp_path):
    args = build_parser().parse_args(
        [
            "--profile",
            "quick",
            "--scenario",
            SCENARIO_ONE_LARGE,
            "--replicates",
            "1",
            "--generations",
            "1",
            "--output-dir",
            str(tmp_path),
            "--prefix",
            "repeat",
        ]
    )
    run_from_namespace(args)

    with pytest.raises(FileExistsError):
        run_from_namespace(args)


def test_select_scenarios_preserves_requested_order_and_removes_duplicates():
    spec = profile_spec("quick")
    scenarios = select_scenarios(
        spec,
        (SCENARIO_EQUAL_MIGRATING, SCENARIO_ONE_LARGE, SCENARIO_EQUAL_MIGRATING),
    )

    assert [scenario.scenario_id for scenario in scenarios] == [
        SCENARIO_EQUAL_MIGRATING,
        SCENARIO_ONE_LARGE,
    ]
