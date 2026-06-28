import json
from dataclasses import replace

import pandas as pd
import pytest

from causal_model.multipatch_criticality_experiments import scenario_one_large, standard_profile
from causal_model.paired_baseline_cli import build_parser, run_from_namespace
from causal_model.paired_baseline_comparisons import (
    BASELINE_FULL_ECO_GENETIC,
    BASELINE_GENETIC_ONLY,
    BASELINE_IDS,
    BASELINE_TRAIT_ONLY,
    baseline_definition,
    baseline_parameters,
    comparison_quick_profile,
    resolved_feedback_weights,
    run_paired_baseline_comparisons,
    write_paired_baseline_artifacts,
)


def test_ablation_replaces_removed_channel_with_interaction_memory():
    parameters = standard_profile().base_parameters

    trait_only = baseline_parameters(parameters, BASELINE_TRAIT_ONLY)
    genetic_only = baseline_parameters(parameters, BASELINE_GENETIC_ONLY)
    full = baseline_parameters(parameters, BASELINE_FULL_ECO_GENETIC)

    assert resolved_feedback_weights(full) == (0.6, 0.3, 0.1)
    assert resolved_feedback_weights(trait_only) == (0.7, 0.3, 0.0)
    assert resolved_feedback_weights(genetic_only) == (0.9, 0.0, 0.1)
    assert trait_only.genotype_trait_recruitment == "resident_trait_only"
    assert genetic_only.genotype_trait_recruitment == "resident_trait_only"
    assert trait_only.high_allele_growth == 0.0
    assert genetic_only.high_allele_growth == parameters.high_allele_growth
    assert full == parameters
    assert sum(resolved_feedback_weights(full)) == pytest.approx(1.0)
    assert sum(resolved_feedback_weights(trait_only)) == pytest.approx(1.0)
    assert sum(resolved_feedback_weights(genetic_only)) == pytest.approx(1.0)


def test_paired_comparison_shares_seed_and_records_full_model_contrasts():
    spec = replace(
        comparison_quick_profile(),
        total_area=1.0,
        patch_count=1,
        generations=2,
        replicates=2,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.4,),
        master_seed=17,
    )
    cells = run_paired_baseline_comparisons(spec, scenarios=(scenario_one_large(spec),))

    assert len(cells) == 1
    cell = cells[0]
    assert [definition.baseline_id for definition in cell.baseline_definitions] == list(BASELINE_IDS)
    assert cell.baseline_definitions[0].high_allele_growth == 0.0
    assert len(cell.replicates) == 2
    for replicate in cell.replicates:
        assert tuple(replicate.outcomes) == BASELINE_IDS
        assert all(outcome.seed == replicate.seed for outcome in replicate.outcomes.values())
    contrast = cell.summary["paired_contrasts"]["full_minus_trait_only"]
    assert contrast["replicate_count"] == 2
    assert "final_h_alpha_difference_mean" in contrast
    assert baseline_definition(spec.base_parameters, BASELINE_FULL_ECO_GENETIC).baseline_id == BASELINE_FULL_ECO_GENETIC


def test_comparison_artifacts_and_cli_manifest_are_written(tmp_path):
    spec = replace(
        comparison_quick_profile(),
        total_area=1.0,
        patch_count=1,
        generations=1,
        replicates=1,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.4,),
    )
    cells = run_paired_baseline_comparisons(spec, scenarios=(scenario_one_large(spec),))
    csv_path = tmp_path / "comparison.csv"
    json_path = tmp_path / "comparison.json"
    write_paired_baseline_artifacts(cells, csv_path=csv_path, json_path=json_path)

    flat = pd.read_csv(csv_path)
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert "models.full_eco_genetic.final_h_alpha_mean" in flat.columns
    assert "paired_contrasts.full_minus_genetic_only.final_h_alpha_difference_mean" in flat.columns
    assert records[0]["baseline_definitions"][0]["baseline_id"] == BASELINE_TRAIT_ONLY
    assert records[0]["baseline_definitions"][0]["high_allele_growth"] == 0.0

    args = build_parser().parse_args(
        [
            "--profile",
            "quick",
            "--scenario",
            "one_large",
            "--replicates",
            "1",
            "--generations",
            "1",
            "--master-seed",
            "21",
            "--output-dir",
            str(tmp_path / "cli"),
        ]
    )
    _, _, manifest_path = run_from_namespace(args)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["baseline_ids"] == list(BASELINE_IDS)
    assert manifest["baseline_definitions"][0]["high_allele_growth"] == 0.0
    assert manifest["shared_seed_rule"].startswith("Each baseline")
