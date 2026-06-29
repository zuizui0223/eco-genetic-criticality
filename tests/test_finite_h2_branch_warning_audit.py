import csv
import json
from dataclasses import replace

from causal_model.finite_h2_branch_warning_audit import (
    run_finite_h2_branch_warning_audit,
    write_finite_h2_branch_warning_artifacts,
)
from causal_model.multipatch_criticality_experiments import scenario_equal_isolated, standard_profile


def _spec(**overrides):
    values = dict(
        total_area=4.0,
        patch_count=4,
        generations=6,
        replicates=3,
        master_seed=41,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(0.5,),
    )
    values.update(overrides)
    return replace(standard_profile(), **values)


def test_h2_results_are_available_only_when_same_replicate_meets_finite_h1_precondition():
    spec = _spec()
    cells = run_finite_h2_branch_warning_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        interaction_separation_threshold=0.0,
        terminal_window=2,
    )

    assert len(cells) == 1
    cell = cells[0]
    assert cell.h1_precondition["canonical_branch_dependent_high_trait_mode"]
    preconditioned = 0
    for replicate in cell.replicates:
        if replicate.finite_h1_mechanism_supported is True:
            preconditioned += 1
            assert replicate.low_start is not None
            assert replicate.high_start is not None
            assert {comparison.warning_id for comparison in replicate.low_start} == {
                "H_alpha", "H_gamma", "allele_loss"
            }
            assert {comparison.warning_id for comparison in replicate.high_start} == {
                "H_alpha", "H_gamma", "allele_loss"
            }
            assert all(comparison.branch_id == "low_start" for comparison in replicate.low_start)
            assert all(comparison.branch_id == "high_start" for comparison in replicate.high_start)
        else:
            assert replicate.low_start is None
            assert replicate.high_start is None

    conditioning = cell.summary["h1_conditioning"]
    assert conditioning["finite_h1_mechanism_precondition_count"] == preconditioned
    assert conditioning["total_replicates"] == 3
    for branch_id in ("low_start", "high_start"):
        branch = cell.summary[branch_id]
        assert branch["h1_preconditioned_replicates"] == preconditioned
        for warning in branch["warning_vs_realised_trait"].values():
            assert warning["valid_pair_count"] + warning["censored_pair_count"] == preconditioned
            if warning["lead_probability_conditional_on_valid_pair"] is not None:
                assert 0.0 <= warning["lead_probability_conditional_on_valid_pair"] <= 1.0


def test_noncanonical_cells_are_not_misreported_as_h2_warning_failures():
    spec = _spec(interaction_feedback_values=(2.0,))
    cell = run_finite_h2_branch_warning_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
    )[0]

    assert not cell.h1_precondition["canonical_branch_dependent_high_trait_mode"]
    assert cell.summary["h1_conditioning"]["finite_h1_mechanism_precondition_count"] == 0
    assert all(replicate.low_start is None and replicate.high_start is None for replicate in cell.replicates)
    for branch_id in ("low_start", "high_start"):
        for warning in cell.summary[branch_id]["warning_vs_realised_trait"].values():
            assert warning["valid_pair_count"] == 0
            assert warning["lead_probability_conditional_on_valid_pair"] is None


def test_artifacts_preserve_h1_conditioning_branch_status_and_censoring(tmp_path):
    spec = _spec(replicates=1)
    cells = run_finite_h2_branch_warning_audit(
        spec,
        scenarios=(scenario_equal_isolated(spec),),
        interaction_separation_threshold=0.0,
        terminal_window=2,
    )
    csv_path = tmp_path / "finite-h2.csv"
    json_path = tmp_path / "finite-h2.json"
    write_finite_h2_branch_warning_artifacts(cells, csv_path=csv_path, json_path=json_path)

    with csv_path.open(encoding="utf-8", newline="") as handle:
        table = tuple(csv.DictReader(handle))
    records = json.loads(json_path.read_text(encoding="utf-8"))
    assert len(table) == len(records) == 1
    assert "h1_precondition.canonical_branch_dependent_high_trait_mode" in table[0]
    assert "low_start.warning_vs_realised_trait.H_alpha.valid_pair_count" in table[0]
    assert "high_start.warning_vs_realised_trait.allele_loss.censored_pair_count" in table[0]
    replicate = records[0]["replicates"][0]
    assert "finite_h1_mechanism_supported" in replicate
    if replicate["finite_h1_mechanism_supported"] is True:
        assert replicate["low_start"][0]["censored"] in {True, False}
        assert replicate["high_start"][0]["warning_id"] == "H_alpha"
