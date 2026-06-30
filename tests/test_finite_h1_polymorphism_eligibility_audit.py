from dataclasses import replace

from causal_model.finite_h1_polymorphism_eligibility_audit import (
    _genetic_baseline,
    run_finite_h1_polymorphism_eligibility_audit,
)
from causal_model.multipatch_criticality_experiments import quick_profile


def test_preexisting_warning_is_not_h2_dynamic_warning_eligible():
    spec = quick_profile()
    baseline = _genetic_baseline(
        "equal_isolated",
        population=(32, 32, 32, 32),
        frequencies=(1.0, 1.0, 1.0, 1.0),
        spec=spec,
        epsilon=1e-12,
    )
    assert baseline.polymorphic is False
    assert baseline.h_alpha == baseline.h_gamma == 0.0
    assert baseline.h_alpha_warning_preexisting is True
    assert baseline.h_gamma_warning_preexisting is True
    assert baseline.h2_dynamic_warning_eligible is False
    assert baseline.h3_genetic_contrast_eligible is False


def test_polymorphic_baseline_above_threshold_is_eligible():
    spec = quick_profile()
    baseline = _genetic_baseline(
        "equal_isolated",
        population=(32, 32, 32, 32),
        frequencies=(0.5, 0.5, 0.5, 0.5),
        spec=spec,
        epsilon=1e-12,
    )
    assert baseline.polymorphic is True
    assert baseline.h_alpha == baseline.h_gamma == 0.5
    assert baseline.h2_dynamic_warning_eligible is True
    assert baseline.h3_genetic_contrast_eligible is True


def test_tiny_eligibility_audit_retains_all_seed_replicates():
    spec = replace(
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
    cells = run_finite_h1_polymorphism_eligibility_audit(
        spec,
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
