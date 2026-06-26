import pytest

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, sigmoid, simulate
from causal_model.multipatch_criticality_ensemble import (
    canonical_interaction_update,
    simulate_ensemble,
)


def test_canonical_reduction_matches_declared_logistic_map():
    q = 0.4
    area = 1.5
    area_reference = 1.0
    feedback = 3.0
    barrier = 0.2
    expected = sigmoid(feedback * ((area / area_reference) * q - barrier))
    assert canonical_interaction_update(
        q,
        area=area,
        area_reference=area_reference,
        feedback_strength=feedback,
        barrier=barrier,
    ) == pytest.approx(expected)


def test_simulator_first_update_has_canonical_reduction_when_density_is_one():
    area = 1.5
    q = 0.4
    feedback = 3.0
    barrier = 0.2
    density_capacity = 40.0
    parameters = DynamicsParameters(
        patch_areas=(area,),
        initial_population=(round(density_capacity * area),),
        initial_interaction=(q,),
        initial_high_allele_frequency=(0.1,),
        density_capacity=density_capacity,
        interaction_memory_weight=1.0,
        interaction_feedback=feedback,
        interaction_barrier=barrier,
        generations=1,
        random_seed=1,
    )
    observed = simulate(parameters).snapshots[1].interaction[0]
    expected = canonical_interaction_update(
        q,
        area=area,
        area_reference=parameters.area_reference,
        feedback_strength=feedback,
        barrier=barrier,
    )
    assert observed == pytest.approx(expected)


def test_ensemble_reports_mean_path_and_lead_probability_separately():
    parameters = DynamicsParameters(
        patch_areas=(1.0, 1.0),
        initial_population=(20, 20),
        initial_interaction=(0.4, 0.4),
        initial_high_allele_frequency=(0.5, 0.5),
        generations=5,
        random_seed=11,
    )
    results, summary = simulate_ensemble(parameters, replicates=4, warning_threshold=0.1)
    assert len(results) == 4
    assert len(summary.mean_h_alpha) == 6
    assert len(summary.mean_h_gamma) == 6
    assert 0.0 <= summary.genetic_lead_probability <= 1.0
    assert len(summary.trait_absence_times) == 4
    assert len(summary.alpha_warning_times) == 4


def test_ensemble_is_reproducible_and_does_not_treat_missing_event_as_lead():
    parameters = DynamicsParameters(patch_areas=(1.0,), generations=3, random_seed=7)
    _, first = simulate_ensemble(parameters, replicates=3, warning_threshold=0.0)
    _, second = simulate_ensemble(parameters, replicates=3, warning_threshold=0.0)
    assert first == second
    assert first.genetic_lead_probability == 0.0


def test_invalid_ensemble_arguments_are_rejected():
    parameters = DynamicsParameters(patch_areas=(1.0,))
    with pytest.raises(ValueError):
        simulate_ensemble(parameters, replicates=0, warning_threshold=0.2)
    with pytest.raises(ValueError):
        simulate_ensemble(parameters, replicates=1, warning_threshold=1.1)
