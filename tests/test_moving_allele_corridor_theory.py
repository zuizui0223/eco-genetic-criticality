import pytest

from causal_model.moving_allele_corridor_theory import (
    moving_corridor_certificate,
    moving_corridor_step,
    selected_frequency,
    selection_shifted_corridor,
)
from causal_model.multipatch_criticality_dynamics import DynamicsParameters


def _parameters() -> DynamicsParameters:
    return DynamicsParameters(
        patch_areas=(1.0, 1.0, 1.0),
        high_base=2.0,
        high_interaction_benefit=0.0,
        viability_threshold=1.0,
        selection_strength=0.5,
        effective_fraction=1.0,
        skew_penalty=0.0,
        migration_rate=0.8,
    )


def test_positive_selection_moves_any_nontrivial_fixed_upper_endpoint_upward() -> None:
    assert selected_frequency(0.7, 1.5) > 0.7
    step = moving_corridor_step(
        _parameters(),
        patches=2,
        lower_before=0.6,
        upper_before=0.7,
        lower_after=0.6,
        upper_after=0.7,
        interaction_lower_bound=0.0,
        interaction_upper_bound=0.0,
        population_lower_bound=1000,
    )
    assert not step.deterministic_corridor_contains_envelope
    assert step.upper_escape_probability_upper_bound == 1.0


def test_selection_shifted_corridor_contains_deterministic_envelope() -> None:
    lower, upper = selection_shifted_corridor(
        _parameters(),
        lower_initial=0.60,
        upper_initial=0.70,
        interaction_lower_bound=0.0,
        interaction_upper_bound=0.0,
        generations=3,
        slack=0.05,
    )
    assert len(lower) == len(upper) == 4
    assert lower[1] > lower[0]
    assert upper[1] > upper[0]

    certificate = moving_corridor_certificate(
        _parameters(),
        patches=3,
        lower_path=lower,
        upper_path=upper,
        interaction_lower_bound=0.0,
        interaction_upper_bound=0.0,
        population_lower_bound=1000,
    )
    assert certificate.certified
    assert certificate.horizon_retention_probability_lower_bound > 0.9


def test_corridor_requires_strict_sampling_room() -> None:
    with pytest.raises(ValueError):
        moving_corridor_certificate(
            _parameters(),
            patches=2,
            lower_path=(0.6,),
            upper_path=(0.7,),
            interaction_lower_bound=0.0,
            interaction_upper_bound=0.0,
            population_lower_bound=1000,
        )


def test_invalid_interval_is_rejected() -> None:
    with pytest.raises(ValueError):
        moving_corridor_step(
            _parameters(),
            patches=2,
            lower_before=0.8,
            upper_before=0.7,
            lower_after=0.8,
            upper_after=0.9,
            interaction_lower_bound=0.0,
            interaction_upper_bound=0.0,
            population_lower_bound=100,
        )
