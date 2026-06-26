import pytest

from causal_model.multipatch_criticality_dynamics import DynamicsParameters
from causal_model.network_allele_floor_theory import (
    common_floor_migration_bound,
    network_allele_floor_horizon,
    network_allele_floor_one_step,
    selected_allele_floor,
)


def _supportive_parameters() -> DynamicsParameters:
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


def test_migration_preserves_a_common_selected_floor() -> None:
    assert common_floor_migration_bound(0.55, 0.0) == pytest.approx(0.55)
    assert common_floor_migration_bound(0.55, 1.0) == pytest.approx(0.55)


def test_selection_and_migration_preserve_common_floor_before_sampling() -> None:
    bound = network_allele_floor_one_step(
        _supportive_parameters(),
        patches=3,
        allele_floor=0.4,
        interaction_lower_bound=0.2,
        population_lower_bound=1000,
    )

    assert bound.selected_floor > 0.4
    assert bound.migrated_floor == pytest.approx(bound.selected_floor)
    assert bound.deterministic_floor_preserved_before_sampling
    assert bound.any_patch_sampling_failure_upper_bound < 0.01


def test_network_horizon_uses_patch_and_time_union_bound() -> None:
    certificate = network_allele_floor_horizon(
        _supportive_parameters(),
        patches=3,
        allele_floor=0.4,
        interaction_lower_bound=0.2,
        population_lower_bound=1000,
        horizon=4,
    )

    assert certificate.certified
    assert certificate.horizon_failure_upper_bound == pytest.approx(
        min(1.0, 4 * certificate.one_step.any_patch_sampling_failure_upper_bound)
    )
    assert certificate.horizon_retention_probability_lower_bound > 0.9


def test_negative_selection_margin_can_fail_the_common_floor() -> None:
    parameters = DynamicsParameters(
        patch_areas=(1.0, 1.0),
        high_base=0.2,
        high_interaction_benefit=0.0,
        viability_threshold=1.0,
        selection_strength=0.5,
    )
    bound = network_allele_floor_one_step(
        parameters,
        patches=2,
        allele_floor=0.4,
        interaction_lower_bound=0.2,
        population_lower_bound=1000,
    )
    assert bound.selected_floor < 0.4
    assert not bound.deterministic_floor_preserved_before_sampling


def test_monotonicity_conditions_are_explicit() -> None:
    with pytest.raises(ValueError):
        selected_allele_floor(
            0.4,
            0.2,
            DynamicsParameters(patch_areas=(1.0,), high_interaction_benefit=-0.1),
        )
    with pytest.raises(ValueError):
        selected_allele_floor(
            0.4,
            0.2,
            DynamicsParameters(patch_areas=(1.0,), selection_strength=-0.1),
        )
