import pytest

from causal_model.multipatch_criticality_dynamics import DynamicsParameters
from causal_model.network_high_trait_floor_theory import (
    network_high_trait_floor_horizon,
    network_high_trait_floor_one_step,
    preselection_high_trait_mass_lower_bound,
    selected_high_trait_probability_lower_bound,
)


def _parameters() -> DynamicsParameters:
    return DynamicsParameters(
        patch_areas=(1.0, 1.0, 1.0),
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        inheritance_weight=0.5,
        high_base=2.0,
        high_interaction_benefit=0.0,
        high_peak_width=0.4,
        viability_threshold=1.0,
        trait_grid_size=21,
        realised_high_trait_abundance_threshold=1,
    )


def test_two_kernel_preselection_mass_has_exact_mixture_lower_envelope() -> None:
    mass = preselection_high_trait_mass_lower_bound(0.4, 0.2, _parameters())
    assert mass == pytest.approx(0.3)


def test_selected_high_trait_probability_exceeds_floor_in_supportive_region() -> None:
    probability = selected_high_trait_probability_lower_bound(
        0.4,
        0.1,
        0.2,
        _parameters(),
    )
    assert probability > 0.1


def test_network_common_high_trait_floor_is_certified() -> None:
    bound = network_high_trait_floor_one_step(
        _parameters(),
        patches=3,
        allele_floor=0.4,
        resident_high_trait_mass_floor=0.05,
        interaction_lower_bound=0.2,
        cohort_lower_bound=1000,
    )

    assert bound.preselection_high_trait_mass_lower_bound == pytest.approx(0.225)
    assert bound.selected_high_trait_probability_lower_bound > 0.05
    assert bound.deterministic_mass_floor_possible
    assert bound.any_patch_failure_upper_bound < 0.01


def test_network_high_trait_horizon_uses_patch_and_time_union_bound() -> None:
    certificate = network_high_trait_floor_horizon(
        _parameters(),
        patches=3,
        allele_floor=0.4,
        resident_high_trait_mass_floor=0.05,
        interaction_lower_bound=0.2,
        cohort_lower_bound=1000,
        horizon=4,
    )

    assert certificate.certified
    assert certificate.horizon_failure_upper_bound == pytest.approx(
        min(1.0, 4 * certificate.one_step.any_patch_failure_upper_bound)
    )
    assert certificate.horizon_retention_probability_lower_bound > 0.9


def test_mass_floor_must_imply_count_based_occupancy() -> None:
    with pytest.raises(ValueError):
        network_high_trait_floor_one_step(
            _parameters(),
            patches=2,
            allele_floor=0.4,
            resident_high_trait_mass_floor=0.001,
            interaction_lower_bound=0.2,
            cohort_lower_bound=100,
        )


def test_theorem_rejects_nonpartitioned_kernel_declaration() -> None:
    bad = DynamicsParameters(
        patch_areas=(1.0,),
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        low_trait_kernel_center=0.9,
    )
    with pytest.raises(ValueError):
        preselection_high_trait_mass_lower_bound(0.4, 0.1, bad)
