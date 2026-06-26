import pytest

from causal_model.finite_bin_closure_bounds import (
    InvariantRegion,
    finite_bin_closure_bound_certificate,
    finite_bin_cohort_size_bound,
    finite_bin_trait_recruitment_bound,
    sampling_multiplier_bound,
)
from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    recruit_trait_distribution,
    trait_fitness,
    trait_grid,
)


def _parameters() -> DynamicsParameters:
    return DynamicsParameters(
        patch_areas=(1.0,),
        trait_grid_size=21,
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        inheritance_weight=0.5,
        trait_kernel_width=0.2,
    )


def _region() -> InvariantRegion:
    return InvariantRegion(
        interaction_lower_bound=0.4,
        allele_frequency_lower_bound=0.4,
        resident_high_trait_mass_lower_bound=0.2,
        population_lower_bound=10,
        population_upper_bound=20,
        next_interaction_lower_bound=0.3,
        selected_allele_frequency_lower_bound=0.2,
    )


def test_recruitment_probability_bound_is_below_declared_selected_mass() -> None:
    parameters = _parameters()
    region = _region()
    bound = finite_bin_trait_recruitment_bound(parameters, region)
    grid = trait_grid(parameters)
    high = [index for index, z in enumerate(grid) if z >= parameters.high_trait_cutoff]
    low = [index for index, z in enumerate(grid) if z < parameters.high_trait_cutoff]
    resident = [0.0] * len(grid)
    resident[low[0]] = 0.8
    resident[high[-1]] = 0.2
    recruit = recruit_trait_distribution(resident, 0.4, parameters)
    weighted = [
        mass * max(parameters.trait_selection_floor, trait_fitness(z, 0.4, parameters))
        for z, mass in zip(grid, recruit)
    ]
    selected_high_mass = sum(weighted[index] for index in high) / sum(weighted)

    assert bound.preselection_high_trait_mass_lower_bound == pytest.approx(0.3)
    assert 0.0 < bound.selected_high_trait_probability_lower_bound <= selected_high_mass


def test_cohort_bound_is_positive_and_conservative() -> None:
    bound = finite_bin_cohort_size_bound(_parameters(), 1.0, _region())
    assert bound.next_population_lower_bound >= 1
    assert bound.growth_exponent_lower_bound < 0.0


def test_sampling_bound_keeps_rho_explicit() -> None:
    bound = sampling_multiplier_bound(_parameters(), _region(), 0.9)
    assert bound.gene_copy_upper_bound >= 2
    assert bound.sampling_multiplier_upper_bound < 1.0
    assert bound.combined_diversity_multiplier_upper_bound < 1.0


def test_combined_certificate_requires_a_contractive_diversity_premise() -> None:
    good = finite_bin_closure_bound_certificate(_parameters(), 1.0, _region(), 0.9)
    not_contractive = finite_bin_closure_bound_certificate(_parameters(), 1.0, _region(), 1.2)

    assert good.l4_ready
    assert not not_contractive.l4_ready


def test_region_and_kernel_assumptions_are_checked() -> None:
    bad_region = InvariantRegion(0.3, 0.2, 0.2, 5, 4, 0.3, 0.2)
    with pytest.raises(ValueError):
        finite_bin_trait_recruitment_bound(_parameters(), bad_region)

    invalid_kernel = DynamicsParameters(
        patch_areas=(1.0,),
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        low_trait_kernel_center=0.9,
    )
    with pytest.raises(ValueError):
        finite_bin_trait_recruitment_bound(invalid_kernel, _region())
