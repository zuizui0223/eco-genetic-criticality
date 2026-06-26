import pytest

from causal_model.multipatch_criticality_dynamics import DynamicsParameters
from causal_model.stochastic_refuge_invariance import (
    RefugeRegion,
    finite_horizon_refuge_certificate,
    one_step_refuge_bound,
)


def _parameters() -> DynamicsParameters:
    return DynamicsParameters(
        patch_areas=(1.0,),
        density_capacity=1000.0,
        baseline_growth=0.8,
        interaction_growth=0.0,
        high_allele_growth=0.0,
        interaction_feedback=1.0,
        interaction_barrier=0.0,
        selection_strength=0.5,
        high_base=2.0,
        high_interaction_benefit=0.0,
        trait_grid_size=21,
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        inheritance_weight=0.5,
        migration_rate=0.0,
    )


def _region() -> RefugeRegion:
    return RefugeRegion(
        interaction_lower=0.5,
        allele_lower=0.5,
        high_trait_mass_lower=0.05,
        population_lower=800,
        population_upper=800,
    )


def test_one_step_refuge_bound_closes_declared_rectangle() -> None:
    bound = one_step_refuge_bound(_parameters(), 1.0, _region())

    assert bound.interaction_next_lower >= 0.5
    assert bound.population_next_lower >= 800
    assert bound.population_next_upper <= 800
    assert bound.selected_allele_lower > 0.5
    assert bound.high_trait_recruit_lower > 0.05
    assert bound.deterministic_region_closed
    assert 0.0 < bound.one_step_retention_probability_lower_bound <= 1.0


def test_horizon_certificate_uses_union_bound() -> None:
    certificate = finite_horizon_refuge_certificate(_parameters(), 1.0, _region(), horizon=3)

    assert certificate.certified
    assert certificate.horizon_failure_upper_bound == pytest.approx(
        min(
            1.0,
            3
            * (
                certificate.one_step.allele_failure_upper_bound
                + certificate.one_step.trait_failure_upper_bound
            ),
        )
    )


def test_certificate_rejects_migration_because_single_patch_formula_is_exact_only_without_it() -> None:
    parameters = DynamicsParameters(
        patch_areas=(1.0,),
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        migration_rate=0.1,
    )
    with pytest.raises(ValueError):
        one_step_refuge_bound(parameters, 1.0, _region())


def test_nonclosed_rectangle_is_not_certified() -> None:
    region = RefugeRegion(0.9, 0.5, 0.05, 800, 800)
    certificate = finite_horizon_refuge_certificate(_parameters(), 1.0, region, horizon=2)
    assert not certificate.certified
    assert certificate.horizon_retention_probability_lower_bound == 0.0
