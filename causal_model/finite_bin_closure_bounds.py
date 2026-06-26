"""Analytic bound ingredients for the declared finite-bin simulator closure.

This module derives only those L4 inputs that follow algebraically from the
finite-bin update equations once an invariant state region is declared. It does
not estimate bounds from simulated trajectories.

For one designated refuge patch, assume the state remains in a region with
lower bounds on interaction, allele frequency, resident high-trait mass, and
population, plus an upper bound on population. The module derives:

* a lower bound on selected high-trait recruitment probability ``pi_min``;
* a lower bound on the next finite recruitment cohort size ``n_min``;
* the exact Wright--Fisher sampling contraction factor conditional on a bounded
  effective population size.

Selection and migration can change H-alpha before allele sampling. Therefore a
pre-sampling diversity expansion factor ``rho`` remains an explicit external
assumption. This is deliberate: it prevents a simulation observation from being
silently promoted to a theorem premise.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, floor

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    trait_fitness,
    trait_grid,
)


@dataclass(frozen=True)
class InvariantRegion:
    """Declared one-patch bounds over a time interval.

    ``interaction_lower_bound`` and ``allele_frequency_lower_bound`` refer to
    the current state used by finite trait recruitment. ``population_lower`` and
    ``population_upper`` bound the current census size. These are assumptions
    that must be established separately for a chosen parameter region.
    """

    interaction_lower_bound: float
    allele_frequency_lower_bound: float
    resident_high_trait_mass_lower_bound: float
    population_lower_bound: int
    population_upper_bound: int
    next_interaction_lower_bound: float
    selected_allele_frequency_lower_bound: float


@dataclass(frozen=True)
class TraitRecruitmentBound:
    """Selected high-trait recruit probability implied by the declared region."""

    preselection_high_trait_mass_lower_bound: float
    high_trait_fitness_lower_bound: float
    total_fitness_upper_bound: float
    selected_high_trait_probability_lower_bound: float


@dataclass(frozen=True)
class CohortSizeBound:
    """Lower bound on one next-generation recruitment cohort size."""

    growth_exponent_lower_bound: float
    next_population_lower_bound: int


@dataclass(frozen=True)
class SamplingMultiplierBound:
    """Wright--Fisher sampling and externally declared pre-sampling expansion."""

    effective_size_upper_bound: float
    gene_copy_upper_bound: int
    sampling_multiplier_upper_bound: float
    pre_sampling_diversity_expansion_upper_bound: float
    combined_diversity_multiplier_upper_bound: float


@dataclass(frozen=True)
class FiniteBinClosureBoundCertificate:
    """L4-ready bounds derived from closure plus explicit invariant assumptions."""

    region: InvariantRegion
    trait_recruitment: TraitRecruitmentBound
    cohort_size: CohortSizeBound
    sampling: SamplingMultiplierBound
    l4_ready: bool


def _validate_probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def validate_invariant_region(region: InvariantRegion) -> InvariantRegion:
    """Validate a declared region without claiming that it is invariant."""
    q = _validate_probability(region.interaction_lower_bound, "interaction_lower_bound")
    p = _validate_probability(region.allele_frequency_lower_bound, "allele_frequency_lower_bound")
    r = _validate_probability(
        region.resident_high_trait_mass_lower_bound,
        "resident_high_trait_mass_lower_bound",
    )
    q_next = _validate_probability(region.next_interaction_lower_bound, "next_interaction_lower_bound")
    p_selected = _validate_probability(
        region.selected_allele_frequency_lower_bound,
        "selected_allele_frequency_lower_bound",
    )
    if region.population_lower_bound < 1:
        raise ValueError("population_lower_bound must be at least one")
    if region.population_upper_bound < region.population_lower_bound:
        raise ValueError("population_upper_bound must be at least population_lower_bound")
    return InvariantRegion(q, p, r, int(region.population_lower_bound), int(region.population_upper_bound), q_next, p_selected)


def finite_bin_trait_recruitment_bound(
    parameters: DynamicsParameters,
    region: InvariantRegion,
) -> TraitRecruitmentBound:
    """Derive a high-trait recruit probability bound for ``two_kernel_recruitment``.

    The derivation requires the declared low and high kernels to be supported on
    opposite sides of ``high_trait_cutoff``. In that case the high-region mass
    before viability selection is at least

        (1-w) p_min + w r_min,

    where ``w`` is ``inheritance_weight``, ``p_min`` the allele-frequency lower
    bound, and ``r_min`` the resident high-trait mass lower bound. Viability
    selection then gives ``pi_min >= r_pre W_H,min / W_max``.
    """
    region = validate_invariant_region(region)
    if parameters.trait_occupancy_mode != "finite_trait_bin_recruitment":
        raise ValueError("finite-bin recruitment bound requires finite_trait_bin_recruitment mode")
    if parameters.genotype_trait_recruitment != "two_kernel_recruitment":
        raise ValueError("bound requires two_kernel_recruitment")
    if parameters.low_trait_kernel_center >= parameters.high_trait_cutoff:
        raise ValueError("low trait kernel must be centred below high_trait_cutoff")
    if parameters.high_trait_kernel_center < parameters.high_trait_cutoff:
        raise ValueError("high trait kernel must be centred in the high-trait region")

    pre_high = (
        (1.0 - parameters.inheritance_weight) * region.allele_frequency_lower_bound
        + parameters.inheritance_weight * region.resident_high_trait_mass_lower_bound
    )
    grid = trait_grid(parameters)
    high_fitness_lower = min(
        max(parameters.trait_selection_floor, trait_fitness(z, region.interaction_lower_bound, parameters))
        for z in grid
        if z >= parameters.high_trait_cutoff
    )
    total_fitness_upper = max(
        max(parameters.trait_selection_floor, trait_fitness(z, 1.0, parameters))
        for z in grid
    )
    selected_lower = min(1.0, pre_high * high_fitness_lower / total_fitness_upper)
    return TraitRecruitmentBound(pre_high, high_fitness_lower, total_fitness_upper, selected_lower)


def finite_bin_cohort_size_bound(
    parameters: DynamicsParameters,
    patch_area: float,
    region: InvariantRegion,
) -> CohortSizeBound:
    """Derive a safe lower bound for the simulator's rounded population update."""
    region = validate_invariant_region(region)
    if patch_area <= 0.0:
        raise ValueError("patch_area must be positive")
    carrying_capacity = parameters.density_capacity * patch_area
    exponent_lower = (
        parameters.baseline_growth
        + parameters.interaction_growth * region.next_interaction_lower_bound
        + parameters.high_allele_growth * region.selected_allele_frequency_lower_bound
        - region.population_upper_bound / carrying_capacity
    )
    raw_lower = region.population_lower_bound * exp(exponent_lower)
    # Python's bankers rounding is always at least ceil(x - 1/2).
    rounded_lower = ceil(raw_lower - 0.5)
    return CohortSizeBound(exponent_lower, max(1, rounded_lower))


def sampling_multiplier_bound(
    parameters: DynamicsParameters,
    region: InvariantRegion,
    pre_sampling_diversity_expansion_upper_bound: float,
) -> SamplingMultiplierBound:
    """Derive the allele-sampling contribution to an expected diversity bound.

    The simulator samples ``M=max(2, round(2 N_e))`` gene copies. Conditional on
    the pre-sampling allele frequency, expected heterozygosity is multiplied by
    ``1-1/M``. An upper census bound and interaction lower bound provide a safe
    upper bound on ``M``. ``rho`` controls all preceding selection/migration
    effects and must be supplied as a mathematical premise.
    """
    region = validate_invariant_region(region)
    rho = float(pre_sampling_diversity_expansion_upper_bound)
    if rho <= 0.0:
        raise ValueError("pre_sampling_diversity_expansion_upper_bound must be positive")
    effective_upper = max(
        1.0,
        parameters.effective_fraction
        * region.population_upper_bound
        * (1.0 - parameters.skew_penalty * region.next_interaction_lower_bound),
    )
    copies_upper = max(2, floor(2.0 * effective_upper + 0.5))
    sampling_upper = 1.0 - 1.0 / copies_upper
    return SamplingMultiplierBound(
        effective_upper,
        copies_upper,
        sampling_upper,
        rho,
        rho * sampling_upper,
    )


def finite_bin_closure_bound_certificate(
    parameters: DynamicsParameters,
    patch_area: float,
    region: InvariantRegion,
    pre_sampling_diversity_expansion_upper_bound: float,
) -> FiniteBinClosureBoundCertificate:
    """Assemble closure-derived L4 ingredients for a declared refuge patch."""
    recruitment = finite_bin_trait_recruitment_bound(parameters, region)
    cohort = finite_bin_cohort_size_bound(parameters, patch_area, region)
    sampling = sampling_multiplier_bound(
        parameters,
        region,
        pre_sampling_diversity_expansion_upper_bound,
    )
    return FiniteBinClosureBoundCertificate(
        region=validate_invariant_region(region),
        trait_recruitment=recruitment,
        cohort_size=cohort,
        sampling=sampling,
        l4_ready=(
            recruitment.selected_high_trait_probability_lower_bound > 0.0
            and cohort.next_population_lower_bound >= 1
            and sampling.combined_diversity_multiplier_upper_bound < 1.0
        ),
    )
