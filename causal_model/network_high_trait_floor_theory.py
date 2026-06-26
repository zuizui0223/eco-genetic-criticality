"""Conditional network high-trait persistence theorem for finite-bin recruitment.

The finite-bin simulator updates realised trait abundance before allele migration.
For the declared two-kernel recruitment closure, the low and high kernels have
exactly disjoint support relative to ``high_trait_cutoff``. Thus a common allele
floor and a common realised high-trait mass floor imply a lower bound on the
pre-selection high-trait recruit mass. Viability selection and finite
multinomial recruitment then give a finite-horizon probability lower bound that
every patch retains a declared high-trait mass floor.

The theorem is conditional on lower bounds for the current allele floor,
current resident high-trait mass, interaction, and next cohort size. It does not
prove those region premises and it does not use migration in the same update,
because migration occurs later in the simulator life cycle.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, trait_fitness, trait_grid


@dataclass(frozen=True)
class NetworkHighTraitFloorOneStep:
    """One-step lower envelopes and failure risk for realised high-trait mass."""

    allele_floor: float
    resident_high_trait_mass_floor: float
    interaction_lower_bound: float
    cohort_lower_bound: int
    preselection_high_trait_mass_lower_bound: float
    selected_high_trait_probability_lower_bound: float
    realised_high_trait_mass_floor: float
    per_patch_failure_upper_bound: float
    any_patch_failure_upper_bound: float
    deterministic_mass_floor_possible: bool


@dataclass(frozen=True)
class NetworkHighTraitFloorHorizon:
    """Finite-horizon lower probability for high-trait mass retention in all patches."""

    patches: int
    horizon: int
    one_step: NetworkHighTraitFloorOneStep
    horizon_failure_upper_bound: float
    horizon_retention_probability_lower_bound: float
    certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _require_two_kernel_partition(parameters: DynamicsParameters) -> None:
    if parameters.trait_occupancy_mode != "finite_trait_bin_recruitment":
        raise ValueError("network high-trait theorem requires finite_trait_bin_recruitment")
    if parameters.genotype_trait_recruitment != "two_kernel_recruitment":
        raise ValueError("network high-trait theorem requires two_kernel_recruitment")
    if parameters.low_trait_kernel_center >= parameters.high_trait_cutoff:
        raise ValueError("low kernel must be declared below high_trait_cutoff")
    if parameters.high_trait_kernel_center < parameters.high_trait_cutoff:
        raise ValueError("high kernel must be declared at or above high_trait_cutoff")
    if parameters.high_interaction_benefit < 0.0:
        raise ValueError("theorem requires non-negative high_interaction_benefit")


def preselection_high_trait_mass_lower_bound(
    allele_floor: float,
    resident_high_trait_mass_floor: float,
    parameters: DynamicsParameters,
) -> float:
    """Return the exact lower envelope for high-trait recruit mass.

    In the two-kernel closure the low kernel has zero high-trait mass and the
    high kernel has unit high-trait mass. Therefore allele-linked recruit mass
    in the high region equals p. Mixing it with the resident distribution gives

        (1-w) p + w r_H.
    """
    _require_two_kernel_partition(parameters)
    p = _probability(allele_floor, "allele_floor")
    r = _probability(resident_high_trait_mass_floor, "resident_high_trait_mass_floor")
    return (1.0 - parameters.inheritance_weight) * p + parameters.inheritance_weight * r


def selected_high_trait_probability_lower_bound(
    allele_floor: float,
    resident_high_trait_mass_floor: float,
    interaction_lower_bound: float,
    parameters: DynamicsParameters,
) -> float:
    """Lower-bound selected high-trait recruit probability after viability selection."""
    pre = preselection_high_trait_mass_lower_bound(
        allele_floor,
        resident_high_trait_mass_floor,
        parameters,
    )
    q = _probability(interaction_lower_bound, "interaction_lower_bound")
    grid = trait_grid(parameters)
    high_fitness_lower = min(
        max(parameters.trait_selection_floor, trait_fitness(z, q, parameters))
        for z in grid
        if z >= parameters.high_trait_cutoff
    )
    total_fitness_upper = max(
        max(parameters.trait_selection_floor, trait_fitness(z, 1.0, parameters))
        for z in grid
    )
    return min(1.0, pre * high_fitness_lower / total_fitness_upper)


def _chernoff_fraction_failure(
    trials_lower: int,
    success_probability_lower: float,
    mass_floor: float,
) -> float:
    if trials_lower < 1:
        raise ValueError("cohort_lower_bound must be positive")
    pi = _probability(success_probability_lower, "success_probability_lower")
    r = _probability(mass_floor, "realised_high_trait_mass_floor")
    if r >= pi:
        return 1.0
    mu = trials_lower * pi
    delta = 1.0 - r / pi
    return min(1.0, exp(-mu * delta * delta / 2.0))


def network_high_trait_floor_one_step(
    parameters: DynamicsParameters,
    patches: int,
    allele_floor: float,
    resident_high_trait_mass_floor: float,
    interaction_lower_bound: float,
    cohort_lower_bound: int,
) -> NetworkHighTraitFloorOneStep:
    """Certify one-step high-trait mass retention conditional on a common region.

    The theorem requires that the declared mass floor is large enough to imply
    the simulator's count-based realised-occupancy threshold at the smallest
    cohort size. It then bounds the chance that any patch falls below that mass
    floor after multinomial recruitment.
    """
    if patches < 1:
        raise ValueError("patches must be at least one")
    if cohort_lower_bound < 1:
        raise ValueError("cohort_lower_bound must be at least one")
    mass_floor = _probability(resident_high_trait_mass_floor, "resident_high_trait_mass_floor")
    minimum_occupancy_mass = (
        parameters.realised_high_trait_abundance_threshold / cohort_lower_bound
    )
    if mass_floor < minimum_occupancy_mass:
        raise ValueError(
            "resident_high_trait_mass_floor must imply the realised abundance threshold at cohort_lower_bound"
        )
    pre = preselection_high_trait_mass_lower_bound(
        allele_floor,
        mass_floor,
        parameters,
    )
    selected = selected_high_trait_probability_lower_bound(
        allele_floor,
        mass_floor,
        interaction_lower_bound,
        parameters,
    )
    per_patch_failure = _chernoff_fraction_failure(
        cohort_lower_bound,
        selected,
        mass_floor,
    )
    any_patch_failure = min(1.0, patches * per_patch_failure)
    return NetworkHighTraitFloorOneStep(
        allele_floor=float(allele_floor),
        resident_high_trait_mass_floor=mass_floor,
        interaction_lower_bound=float(interaction_lower_bound),
        cohort_lower_bound=int(cohort_lower_bound),
        preselection_high_trait_mass_lower_bound=pre,
        selected_high_trait_probability_lower_bound=selected,
        realised_high_trait_mass_floor=mass_floor,
        per_patch_failure_upper_bound=per_patch_failure,
        any_patch_failure_upper_bound=any_patch_failure,
        deterministic_mass_floor_possible=selected > mass_floor,
    )


def network_high_trait_floor_horizon(
    parameters: DynamicsParameters,
    patches: int,
    allele_floor: float,
    resident_high_trait_mass_floor: float,
    interaction_lower_bound: float,
    cohort_lower_bound: int,
    horizon: int,
) -> NetworkHighTraitFloorHorizon:
    """Union-bound all-patch high-trait mass failure through a finite horizon."""
    if horizon < 1:
        raise ValueError("horizon must be at least one")
    one_step = network_high_trait_floor_one_step(
        parameters,
        patches,
        allele_floor,
        resident_high_trait_mass_floor,
        interaction_lower_bound,
        cohort_lower_bound,
    )
    failure = min(1.0, horizon * one_step.any_patch_failure_upper_bound)
    retention = max(0.0, 1.0 - failure) if one_step.deterministic_mass_floor_possible else 0.0
    return NetworkHighTraitFloorHorizon(
        patches=patches,
        horizon=horizon,
        one_step=one_step,
        horizon_failure_upper_bound=failure,
        horizon_retention_probability_lower_bound=retention,
        certified=one_step.deterministic_mass_floor_possible and retention > 0.0,
    )
