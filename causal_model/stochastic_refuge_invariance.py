"""High-probability invariant-region certificate for a finite-bin refuge patch.

This is an exact specialisation of the declared simulator to one refuge patch
with ``migration_rate == 0``. The certificate is deliberately conditional on a
rectangular state region R. It derives deterministic lower and upper envelopes
for q and N, and Chernoff lower-tail envelopes for the two finite samples:
allele frequency and realised high-trait recruitment.

A successful certificate proves a lower bound on the probability that the patch
remains in R through a finite horizon. Such a refuge is sufficient to prevent
the simulator's all-patch realised high-trait-loss event over that horizon.
It does not prove invariance for arbitrary migration, arbitrary landscape, or
unbounded time.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, floor

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, sigmoid, trait_fitness, trait_grid


@dataclass(frozen=True)
class RefugeRegion:
    """One-patch state rectangle used by the certificate."""

    interaction_lower: float
    allele_lower: float
    high_trait_mass_lower: float
    population_lower: int
    population_upper: int


@dataclass(frozen=True)
class RefugeOneStepBound:
    """Deterministic envelopes and stochastic failure bounds for one update."""

    interaction_next_lower: float
    selected_allele_lower: float
    population_next_lower: int
    population_next_upper: int
    high_trait_recruit_lower: float
    gene_copy_lower: int
    allele_failure_upper_bound: float
    trait_failure_upper_bound: float
    deterministic_region_closed: bool
    one_step_retention_probability_lower_bound: float


@dataclass(frozen=True)
class RefugeHorizonCertificate:
    """Finite-horizon lower bound for remaining inside the refuge region."""

    horizon: int
    one_step: RefugeOneStepBound
    horizon_failure_upper_bound: float
    horizon_retention_probability_lower_bound: float
    certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def validate_refuge_region(region: RefugeRegion) -> RefugeRegion:
    """Validate the rectangle but do not claim it is invariant."""
    q = _probability(region.interaction_lower, "interaction_lower")
    p = _probability(region.allele_lower, "allele_lower")
    r = _probability(region.high_trait_mass_lower, "high_trait_mass_lower")
    if region.population_lower < 1:
        raise ValueError("population_lower must be at least one")
    if region.population_upper < region.population_lower:
        raise ValueError("population_upper must be at least population_lower")
    return RefugeRegion(q, p, r, int(region.population_lower), int(region.population_upper))


def _feedback_weights(parameters: DynamicsParameters) -> tuple[float, float, float]:
    if parameters.q_feedback_alpha is None and parameters.q_feedback_gamma_allele is None:
        return (
            parameters.interaction_memory_weight,
            parameters.q_feedback_beta_trait,
            1.0 - parameters.interaction_memory_weight,
        )
    alpha = parameters.interaction_memory_weight if parameters.q_feedback_alpha is None else parameters.q_feedback_alpha
    gamma = 0.0 if parameters.q_feedback_gamma_allele is None else parameters.q_feedback_gamma_allele
    return alpha, parameters.q_feedback_beta_trait, gamma


def _round_lower(value: float) -> int:
    """Safe lower envelope for Python's round(x), including half-even ties."""
    return ceil(value - 0.5)


def _round_upper(value: float) -> int:
    """Safe upper envelope for Python's round(x), including half-even ties."""
    return floor(value + 0.5)


def _selected_high_allele_lower(parameters: DynamicsParameters, q_next_lower: float, p_lower: float) -> float:
    margin = trait_fitness(1.0, q_next_lower, parameters) - parameters.viability_threshold
    fitness = max(1e-12, 1.0 + parameters.selection_strength * margin)
    return p_lower * fitness / (p_lower * fitness + (1.0 - p_lower))


def _selected_trait_high_probability_lower(
    parameters: DynamicsParameters,
    q_lower: float,
    p_lower: float,
    r_lower: float,
) -> float:
    if parameters.genotype_trait_recruitment != "two_kernel_recruitment":
        raise ValueError("certificate requires two_kernel_recruitment")
    if parameters.low_trait_kernel_center >= parameters.high_trait_cutoff:
        raise ValueError("low trait kernel must be centred below high_trait_cutoff")
    if parameters.high_trait_kernel_center < parameters.high_trait_cutoff:
        raise ValueError("high trait kernel must be centred in high-trait region")
    pre_high = (1.0 - parameters.inheritance_weight) * p_lower + parameters.inheritance_weight * r_lower
    grid = trait_grid(parameters)
    high_fitness = min(
        max(parameters.trait_selection_floor, trait_fitness(z, q_lower, parameters))
        for z in grid
        if z >= parameters.high_trait_cutoff
    )
    total_upper = max(
        max(parameters.trait_selection_floor, trait_fitness(z, 1.0, parameters))
        for z in grid
    )
    return min(1.0, pre_high * high_fitness / total_upper)


def _chernoff_lower_tail_failure(trials_lower: int, probability_lower: float, realised_fraction_lower: float) -> float:
    """Upper bound P[X/trials < realised_fraction_lower] for X>=Binomial(n,p)."""
    if trials_lower < 1:
        raise ValueError("trials_lower must be positive")
    probability_lower = _probability(probability_lower, "probability_lower")
    realised_fraction_lower = _probability(realised_fraction_lower, "realised_fraction_lower")
    if realised_fraction_lower >= probability_lower:
        return 1.0
    mu = trials_lower * probability_lower
    delta = 1.0 - realised_fraction_lower / probability_lower
    return min(1.0, exp(-mu * delta * delta / 2.0))


def one_step_refuge_bound(
    parameters: DynamicsParameters,
    patch_area: float,
    region: RefugeRegion,
) -> RefugeOneStepBound:
    """Derive a one-step high-probability refuge certificate.

    The certificate exactly matches the simulator only for a one-patch,
    no-migration submodel. For a multipatch system it can be used as a declared
    refuge-patch abstraction only after separately justifying that migration does
    not lower the patch's incoming allele frequency below the stated region.
    """
    region = validate_refuge_region(region)
    if patch_area <= 0.0:
        raise ValueError("patch_area must be positive")
    if parameters.trait_occupancy_mode != "finite_trait_bin_recruitment":
        raise ValueError("certificate requires finite_trait_bin_recruitment")
    if parameters.migration_rate != 0.0:
        raise ValueError("exact one-patch certificate requires migration_rate == 0")

    carrying = parameters.density_capacity * patch_area
    density_lower = min(1.0, region.population_lower / carrying)
    alpha, beta, gamma = _feedback_weights(parameters)
    support_lower = alpha * region.interaction_lower + beta * region.high_trait_mass_lower + gamma * region.allele_lower
    q_next_lower = sigmoid(
        parameters.interaction_feedback
        * ((patch_area / parameters.area_reference) * density_lower * support_lower - parameters.interaction_barrier)
    )
    p_selected_lower = _selected_high_allele_lower(parameters, q_next_lower, region.allele_lower)

    exponent_lower = (
        parameters.baseline_growth
        + parameters.interaction_growth * q_next_lower
        + parameters.high_allele_growth * p_selected_lower
        - region.population_upper / carrying
    )
    exponent_upper = (
        parameters.baseline_growth
        + parameters.interaction_growth
        + parameters.high_allele_growth
        - region.population_lower / carrying
    )
    n_next_lower = max(1, _round_lower(region.population_lower * exp(exponent_lower)))
    n_next_upper = max(1, _round_upper(region.population_upper * exp(exponent_upper)))

    trait_pi_lower = _selected_trait_high_probability_lower(
        parameters,
        region.interaction_lower,
        region.allele_lower,
        region.high_trait_mass_lower,
    )
    trait_failure = _chernoff_lower_tail_failure(
        n_next_lower,
        trait_pi_lower,
        region.high_trait_mass_lower,
    )

    effective_lower = max(
        1.0,
        parameters.effective_fraction * n_next_lower * (1.0 - parameters.skew_penalty),
    )
    gene_copy_lower = max(2, _round_lower(2.0 * effective_lower))
    allele_failure = _chernoff_lower_tail_failure(
        gene_copy_lower,
        p_selected_lower,
        region.allele_lower,
    )

    deterministic_closed = (
        q_next_lower >= region.interaction_lower
        and n_next_lower >= region.population_lower
        and n_next_upper <= region.population_upper
    )
    retention = max(0.0, 1.0 - allele_failure - trait_failure) if deterministic_closed else 0.0
    return RefugeOneStepBound(
        interaction_next_lower=q_next_lower,
        selected_allele_lower=p_selected_lower,
        population_next_lower=n_next_lower,
        population_next_upper=n_next_upper,
        high_trait_recruit_lower=trait_pi_lower,
        gene_copy_lower=gene_copy_lower,
        allele_failure_upper_bound=allele_failure,
        trait_failure_upper_bound=trait_failure,
        deterministic_region_closed=deterministic_closed,
        one_step_retention_probability_lower_bound=retention,
    )


def finite_horizon_refuge_certificate(
    parameters: DynamicsParameters,
    patch_area: float,
    region: RefugeRegion,
    horizon: int,
) -> RefugeHorizonCertificate:
    """Union-bound the one-step refuge failure probability over a finite horizon."""
    if horizon < 1:
        raise ValueError("horizon must be at least one")
    one_step = one_step_refuge_bound(parameters, patch_area, region)
    failure = min(1.0, horizon * (one_step.allele_failure_upper_bound + one_step.trait_failure_upper_bound))
    retention = max(0.0, 1.0 - failure) if one_step.deterministic_region_closed else 0.0
    return RefugeHorizonCertificate(
        horizon=horizon,
        one_step=one_step,
        horizon_failure_upper_bound=failure,
        horizon_retention_probability_lower_bound=retention,
        certified=one_step.deterministic_region_closed and retention > 0.0,
    )
