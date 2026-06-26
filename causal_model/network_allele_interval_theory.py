"""High-probability allele-interval retention for the finite-bin network closure.

This module closes the interval premise required by a restricted H-alpha
multiplier theorem.  It follows the simulator life cycle exactly:

    allele selection -> census-weighted migration -> finite gene-copy sampling.

Selection maps an input frequency interval to a selected interval.  Migration
is a convex combination and cannot leave that common selected interval.  The
only interval-breaking step is finite sampling, which is bounded using binary
relative-entropy Chernoff inequalities and a union bound across patches and
finite time.

The result is conditional on declared bounds for interaction and next census
size.  It does not claim an invariant region for the complete q/N/trait system.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, log

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, trait_fitness


@dataclass(frozen=True)
class AlleleIntervalOneStep:
    """Deterministic selection/migration envelope and sampling escape risk."""

    allele_lower_bound: float
    allele_upper_bound: float
    interaction_lower_bound: float
    interaction_upper_bound: float
    population_lower_bound: int
    selected_lower_bound: float
    selected_upper_bound: float
    migrated_lower_bound: float
    migrated_upper_bound: float
    gene_copy_lower_bound: int
    lower_escape_probability_upper_bound: float
    upper_escape_probability_upper_bound: float
    per_patch_escape_probability_upper_bound: float
    any_patch_escape_probability_upper_bound: float
    deterministic_inner_interval: bool


@dataclass(frozen=True)
class AlleleIntervalHorizon:
    """Finite-horizon common-interval retention certificate."""

    patches: int
    horizon: int
    one_step: AlleleIntervalOneStep
    horizon_escape_probability_upper_bound: float
    horizon_retention_probability_lower_bound: float
    certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _validate_interval(lower: float, upper: float) -> tuple[float, float]:
    lower = _probability(lower, "allele_lower_bound")
    upper = _probability(upper, "allele_upper_bound")
    if not lower < upper:
        raise ValueError("allele_lower_bound must be strictly smaller than allele_upper_bound")
    return lower, upper


def _binary_kl(x: float, p: float) -> float:
    """Bernoulli relative entropy D(x || p), with stable boundary limits."""
    x = _probability(x, "x")
    p = _probability(p, "p")
    if x == p:
        return 0.0
    first = 0.0 if x == 0.0 else x * log(x / p)
    second = 0.0 if x == 1.0 else (1.0 - x) * log((1.0 - x) / (1.0 - p))
    return first + second


def _require_monotone_high_selection(parameters: DynamicsParameters) -> None:
    if parameters.high_interaction_benefit < 0.0:
        raise ValueError("theorem requires non-negative high_interaction_benefit")
    if parameters.selection_strength < 0.0:
        raise ValueError("theorem requires non-negative selection_strength")
    if parameters.skew_penalty < 0.0:
        raise ValueError("theorem requires non-negative skew_penalty")


def _high_allele_fitness(interaction: float, parameters: DynamicsParameters) -> float:
    margin = trait_fitness(1.0, interaction, parameters) - parameters.viability_threshold
    return max(1e-12, 1.0 + parameters.selection_strength * margin)


def selected_frequency(frequency: float, high_allele_fitness: float) -> float:
    """The simulator's deterministic one-allele selection map."""
    frequency = _probability(frequency, "frequency")
    if high_allele_fitness <= 0.0:
        raise ValueError("high_allele_fitness must be positive")
    return frequency * high_allele_fitness / (
        frequency * high_allele_fitness + (1.0 - frequency)
    )


def _gene_copy_lower_bound(population_lower_bound: int, parameters: DynamicsParameters) -> int:
    if population_lower_bound < 1:
        raise ValueError("population_lower_bound must be at least one")
    # n_eff = N*effective_fraction*(1-skew*q).  With non-negative skew,
    # q=1 is a safe lower envelope over the simulator's [0,1] state domain.
    effective_lower = max(
        1.0,
        parameters.effective_fraction
        * population_lower_bound
        * (1.0 - parameters.skew_penalty),
    )
    # ceil(x - 1/2) is a lower envelope for Python's half-even round(x).
    return max(2, ceil(2.0 * effective_lower - 0.5))


def _lower_tail_bound(copies: int, success_lower: float, threshold: float) -> float:
    """Bound P(X/copies < threshold) when X stochastically dominates Bin(copies,p)."""
    if threshold >= success_lower:
        return 1.0
    return min(1.0, exp(-copies * _binary_kl(threshold, success_lower)))


def _upper_tail_bound(copies: int, success_upper: float, threshold: float) -> float:
    """Bound P(X/copies > threshold) when X is stochastically below Bin(copies,p)."""
    if threshold <= success_upper:
        return 1.0
    return min(1.0, exp(-copies * _binary_kl(threshold, success_upper)))


def network_allele_interval_one_step(
    parameters: DynamicsParameters,
    patches: int,
    allele_lower_bound: float,
    allele_upper_bound: float,
    interaction_lower_bound: float,
    interaction_upper_bound: float,
    population_lower_bound: int,
) -> AlleleIntervalOneStep:
    """Bound one-step escape from a common allele interval.

    Assume every patch begins in ``[allele_lower_bound, allele_upper_bound]``,
    interaction lies in the declared interval, and each next census is at least
    ``population_lower_bound``.  The returned bound is exact for the deterministic
    selection/migration envelope and conservative for sampling.
    """
    if patches < 1:
        raise ValueError("patches must be at least one")
    _require_monotone_high_selection(parameters)
    lower, upper = _validate_interval(allele_lower_bound, allele_upper_bound)
    q_lower, q_upper = _validate_interval(interaction_lower_bound, interaction_upper_bound)
    fitness_lower = _high_allele_fitness(q_lower, parameters)
    fitness_upper = _high_allele_fitness(q_upper, parameters)
    selected_lower = selected_frequency(lower, fitness_lower)
    selected_upper = selected_frequency(upper, fitness_upper)
    copies = _gene_copy_lower_bound(population_lower_bound, parameters)
    lower_escape = _lower_tail_bound(copies, selected_lower, lower)
    upper_escape = _upper_tail_bound(copies, selected_upper, upper)
    per_patch = min(1.0, lower_escape + upper_escape)
    any_patch = min(1.0, patches * per_patch)
    return AlleleIntervalOneStep(
        allele_lower_bound=lower,
        allele_upper_bound=upper,
        interaction_lower_bound=q_lower,
        interaction_upper_bound=q_upper,
        population_lower_bound=int(population_lower_bound),
        selected_lower_bound=selected_lower,
        selected_upper_bound=selected_upper,
        migrated_lower_bound=selected_lower,
        migrated_upper_bound=selected_upper,
        gene_copy_lower_bound=copies,
        lower_escape_probability_upper_bound=lower_escape,
        upper_escape_probability_upper_bound=upper_escape,
        per_patch_escape_probability_upper_bound=per_patch,
        any_patch_escape_probability_upper_bound=any_patch,
        deterministic_inner_interval=(selected_lower > lower and selected_upper < upper),
    )


def network_allele_interval_horizon(
    parameters: DynamicsParameters,
    patches: int,
    allele_lower_bound: float,
    allele_upper_bound: float,
    interaction_lower_bound: float,
    interaction_upper_bound: float,
    population_lower_bound: int,
    horizon: int,
) -> AlleleIntervalHorizon:
    """Union-bound common-interval escape through a finite horizon."""
    if horizon < 1:
        raise ValueError("horizon must be at least one")
    one_step = network_allele_interval_one_step(
        parameters,
        patches,
        allele_lower_bound,
        allele_upper_bound,
        interaction_lower_bound,
        interaction_upper_bound,
        population_lower_bound,
    )
    escape = min(1.0, horizon * one_step.any_patch_escape_probability_upper_bound)
    retention = max(0.0, 1.0 - escape) if one_step.deterministic_inner_interval else 0.0
    return AlleleIntervalHorizon(
        patches=patches,
        horizon=horizon,
        one_step=one_step,
        horizon_escape_probability_upper_bound=escape,
        horizon_retention_probability_lower_bound=retention,
        certified=one_step.deterministic_inner_interval and retention > 0.0,
    )
