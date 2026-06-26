"""Finite-horizon moving allele-corridor theorem.

A fixed interval below one cannot be invariant under uniformly positive
high-allele selection: the deterministic selection map pushes its upper endpoint
upward.  This module supplies the correct replacement for the finite-bin
simulator: a time-indexed corridor whose endpoints may move with selection.

At generation t, assume every patch frequency lies in [L_t, U_t], interaction
lies in [q_min, q_max], and next census size is at least N_min.  Local selection
and census-weighted migration place every pre-sampling frequency in

    [s(L_t, f_min), s(U_t, f_max)].

A declared next corridor [L_(t+1), U_(t+1)] that strictly contains this envelope
has a two-sided binomial Chernoff/KL sampling-escape bound.  Union bounds across
patches and generations yield a finite-horizon retention probability.

The theorem is conditional on q/N bounds.  It does not derive those ecological
bounds from the complete coupled system.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, log
from typing import Sequence

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, trait_fitness


@dataclass(frozen=True)
class CorridorStep:
    """One transition between two declared frequency corridors."""

    generation: int
    lower_before: float
    upper_before: float
    lower_after: float
    upper_after: float
    selected_lower_envelope: float
    selected_upper_envelope: float
    gene_copy_lower_bound: int
    lower_escape_probability_upper_bound: float
    upper_escape_probability_upper_bound: float
    any_patch_escape_probability_upper_bound: float
    deterministic_corridor_contains_envelope: bool


@dataclass(frozen=True)
class MovingCorridorCertificate:
    """A finite-horizon high-probability common-corridor certificate."""

    patches: int
    steps: tuple[CorridorStep, ...]
    horizon_escape_probability_upper_bound: float
    horizon_retention_probability_lower_bound: float
    certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _closed_interval(lower: float, upper: float, name: str) -> tuple[float, float]:
    lower = _probability(lower, f"{name}_lower")
    upper = _probability(upper, f"{name}_upper")
    if lower > upper:
        raise ValueError(f"{name}_lower must not exceed {name}_upper")
    return lower, upper


def _binary_kl(x: float, p: float) -> float:
    """Bernoulli relative entropy D(x || p)."""
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


def selected_frequency(frequency: float, high_allele_fitness: float) -> float:
    """The simulator's deterministic one-allele selection map."""
    frequency = _probability(frequency, "frequency")
    if high_allele_fitness <= 0.0:
        raise ValueError("high_allele_fitness must be positive")
    return frequency * high_allele_fitness / (
        frequency * high_allele_fitness + (1.0 - frequency)
    )


def _fitness(interaction: float, parameters: DynamicsParameters) -> float:
    margin = trait_fitness(1.0, interaction, parameters) - parameters.viability_threshold
    return max(1e-12, 1.0 + parameters.selection_strength * margin)


def _gene_copy_lower_bound(population_lower_bound: int, parameters: DynamicsParameters) -> int:
    if population_lower_bound < 1:
        raise ValueError("population_lower_bound must be at least one")
    effective_lower = max(
        1.0,
        parameters.effective_fraction
        * population_lower_bound
        * (1.0 - parameters.skew_penalty),
    )
    return max(2, ceil(2.0 * effective_lower - 0.5))


def _lower_tail(copies: int, probability_lower: float, threshold: float) -> float:
    if threshold >= probability_lower:
        return 1.0
    return min(1.0, exp(-copies * _binary_kl(threshold, probability_lower)))


def _upper_tail(copies: int, probability_upper: float, threshold: float) -> float:
    if threshold <= probability_upper:
        return 1.0
    return min(1.0, exp(-copies * _binary_kl(threshold, probability_upper)))


def moving_corridor_step(
    parameters: DynamicsParameters,
    patches: int,
    lower_before: float,
    upper_before: float,
    lower_after: float,
    upper_after: float,
    interaction_lower_bound: float,
    interaction_upper_bound: float,
    population_lower_bound: int,
    generation: int = 0,
) -> CorridorStep:
    """Certify one declared moving-corridor transition.

    The next corridor must lie strictly outside the deterministic
    selection/migration envelope to obtain a nontrivial sampling bound.
    """
    if patches < 1:
        raise ValueError("patches must be at least one")
    _require_monotone_high_selection(parameters)
    lower_before, upper_before = _closed_interval(lower_before, upper_before, "before")
    lower_after, upper_after = _closed_interval(lower_after, upper_after, "after")
    q_lower, q_upper = _closed_interval(
        interaction_lower_bound,
        interaction_upper_bound,
        "interaction",
    )
    selected_lower = selected_frequency(lower_before, _fitness(q_lower, parameters))
    selected_upper = selected_frequency(upper_before, _fitness(q_upper, parameters))
    copies = _gene_copy_lower_bound(population_lower_bound, parameters)
    lower_escape = _lower_tail(copies, selected_lower, lower_after)
    upper_escape = _upper_tail(copies, selected_upper, upper_after)
    per_patch = min(1.0, lower_escape + upper_escape)
    return CorridorStep(
        generation=int(generation),
        lower_before=lower_before,
        upper_before=upper_before,
        lower_after=lower_after,
        upper_after=upper_after,
        selected_lower_envelope=selected_lower,
        selected_upper_envelope=selected_upper,
        gene_copy_lower_bound=copies,
        lower_escape_probability_upper_bound=lower_escape,
        upper_escape_probability_upper_bound=upper_escape,
        any_patch_escape_probability_upper_bound=min(1.0, patches * per_patch),
        deterministic_corridor_contains_envelope=(
            lower_after < selected_lower and selected_upper < upper_after
        ),
    )


def moving_corridor_certificate(
    parameters: DynamicsParameters,
    patches: int,
    lower_path: Sequence[float],
    upper_path: Sequence[float],
    interaction_lower_bound: float,
    interaction_upper_bound: float,
    population_lower_bound: int,
) -> MovingCorridorCertificate:
    """Union-bound retention of a declared corridor path.

    ``lower_path`` and ``upper_path`` each contain the initial endpoint followed
    by one endpoint per transition.  No independence between patches or times is
    assumed.
    """
    if len(lower_path) != len(upper_path):
        raise ValueError("lower_path and upper_path must have equal length")
    if len(lower_path) < 2:
        raise ValueError("a moving corridor needs at least one transition")
    steps = tuple(
        moving_corridor_step(
            parameters,
            patches,
            lower_path[t],
            upper_path[t],
            lower_path[t + 1],
            upper_path[t + 1],
            interaction_lower_bound,
            interaction_upper_bound,
            population_lower_bound,
            generation=t + 1,
        )
        for t in range(len(lower_path) - 1)
    )
    escape = min(1.0, sum(step.any_patch_escape_probability_upper_bound for step in steps))
    all_contain = all(step.deterministic_corridor_contains_envelope for step in steps)
    retention = max(0.0, 1.0 - escape) if all_contain else 0.0
    return MovingCorridorCertificate(
        patches=patches,
        steps=steps,
        horizon_escape_probability_upper_bound=escape,
        horizon_retention_probability_lower_bound=retention,
        certified=all_contain and retention > 0.0,
    )


def selection_shifted_corridor(
    parameters: DynamicsParameters,
    lower_initial: float,
    upper_initial: float,
    interaction_lower_bound: float,
    interaction_upper_bound: float,
    generations: int,
    slack: float,
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Build a simple upward-moving corridor with a declared sampling slack.

    The corridor is a construction aid, not a proof by itself.  Positive ``slack``
    puts the next lower endpoint below its deterministic envelope and the next
    upper endpoint above it, leaving room for finite-sampling deviations.
    """
    if generations < 1:
        raise ValueError("generations must be at least one")
    slack = _probability(slack, "slack")
    _require_monotone_high_selection(parameters)
    lower, upper = _closed_interval(lower_initial, upper_initial, "initial")
    q_lower, q_upper = _closed_interval(interaction_lower_bound, interaction_upper_bound, "interaction")
    lower_values = [lower]
    upper_values = [upper]
    for _ in range(generations):
        selected_lower = selected_frequency(lower_values[-1], _fitness(q_lower, parameters))
        selected_upper = selected_frequency(upper_values[-1], _fitness(q_upper, parameters))
        lower_values.append(max(0.0, selected_lower - slack))
        upper_values.append(min(1.0, selected_upper + slack))
    return tuple(lower_values), tuple(upper_values)
