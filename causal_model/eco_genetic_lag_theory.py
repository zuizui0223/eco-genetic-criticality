"""First-passage theory for conditional eco-genetic lag.

This module does not assume that a genetic warning necessarily precedes trait-mode
collapse. It formalises the opposite: a genetic lead is a first-passage property
of a declared diversity recursion and warning threshold.

Let h_t denote expected local (alpha) gene diversity after generation t. For a
specified life cycle, suppose

    h_{t+1} = lambda_t h_t,

where lambda_t>0 is the combined selection-and-transmission diversity multiplier.
The finite-transmission theorem in ``eco_genetic_principles`` supplies the
transmission part; the selection part must be obtained from the model's allele and
trait dynamics.

For a declared trait-collapse time T and warning threshold h_warn, the genetic
lead event is

    tau_H < T,

where tau_H is the first t with h_t <= h_warn. It is equivalent to a cumulative
product inequality. Hence genetic lag is not universal: both lead and no-lead
trajectories exist in the same mathematical family.

The proofs are in ``docs/eco_genetic_criticality_proofs.md``. These utilities
perform exact finite-sequence bookkeeping; they do not infer a life cycle from
patch area.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, log
from typing import Sequence


@dataclass(frozen=True)
class DiversityFirstPassage:
    """A declared diversity trajectory and its relation to trait collapse."""

    initial_diversity: float
    warning_threshold: float
    diversity_multipliers: tuple[float, ...]
    diversity_trajectory: tuple[float, ...]
    trait_collapse_time: int
    warning_time: int | None
    genetic_lead: bool


@dataclass(frozen=True)
class UniformLagBound:
    """A sufficient lead bound under a uniform upper bound on multipliers."""

    initial_diversity: float
    warning_threshold: float
    multiplier_upper_bound: float
    trait_collapse_time: int
    latest_guaranteed_warning_time: int | None
    lead_guaranteed: bool


@dataclass(frozen=True)
class TraitPersistenceBound:
    """A lower-bound certificate for realised high-trait persistence.

    ``certified_trait_loss_lower_bound=T`` means the realised high-trait quantity
    is certified above the extinction threshold for every ``t<T``. Therefore the
    first possible realised trait-loss time is at least ``T``.
    """

    trait_lower_bounds: tuple[float, ...]
    extinction_threshold: float
    certified_trait_loss_lower_bound: int


@dataclass(frozen=True)
class ConditionalLeadCertificate:
    """A sufficient-condition certificate for genetic warning before trait loss."""

    diversity_bound: UniformLagBound
    trait_persistence_bound: TraitPersistenceBound
    lead_guaranteed: bool


def _validate_diversity(value: float, name: str) -> None:
    if not 0.0 < value <= 1.0:
        raise ValueError(f"{name} must lie in (0, 1]")


def diversity_trajectory(initial_diversity: float, multipliers: Sequence[float]) -> tuple[float, ...]:
    """Return h_0,...,h_T for h_{t+1}=lambda_t h_t."""
    _validate_diversity(initial_diversity, "initial_diversity")
    result = [float(initial_diversity)]
    for multiplier in multipliers:
        multiplier = float(multiplier)
        if multiplier <= 0.0:
            raise ValueError("diversity multipliers must be positive")
        next_value = result[-1] * multiplier
        if next_value > 1.0 + 1e-12:
            raise ValueError("multipliers imply diversity above one")
        result.append(min(1.0, next_value))
    return tuple(result)


def first_warning_time(trajectory: Sequence[float], warning_threshold: float) -> int | None:
    """Return the first index at which diversity is at or below the declared threshold."""
    _validate_diversity(warning_threshold, "warning_threshold")
    if not trajectory:
        raise ValueError("trajectory must be nonempty")
    for index, value in enumerate(trajectory):
        _validate_diversity(float(value), "trajectory value")
        if value <= warning_threshold:
            return index
    return None


def assess_genetic_lag(
    initial_diversity: float,
    warning_threshold: float,
    multipliers: Sequence[float],
    trait_collapse_time: int,
) -> DiversityFirstPassage:
    """Assess the exact first-passage criterion for a declared collapse time.

    ``trait_collapse_time=T`` means that the high-trait mode is first absent at
    time T. To assess a lead, the supplied multiplier sequence must define the
    diversity path through time T at minimum.
    """
    trajectory = diversity_trajectory(initial_diversity, multipliers)
    if trait_collapse_time < 0 or trait_collapse_time >= len(trajectory):
        raise ValueError("trait_collapse_time must be an index in the diversity trajectory")
    warning = first_warning_time(trajectory, warning_threshold)
    return DiversityFirstPassage(
        initial_diversity=initial_diversity,
        warning_threshold=warning_threshold,
        diversity_multipliers=tuple(float(value) for value in multipliers),
        diversity_trajectory=trajectory,
        trait_collapse_time=trait_collapse_time,
        warning_time=warning,
        genetic_lead=warning is not None and warning < trait_collapse_time,
    )


def cumulative_multiplier(multipliers: Sequence[float], time: int) -> float:
    """Return product_{s=0}^{time-1} lambda_s, with empty product one."""
    if time < 0 or time > len(multipliers):
        raise ValueError("time must lie between zero and the multiplier length")
    result = 1.0
    for multiplier in multipliers[:time]:
        multiplier = float(multiplier)
        if multiplier <= 0.0:
            raise ValueError("diversity multipliers must be positive")
        result *= multiplier
    return result


def exact_lead_condition(
    initial_diversity: float,
    warning_threshold: float,
    multipliers: Sequence[float],
    trait_collapse_time: int,
) -> bool:
    """Evaluate the exact product condition for tau_H < tau_trait.

    The condition is true iff there exists t<T such that

        product_{s=0}^{t-1} lambda_s <= h_warn/h_0.
    """
    assessment = assess_genetic_lag(
        initial_diversity, warning_threshold, multipliers, trait_collapse_time
    )
    return assessment.genetic_lead


def uniform_upper_multiplier_bound(
    initial_diversity: float,
    warning_threshold: float,
    multiplier_upper_bound: float,
    trait_collapse_time: int,
) -> UniformLagBound:
    """Give a sufficient lead certificate when every lambda_t <= lambda_bar < 1.

    If h_t <= h_0 lambda_bar^t, then the first integer t satisfying

        h_0 lambda_bar^t <= h_warn

    is a latest guaranteed warning time. A strict lead is guaranteed only when
    that time is smaller than the declared trait-collapse time.
    """
    _validate_diversity(initial_diversity, "initial_diversity")
    _validate_diversity(warning_threshold, "warning_threshold")
    if not 0.0 < multiplier_upper_bound < 1.0:
        raise ValueError("multiplier_upper_bound must lie strictly between zero and one")
    if trait_collapse_time < 0:
        raise ValueError("trait_collapse_time must be non-negative")
    if initial_diversity <= warning_threshold:
        time = 0
    else:
        raw = log(warning_threshold / initial_diversity) / log(multiplier_upper_bound)
        time = max(0, ceil(raw))
    return UniformLagBound(
        initial_diversity=initial_diversity,
        warning_threshold=warning_threshold,
        multiplier_upper_bound=multiplier_upper_bound,
        trait_collapse_time=trait_collapse_time,
        latest_guaranteed_warning_time=time,
        lead_guaranteed=time < trait_collapse_time,
    )


def certify_trait_persistence_bound(
    trait_lower_bounds: Sequence[float],
    extinction_threshold: float = 0.0,
) -> TraitPersistenceBound:
    """Convert realised high-trait lower bounds into a trait-loss time bound.

    The inputs are deterministic lower bounds on a non-negative realised
    high-trait quantity, indexed by time. If the bounds are positive above the
    declared extinction threshold through times ``0,...,T-1``, then
    ``tau_trait_realised >= T``. The function returns the largest such prefix
    length. It deliberately does not infer the lower bounds from a simulator.
    """
    if extinction_threshold < 0.0:
        raise ValueError("extinction_threshold must be non-negative")
    observed = tuple(float(value) for value in trait_lower_bounds)
    if not observed:
        raise ValueError("trait_lower_bounds must be nonempty")
    certified_until = 0
    for value in observed:
        if value < 0.0:
            raise ValueError("trait lower bounds must be non-negative")
        if value <= extinction_threshold:
            break
        certified_until += 1
    return TraitPersistenceBound(
        trait_lower_bounds=observed,
        extinction_threshold=float(extinction_threshold),
        certified_trait_loss_lower_bound=certified_until,
    )


def conditional_lead_certificate(
    initial_diversity: float,
    warning_threshold: float,
    multiplier_upper_bound: float,
    trait_lower_bounds: Sequence[float],
    extinction_threshold: float = 0.0,
) -> ConditionalLeadCertificate:
    """Combine genetic decay and trait-persistence bounds into a lead theorem.

    If ``h_t <= h_0 lambda_bar^t`` with ``lambda_bar<1`` and the realised
    high-trait lower bounds certify ``tau_trait_realised >= T``, then
    ``tau_H < tau_trait_realised`` is guaranteed whenever the L2 warning-time
    upper bound is strictly smaller than ``T``.
    """
    persistence = certify_trait_persistence_bound(trait_lower_bounds, extinction_threshold)
    diversity = uniform_upper_multiplier_bound(
        initial_diversity=initial_diversity,
        warning_threshold=warning_threshold,
        multiplier_upper_bound=multiplier_upper_bound,
        trait_collapse_time=persistence.certified_trait_loss_lower_bound,
    )
    return ConditionalLeadCertificate(
        diversity_bound=diversity,
        trait_persistence_bound=persistence,
        lead_guaranteed=diversity.lead_guaranteed,
    )
