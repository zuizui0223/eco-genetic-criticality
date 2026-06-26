"""Ensemble summaries for the finite multi-patch criticality simulator.

Theorems L1/L2 concern expected diversity paths. A single simulator run is a
realised stochastic path. This module keeps those levels separate by summarising
replicated runs, including lead probability and mean census-weighted H_alpha.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from statistics import median
from typing import Sequence

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    SimulationResult,
    first_alpha_warning,
    first_high_trait_absence,
    simulate,
)


@dataclass(frozen=True)
class EnsembleSummary:
    """Replicate-level and mean-path summaries with predeclared first passages."""

    replicates: int
    warning_threshold: float
    mean_h_alpha: tuple[float, ...]
    mean_h_gamma: tuple[float, ...]
    trait_absence_times: tuple[int | None, ...]
    alpha_warning_times: tuple[int | None, ...]
    genetic_lead_probability: float
    lead_time_differences: tuple[int, ...]
    median_lead_time_difference: float | None


def simulate_ensemble(
    parameters: DynamicsParameters,
    *,
    replicates: int,
    warning_threshold: float,
) -> tuple[tuple[SimulationResult, ...], EnsembleSummary]:
    """Run seed-indexed replicates and report the correct stochastic lag summaries.

    A lead in replicate r is defined only when both first-passage times exist and
    `tau_H^(r) < tau_trait^(r)`. Runs where either event is absent remain in the
    reported first-passage vectors but are not silently counted as leads.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least one")
    if not 0.0 <= warning_threshold <= 1.0:
        raise ValueError("warning_threshold must lie in [0, 1]")

    results = tuple(
        simulate(replace(parameters, random_seed=parameters.random_seed + replicate))
        for replicate in range(replicates)
    )
    generation_count = len(results[0].snapshots)
    mean_alpha = tuple(
        sum(result.snapshots[generation].h_alpha for result in results) / replicates
        for generation in range(generation_count)
    )
    mean_gamma = tuple(
        sum(result.snapshots[generation].h_gamma for result in results) / replicates
        for generation in range(generation_count)
    )
    trait_times = tuple(first_high_trait_absence(result) for result in results)
    warning_times = tuple(first_alpha_warning(result, warning_threshold) for result in results)
    differences = tuple(
        warning - trait
        for warning, trait in zip(warning_times, trait_times)
        if warning is not None and trait is not None
    )
    lead_count = sum(value < 0 for value in differences)
    return results, EnsembleSummary(
        replicates=replicates,
        warning_threshold=warning_threshold,
        mean_h_alpha=mean_alpha,
        mean_h_gamma=mean_gamma,
        trait_absence_times=trait_times,
        alpha_warning_times=warning_times,
        genetic_lead_probability=lead_count / replicates,
        lead_time_differences=differences,
        median_lead_time_difference=None if not differences else median(differences),
    )


def canonical_interaction_update(
    q: float,
    *,
    area: float,
    area_reference: float,
    feedback_strength: float,
    barrier: float,
) -> float:
    """Return the C1 canonical reduction q_next=sigmoid[kappa((A/Aref)q-theta)]."""
    from causal_model.multipatch_criticality_dynamics import sigmoid

    if not 0.0 <= q <= 1.0:
        raise ValueError("q must lie in [0, 1]")
    if area <= 0.0 or area_reference <= 0.0 or feedback_strength <= 0.0:
        raise ValueError("area, area_reference, and feedback_strength must be positive")
    return sigmoid(feedback_strength * ((area / area_reference) * q - barrier))
