"""Ensemble summaries for the finite multi-patch criticality simulator.

Theorems L1/L2 concern expected diversity paths. A single simulator run is a
realised stochastic path. This module keeps those levels separate and reports
both the conditional lead probability among valid event pairs and the
unconditional observed-lead fraction across all replicates.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Sequence

from causal_model.first_passage_reporting import (
    FirstPassageComparison,
    compare_first_passage_times,
)
from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    SimulationResult,
    first_alpha_warning,
    first_high_trait_absence,
    simulate,
)


@dataclass(frozen=True)
class EnsembleSummary:
    """Replicate-level and mean-path summaries with explicit censoring denominators."""

    replicates: int
    warning_threshold: float
    mean_h_alpha: tuple[float, ...]
    mean_h_gamma: tuple[float, ...]
    trait_absence_times: tuple[int | None, ...]
    alpha_warning_times: tuple[int | None, ...]
    alpha_warning_vs_trait_absence: FirstPassageComparison

    @property
    def valid_event_pair_count(self) -> int:
        """Number of replicates where both alpha warning and trait loss occurred."""
        return self.alpha_warning_vs_trait_absence.valid_pair_count

    @property
    def event_pair_observability(self) -> float:
        """Fraction of all replicates that support the ordered-event comparison."""
        return self.alpha_warning_vs_trait_absence.valid_pair_probability

    @property
    def conditional_genetic_lead_probability(self) -> float | None:
        """P(warning precedes trait loss | both events observed)."""
        return self.alpha_warning_vs_trait_absence.conditional_lead_probability

    @property
    def unconditional_observed_lead_fraction(self) -> float:
        """P(observed warning lead) across every replicate, including censored runs."""
        return self.alpha_warning_vs_trait_absence.unconditional_observed_lead_fraction

    @property
    def genetic_lead_probability(self) -> float:
        """Backward-compatible alias for the unconditional observed-lead fraction.

        Use ``conditional_genetic_lead_probability`` for the estimate conditional
        on valid event pairs.  The older unqualified name had an ambiguous
        denominator and is retained only to avoid breaking callers.
        """
        return self.unconditional_observed_lead_fraction

    @property
    def lead_time_differences(self) -> tuple[int, ...]:
        """Warning-minus-trait times among valid event pairs only."""
        return self.alpha_warning_vs_trait_absence.time_differences

    @property
    def median_lead_time_difference(self) -> float | None:
        return self.alpha_warning_vs_trait_absence.median_time_difference


def simulate_ensemble(
    parameters: DynamicsParameters,
    *,
    replicates: int,
    warning_threshold: float,
) -> tuple[tuple[SimulationResult, ...], EnsembleSummary]:
    """Run seed-indexed replicates with a censoring-aware stochastic lag summary.

    A lead in replicate ``r`` is defined only when both first-passage times exist
    and ``tau_H^(r) < tau_trait^(r)``.  The returned comparison exposes the
    valid-pair denominator, conditional lead probability, and all-replicate
    observed-lead fraction separately.
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
    comparison = compare_first_passage_times(warning_times, trait_times)
    return results, EnsembleSummary(
        replicates=replicates,
        warning_threshold=warning_threshold,
        mean_h_alpha=mean_alpha,
        mean_h_gamma=mean_gamma,
        trait_absence_times=trait_times,
        alpha_warning_times=warning_times,
        alpha_warning_vs_trait_absence=comparison,
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
