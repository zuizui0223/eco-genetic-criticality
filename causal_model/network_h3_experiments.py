"""Replicated summaries for the finite H3 extinction–recolonisation lifecycle."""
from __future__ import annotations

from dataclasses import dataclass, replace
from statistics import median
from typing import Sequence

from causal_model.network_h3_lifecycle import (
    NetworkLifecycleParameters,
    NetworkLifecycleResult,
    PatchState,
    first_metapopulation_extinction,
    first_realised_high_trait_loss,
    first_recolonisation,
    simulate_network_lifecycle,
)


@dataclass(frozen=True)
class H3ReplicateSummary:
    """Event times and terminal state from one finite H3 replicate."""

    replicate_index: int
    seed: int
    metapopulation_extinction_time: int | None
    realised_high_trait_loss_time: int | None
    recolonisation_time: int | None
    rescue_occurred: bool
    final_occupied_patch_count: int
    final_high_trait_patch_count: int
    final_h_alpha: float
    final_h_gamma: float
    final_fst: float | None


@dataclass(frozen=True)
class H3EnsembleSummary:
    """Censoring-visible aggregate for a declared H3 landscape experiment."""

    replicates: tuple[H3ReplicateSummary, ...]
    metapopulation_extinction_probability: float
    realised_high_trait_loss_probability: float
    recolonisation_probability: float
    rescue_probability: float
    final_occupied_patch_count_median: float
    final_high_trait_patch_count_median: float
    final_h_alpha_median: float
    final_h_gamma_median: float
    final_fst_median: float | None


def _rescue_occurred(result: NetworkLifecycleResult) -> bool:
    return any(
        snapshot.transitions is not None and any(transition.status == "rescued" for transition in snapshot.transitions)
        for snapshot in result.snapshots[1:]
    )


def _summary(result: NetworkLifecycleResult, replicate_index: int, seed: int) -> H3ReplicateSummary:
    final = result.snapshots[-1]
    return H3ReplicateSummary(
        replicate_index=replicate_index,
        seed=seed,
        metapopulation_extinction_time=first_metapopulation_extinction(result),
        realised_high_trait_loss_time=first_realised_high_trait_loss(result),
        recolonisation_time=first_recolonisation(result),
        rescue_occurred=_rescue_occurred(result),
        final_occupied_patch_count=sum(state.population > 0 for state in final.states),
        final_high_trait_patch_count=sum(state.high_trait_abundance > 0 for state in final.states),
        final_h_alpha=final.h_alpha,
        final_h_gamma=final.h_gamma,
        final_fst=final.fst,
    )


def simulate_h3_ensemble(
    parameters: NetworkLifecycleParameters,
    initial_states: Sequence[PatchState],
    *,
    replicates: int,
) -> tuple[tuple[NetworkLifecycleResult, ...], H3EnsembleSummary]:
    """Run seed-indexed lifecycle replicates and preserve event censoring.

    Event probabilities mean the fraction of finite-horizon replicates in which
    the event occurred.  The per-replicate event times remain available so
    censored runs are not converted into terminal-generation times.
    """
    if replicates < 1:
        raise ValueError("replicates must be at least one")
    results = tuple(
        simulate_network_lifecycle(
            replace(parameters, random_seed=parameters.random_seed + replicate_index),
            initial_states,
        )
        for replicate_index in range(replicates)
    )
    summaries = tuple(
        _summary(result, replicate_index, parameters.random_seed + replicate_index)
        for replicate_index, result in enumerate(results)
    )
    observed_fst = tuple(summary.final_fst for summary in summaries if summary.final_fst is not None)
    return results, H3EnsembleSummary(
        replicates=summaries,
        metapopulation_extinction_probability=sum(
            summary.metapopulation_extinction_time is not None for summary in summaries
        ) / replicates,
        realised_high_trait_loss_probability=sum(
            summary.realised_high_trait_loss_time is not None for summary in summaries
        ) / replicates,
        recolonisation_probability=sum(summary.recolonisation_time is not None for summary in summaries) / replicates,
        rescue_probability=sum(summary.rescue_occurred for summary in summaries) / replicates,
        final_occupied_patch_count_median=median(summary.final_occupied_patch_count for summary in summaries),
        final_high_trait_patch_count_median=median(summary.final_high_trait_patch_count for summary in summaries),
        final_h_alpha_median=median(summary.final_h_alpha for summary in summaries),
        final_h_gamma_median=median(summary.final_h_gamma for summary in summaries),
        final_fst_median=None if not observed_fst else median(observed_fst),
    )
