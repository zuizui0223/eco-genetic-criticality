"""Finite simulation closure with symmetric high/low allele-state mutation.

This module leaves the legacy simulator unchanged.  It reproduces the declared
life cycle exactly, except that after selection and migration and before finite
genetic drift each patch's high-allele frequency is transformed by

    p_mut = mu + (1 - 2 mu) p.

``mu`` is a per-generation *symmetric allele-state transition* probability:
high-to-low and low-to-high transitions occur at the same rate.  It is an
explicit new numerical closure, not a biological estimate and not a change to
the theorem layer.

At ``mu = 0`` the function delegates to ``simulate`` directly.  Thus legacy
runs retain identical stochastic trajectories and results.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from random import Random
from typing import Iterator

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    SimulationResult,
    _binomial,
    _effective_size,
    _initial_values,
    _normalise_distribution,
    _snapshot,
    interaction_support_signal,
    sigmoid,
    simulate,
    trait_fitness,
    update_trait_abundance,
    update_trait_distribution,
)


def validate_symmetric_allele_mutation_rate(mutation_rate: float) -> float:
    """Validate and return a symmetric allele-state mutation probability."""
    rate = float(mutation_rate)
    if not 0.0 <= rate < 0.5:
        raise ValueError("symmetric_allele_mutation_rate must lie in [0, 0.5)")
    return rate


def apply_symmetric_allele_mutation(frequency: float, mutation_rate: float) -> float:
    """Apply p -> mu + (1 - 2 mu) p, preserving p when mu is zero."""
    rate = validate_symmetric_allele_mutation_rate(mutation_rate)
    if not 0.0 <= frequency <= 1.0:
        raise ValueError("frequency must lie in [0, 1]")
    return rate + (1.0 - 2.0 * rate) * frequency


def simulate_with_symmetric_allele_mutation(
    parameters: DynamicsParameters,
    *,
    mutation_rate: float,
) -> SimulationResult:
    """Run the finite life cycle with post-migration, pre-drift mutation.

    Mutation is placed immediately after allele-frequency migration and before
    the existing binomial drift draw.  Trait recruitment at generation t uses
    the resident p_t exactly as in the base closure; mutation affects p_{t+1}
    and therefore subsequent recruitment and q-feedback.
    """
    rate = validate_symmetric_allele_mutation_rate(mutation_rate)
    if rate == 0.0:
        return simulate(parameters)

    rng = Random(parameters.random_seed)
    population, interaction, frequency, trait_distribution, trait_abundance = _initial_values(parameters)
    snapshots = [_snapshot(0, population, interaction, frequency, trait_distribution, trait_abundance, parameters)]

    for generation in range(1, parameters.generations + 1):
        current_occupancy = snapshots[-1].trait_occupancy
        current_high_mass = tuple(summary.high_trait_mass for summary in current_occupancy)
        carrying = tuple(parameters.density_capacity * area for area in parameters.patch_areas)
        density = tuple(min(1.0, n / k) for n, k in zip(population, carrying))
        support = tuple(
            interaction_support_signal(q, x_h, p, parameters)
            for q, x_h, p in zip(interaction, current_high_mass, frequency)
        )
        q_next = tuple(
            sigmoid(
                parameters.interaction_feedback
                * ((area / parameters.area_reference) * dens * signal - parameters.interaction_barrier)
            )
            for area, dens, signal in zip(parameters.patch_areas, density, support)
        )

        selected: list[float] = []
        for q, p in zip(q_next, frequency):
            high_margin = trait_fitness(1.0, q, parameters) - parameters.viability_threshold
            high_fitness = max(1e-12, 1.0 + parameters.selection_strength * high_margin)
            mean_fitness = p * high_fitness + (1.0 - p)
            selected.append(p * high_fitness / mean_fitness)

        weights = tuple(float(n) for n in population)
        selected_mean = sum(weight * p for weight, p in zip(weights, selected)) / sum(weights)
        migrated = tuple(
            (1.0 - parameters.migration_rate) * p + parameters.migration_rate * selected_mean
            for p in selected
        )
        mutated = tuple(apply_symmetric_allele_mutation(p, rate) for p in migrated)

        next_population: list[int] = []
        for n, k, q, p in zip(population, carrying, q_next, selected):
            exponent = parameters.baseline_growth + parameters.interaction_growth * q + parameters.high_allele_growth * p - n / k
            next_population.append(max(1, round(n * exp(exponent))))

        if parameters.trait_occupancy_mode == "finite_trait_bin_recruitment":
            next_trait_abundance = tuple(
                update_trait_abundance(abundance, q, p, n_next, parameters, rng)
                for abundance, q, p, n_next in zip(trait_abundance, interaction, frequency, next_population)
            )
            next_trait_distribution = tuple(_normalise_distribution(row) for row in next_trait_abundance)
        else:
            next_trait_distribution = tuple(
                update_trait_distribution(mu, q, parameters, p)
                for mu, q, p in zip(trait_distribution, interaction, frequency)
            )
            next_trait_abundance = tuple(
                _abundance_from_distribution(distribution, n_next)
                for distribution, n_next in zip(next_trait_distribution, next_population)
            )

        next_frequency: list[float] = []
        for n, q, p in zip(next_population, q_next, mutated):
            n_eff = _effective_size(n, q, parameters)
            gene_copies = max(2, round(2.0 * n_eff))
            next_frequency.append(_binomial(rng, gene_copies, p) / gene_copies)

        population = tuple(next_population)
        interaction = q_next
        frequency = tuple(next_frequency)
        trait_distribution = next_trait_distribution
        trait_abundance = next_trait_abundance
        snapshots.append(_snapshot(generation, population, interaction, frequency, trait_distribution, trait_abundance, parameters))

    return SimulationResult(parameters, tuple(snapshots))


@contextmanager
def patched_h1_mutation_runner(mutation_rate: float) -> Iterator[None]:
    """Temporarily route finite H1 continuation code through the mutation closure.

    The existing H1 resolution and full-state-hold audits are parameterised around
    the module-level simulator.  This scoped patch preserves their predeclared
    endpoint/grid/seed logic while making the numerical closure explicit in the
    caller's manifest.  The original functions are restored even if an audit
    raises.
    """
    rate = validate_symmetric_allele_mutation_rate(mutation_rate)
    if rate == 0.0:
        yield
        return
    import causal_model.finite_h1_continuation_state_audit as continuation
    import causal_model.finite_h1_hysteresis_audit as hysteresis

    original_continuation = continuation.simulate
    original_hysteresis = hysteresis.simulate

    def runner(parameters: DynamicsParameters) -> SimulationResult:
        return simulate_with_symmetric_allele_mutation(parameters, mutation_rate=rate)

    continuation.simulate = runner
    hysteresis.simulate = runner
    try:
        yield
    finally:
        continuation.simulate = original_continuation
        hysteresis.simulate = original_hysteresis


# Imported lazily above in the legacy-style deterministic trait path.  Keeping
# the import after public definitions makes the dependency explicit without
# changing the existing simulator module.
from causal_model.multipatch_criticality_dynamics import _abundance_from_distribution
from math import exp
