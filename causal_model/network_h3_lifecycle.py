"""Finite-population H3 lifecycle with individual dispersal and recolonisation.

This is a declared stochastic closure for the part of H3 that cannot be
represented by allele-frequency mixing alone.  It carries three finite states
per patch:

- census population size;
- realised high-trait abundance; and
- diploid high-allele copy count.

Each generation performs adult survival, individual emigration, directed
migration through a source-to-destination kernel, thresholded local persistence
or recolonisation, and density-limited recruitment.  It is intentionally a
model-specific simulation layer: the update makes no universal claim about
which networks rescue real populations.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from random import Random
from typing import Literal, Sequence


PatchStatus = Literal["persisted", "rescued", "extinct", "recolonised", "empty"]


@dataclass(frozen=True)
class PatchState:
    """Finite state of one diploid patch at a census time."""

    population: int
    high_trait_abundance: int
    high_allele_copies: int

    def __post_init__(self) -> None:
        if any(not isinstance(value, int) for value in (
            self.population,
            self.high_trait_abundance,
            self.high_allele_copies,
        )):
            raise ValueError("patch state values must be integers")
        if self.population < 0:
            raise ValueError("population must be non-negative")
        if not 0 <= self.high_trait_abundance <= self.population:
            raise ValueError("high_trait_abundance must lie between zero and population")
        if not 0 <= self.high_allele_copies <= 2 * self.population:
            raise ValueError("high_allele_copies must lie between zero and 2*population")

    @property
    def high_allele_frequency(self) -> float | None:
        """Return the allele frequency, or None for an extinct patch."""
        if self.population == 0:
            return None
        return self.high_allele_copies / (2.0 * self.population)


@dataclass(frozen=True)
class NetworkLifecycleParameters:
    """Declared life-cycle and landscape parameters for finite H3 experiments.

    ``source_to_destination_kernel[source][destination]`` is the conditional
    destination distribution of an emigrant from ``source``.  Each source row is
    non-negative and sums to one.  This individual transport kernel is distinct
    from the destination-by-source frequency-mixing matrix used in the H3
    allele-floor theorem.
    """

    capacities: tuple[int, ...]
    source_to_destination_kernel: tuple[tuple[float, ...], ...]
    generations: int = 20
    adult_survival_probability: float = 0.8
    emigration_probability: float = 0.2
    recruitment_per_adult: float = 0.5
    high_trait_recruitment_multiplier: float = 1.0
    persistence_threshold: int = 1
    colonisation_threshold: int = 1
    random_seed: int = 1

    def __post_init__(self) -> None:
        if not self.capacities or any(not isinstance(value, int) or value < 1 for value in self.capacities):
            raise ValueError("capacities must be nonempty positive integers")
        patch_count = len(self.capacities)
        if self.generations < 1:
            raise ValueError("generations must be positive")
        for name, value in (
            ("adult_survival_probability", self.adult_survival_probability),
            ("emigration_probability", self.emigration_probability),
        ):
            _validate_probability(value, name)
        for name, value in (
            ("recruitment_per_adult", self.recruitment_per_adult),
            ("high_trait_recruitment_multiplier", self.high_trait_recruitment_multiplier),
        ):
            if not isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and non-negative")
        if self.persistence_threshold < 1 or self.colonisation_threshold < 1:
            raise ValueError("persistence and colonisation thresholds must be positive")
        if self.persistence_threshold > max(self.capacities):
            raise ValueError("persistence_threshold cannot exceed every patch capacity")
        _validate_source_to_destination_kernel(self.source_to_destination_kernel, patch_count)


@dataclass(frozen=True)
class PatchTransition:
    """One patch's demographic, trait, and allele transition in one generation."""

    status: PatchStatus
    resident_adults_after_emigration: int
    inbound_migrants: int
    candidate_adults: int
    high_trait_adults_after_dispersal: int
    high_allele_copies_after_dispersal: int


@dataclass(frozen=True)
class NetworkSnapshot:
    """State and transition record for one generation."""

    generation: int
    states: tuple[PatchState, ...]
    transitions: tuple[PatchTransition, ...] | None
    h_alpha: float
    h_gamma: float
    fst: float | None


@dataclass(frozen=True)
class NetworkLifecycleResult:
    """Full reproducible trajectory for the declared H3 closure."""

    parameters: NetworkLifecycleParameters
    snapshots: tuple[NetworkSnapshot, ...]


def _validate_probability(value: float, name: str) -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _validate_source_to_destination_kernel(
    kernel: Sequence[Sequence[float]],
    patch_count: int,
    *,
    tolerance: float = 1e-12,
) -> tuple[tuple[float, ...], ...]:
    """Validate the source-row stochastic individual transport kernel."""
    rows = tuple(tuple(float(value) for value in row) for row in kernel)
    if len(rows) != patch_count or any(len(row) != patch_count for row in rows):
        raise ValueError("source_to_destination_kernel must be square and match capacities")
    for row in rows:
        if any(not isfinite(value) or value < -tolerance for value in row):
            raise ValueError("migration-kernel weights must be finite and non-negative")
        if abs(sum(row) - 1.0) > tolerance:
            raise ValueError("each source migration-kernel row must sum to one")
    return tuple(tuple(max(0.0, value) for value in row) for row in rows)


def isolated_kernel(patch_count: int) -> tuple[tuple[float, ...], ...]:
    """Return an individual-transport kernel with no between-patch dispersal."""
    if patch_count < 1:
        raise ValueError("patch_count must be at least one")
    return tuple(
        tuple(1.0 if source == destination else 0.0 for destination in range(patch_count))
        for source in range(patch_count)
    )


def _binomial(rng: Random, trials: int, probability: float) -> int:
    if trials < 0:
        raise ValueError("trials must be non-negative")
    probability = min(1.0, max(0.0, probability))
    return sum(rng.random() < probability for _ in range(trials))


def _multinomial(rng: Random, trials: int, probabilities: Sequence[float]) -> tuple[int, ...]:
    if trials < 0:
        raise ValueError("trials must be non-negative")
    values = tuple(float(value) for value in probabilities)
    if not values or any(value < 0.0 or not isfinite(value) for value in values):
        raise ValueError("probabilities must be nonempty, finite, and non-negative")
    total = sum(values)
    if total <= 0.0:
        raise ValueError("probabilities must have positive total")
    counts = [0 for _ in values]
    if trials == 0:
        return tuple(counts)
    cumulative: list[float] = []
    running = 0.0
    for value in values:
        running += value / total
        cumulative.append(running)
    cumulative[-1] = 1.0
    for _ in range(trials):
        draw = rng.random()
        for index, boundary in enumerate(cumulative):
            if draw <= boundary:
                counts[index] += 1
                break
    return tuple(counts)


def _state_diversity(states: Sequence[PatchState]) -> tuple[float, float, float | None]:
    total_population = sum(state.population for state in states)
    if total_population == 0:
        return 0.0, 0.0, None
    weights = tuple(state.population / total_population for state in states)
    frequencies = tuple(
        0.0 if state.population == 0 else state.high_allele_copies / (2.0 * state.population)
        for state in states
    )
    h_alpha = sum(weight * 2.0 * frequency * (1.0 - frequency) for weight, frequency in zip(weights, frequencies))
    pooled = sum(weight * frequency for weight, frequency in zip(weights, frequencies))
    h_gamma = 2.0 * pooled * (1.0 - pooled)
    fst = None if h_gamma <= 0.0 else 1.0 - h_alpha / h_gamma
    return h_alpha, h_gamma, fst


def _snapshot(
    generation: int,
    states: Sequence[PatchState],
    transitions: Sequence[PatchTransition] | None,
) -> NetworkSnapshot:
    h_alpha, h_gamma, fst = _state_diversity(states)
    return NetworkSnapshot(
        generation=generation,
        states=tuple(states),
        transitions=None if transitions is None else tuple(transitions),
        h_alpha=h_alpha,
        h_gamma=h_gamma,
        fst=fst,
    )


def _cap_state(
    population: int,
    high_trait_abundance: int,
    high_allele_copies: int,
    capacity: int,
    rng: Random,
) -> PatchState:
    """Thin a candidate state to capacity while retaining composition in expectation."""
    if population <= capacity:
        return PatchState(population, high_trait_abundance, high_allele_copies)
    high_probability = high_trait_abundance / population
    allele_probability = high_allele_copies / (2.0 * population)
    return PatchState(
        population=capacity,
        high_trait_abundance=_binomial(rng, capacity, high_probability),
        high_allele_copies=_binomial(rng, 2 * capacity, allele_probability),
    )


def _recruit(adults: PatchState, parameters: NetworkLifecycleParameters, capacity: int, rng: Random) -> PatchState:
    """Apply density-limited recruitment after survival, dispersal, and thresholding."""
    if adults.population == 0 or adults.population >= capacity:
        return adults
    expected_births = parameters.recruitment_per_adult * adults.population * (1.0 - adults.population / capacity)
    births = min(capacity - adults.population, max(0, round(expected_births)))
    if births == 0:
        return adults
    high_weight = adults.high_trait_abundance * parameters.high_trait_recruitment_multiplier
    low_weight = adults.population - adults.high_trait_abundance
    high_probability = 0.0 if high_weight + low_weight <= 0.0 else high_weight / (high_weight + low_weight)
    allele_probability = adults.high_allele_copies / (2.0 * adults.population)
    return PatchState(
        population=adults.population + births,
        high_trait_abundance=adults.high_trait_abundance + _binomial(rng, births, high_probability),
        high_allele_copies=adults.high_allele_copies + _binomial(rng, 2 * births, allele_probability),
    )


def _one_generation(
    states: Sequence[PatchState],
    parameters: NetworkLifecycleParameters,
    rng: Random,
) -> tuple[tuple[PatchState, ...], tuple[PatchTransition, ...]]:
    patch_count = len(states)
    resident_population = [0 for _ in range(patch_count)]
    resident_high_trait = [0 for _ in range(patch_count)]
    resident_allele_copies = [0 for _ in range(patch_count)]
    inbound_population = [0 for _ in range(patch_count)]
    inbound_high_trait = [0 for _ in range(patch_count)]
    inbound_allele_copies = [0 for _ in range(patch_count)]

    for source, state in enumerate(states):
        high_survivors = _binomial(rng, state.high_trait_abundance, parameters.adult_survival_probability)
        low_survivors = _binomial(
            rng,
            state.population - state.high_trait_abundance,
            parameters.adult_survival_probability,
        )
        survivor_population = high_survivors + low_survivors
        source_allele_frequency = 0.0 if state.population == 0 else state.high_allele_copies / (2.0 * state.population)

        high_emigrants = _binomial(rng, high_survivors, parameters.emigration_probability)
        low_emigrants = _binomial(rng, low_survivors, parameters.emigration_probability)
        emigrant_population = high_emigrants + low_emigrants
        resident_population[source] = survivor_population - emigrant_population
        resident_high_trait[source] = high_survivors - high_emigrants
        resident_allele_copies[source] = _binomial(
            rng,
            2 * resident_population[source],
            source_allele_frequency,
        )

        high_destinations = _multinomial(
            rng,
            high_emigrants,
            parameters.source_to_destination_kernel[source],
        )
        low_destinations = _multinomial(
            rng,
            low_emigrants,
            parameters.source_to_destination_kernel[source],
        )
        for destination, (incoming_high, incoming_low) in enumerate(zip(high_destinations, low_destinations)):
            incoming_total = incoming_high + incoming_low
            inbound_population[destination] += incoming_total
            inbound_high_trait[destination] += incoming_high
            inbound_allele_copies[destination] += _binomial(
                rng,
                2 * incoming_total,
                source_allele_frequency,
            )

    next_states: list[PatchState] = []
    transitions: list[PatchTransition] = []
    for patch, previous_state in enumerate(states):
        candidate_population = resident_population[patch] + inbound_population[patch]
        candidate_high_trait = resident_high_trait[patch] + inbound_high_trait[patch]
        candidate_allele_copies = resident_allele_copies[patch] + inbound_allele_copies[patch]
        residents = resident_population[patch]
        immigrants = inbound_population[patch]

        if candidate_population < parameters.persistence_threshold:
            status: PatchStatus = "empty" if previous_state.population == 0 else "extinct"
            adults = PatchState(0, 0, 0)
        elif previous_state.population == 0:
            if immigrants >= parameters.colonisation_threshold:
                status = "recolonised"
                adults = _cap_state(
                    candidate_population,
                    candidate_high_trait,
                    candidate_allele_copies,
                    parameters.capacities[patch],
                    rng,
                )
            else:
                status = "empty"
                adults = PatchState(0, 0, 0)
        else:
            adults = _cap_state(
                candidate_population,
                candidate_high_trait,
                candidate_allele_copies,
                parameters.capacities[patch],
                rng,
            )
            status = "rescued" if residents < parameters.persistence_threshold and immigrants > 0 else "persisted"

        next_state = _recruit(adults, parameters, parameters.capacities[patch], rng)
        next_states.append(next_state)
        transitions.append(
            PatchTransition(
                status=status,
                resident_adults_after_emigration=residents,
                inbound_migrants=immigrants,
                candidate_adults=candidate_population,
                high_trait_adults_after_dispersal=0 if adults.population == 0 else adults.high_trait_abundance,
                high_allele_copies_after_dispersal=0 if adults.population == 0 else adults.high_allele_copies,
            )
        )
    return tuple(next_states), tuple(transitions)


def simulate_network_lifecycle(
    parameters: NetworkLifecycleParameters,
    initial_states: Sequence[PatchState],
) -> NetworkLifecycleResult:
    """Run the declared H3 lifecycle from explicit finite patch states."""
    states = tuple(initial_states)
    if len(states) != len(parameters.capacities):
        raise ValueError("initial_states length must match capacities")
    if any(state.population > capacity for state, capacity in zip(states, parameters.capacities)):
        raise ValueError("initial patch population cannot exceed capacity")

    rng = Random(parameters.random_seed)
    snapshots = [_snapshot(0, states, None)]
    current = states
    for generation in range(1, parameters.generations + 1):
        current, transitions = _one_generation(current, parameters, rng)
        snapshots.append(_snapshot(generation, current, transitions))
    return NetworkLifecycleResult(parameters=parameters, snapshots=tuple(snapshots))


def first_metapopulation_extinction(result: NetworkLifecycleResult) -> int | None:
    """Return the first generation with no occupied patch."""
    for snapshot in result.snapshots:
        if all(state.population == 0 for state in snapshot.states):
            return snapshot.generation
    return None


def first_realised_high_trait_loss(result: NetworkLifecycleResult) -> int | None:
    """Return the first generation with no realised high-trait individual."""
    for snapshot in result.snapshots:
        if all(state.high_trait_abundance == 0 for state in snapshot.states):
            return snapshot.generation
    return None


def first_recolonisation(result: NetworkLifecycleResult) -> int | None:
    """Return the first generation in which at least one empty patch is recolonised."""
    for snapshot in result.snapshots[1:]:
        if snapshot.transitions is not None and any(transition.status == "recolonised" for transition in snapshot.transitions):
            return snapshot.generation
    return None
