"""Declared multi-patch simulator for eco-genetic criticality hypotheses.

This is intentionally a *simulation layer*. Its equations are documented in
``docs/multipatch_criticality_dynamics_contract.md`` and should not be read as a
proof of the general theorem layer. It exposes every life-cycle assumption needed
to examine H_critical, H_genetic_lag, and H_fragmentation under finite populations.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp
from random import Random
from typing import Sequence


@dataclass(frozen=True)
class DynamicsParameters:
    patch_areas: tuple[float, ...]
    generations: int = 30
    initial_population: tuple[int, ...] = ()
    initial_interaction: tuple[float, ...] = ()
    initial_high_allele_frequency: tuple[float, ...] = ()
    initial_trait_distribution: tuple[tuple[float, ...], ...] = ()
    initial_trait_abundance: tuple[tuple[int, ...], ...] = ()
    density_capacity: float = 40.0
    area_reference: float = 1.0
    interaction_feedback: float = 6.0
    interaction_memory_weight: float = 0.6
    interaction_barrier: float = 0.5
    low_base: float = 1.1
    low_cost: float = 0.8
    high_base: float = 0.2
    high_interaction_benefit: float = 0.8
    high_peak_width: float = 0.15
    viability_threshold: float = 1.0
    trait_grid_size: int = 101
    high_trait_cutoff: float = 0.7
    realised_high_trait_threshold: float = 1e-3
    realised_bin_abundance_threshold: int = 1
    realised_high_trait_abundance_threshold: int = 1
    trait_occupancy_mode: str = "deterministic_viability_selection"
    trait_occupancy_model: str = "viability_selection_local_recruitment"
    trait_selection_floor: float = 1e-12
    genotype_trait_recruitment: str = "resident_trait_only"
    inheritance_weight: float = 1.0
    low_trait_kernel_center: float = 0.0
    high_trait_kernel_center: float = 1.0
    trait_kernel_width: float = 0.2
    q_feedback_alpha: float | None = None
    q_feedback_beta_trait: float = 0.0
    q_feedback_gamma_allele: float | None = None
    selection_strength: float = 0.5
    baseline_growth: float = 0.3
    interaction_growth: float = 0.4
    high_allele_growth: float = 0.1
    effective_fraction: float = 0.6
    skew_penalty: float = 0.0
    migration_rate: float = 0.0
    random_seed: int = 1

    def __post_init__(self) -> None:
        if not self.patch_areas or any(area <= 0.0 for area in self.patch_areas):
            raise ValueError("patch_areas must be nonempty and positive")
        n = len(self.patch_areas)
        if self.generations < 1:
            raise ValueError("generations must be positive")
        for values, label, lower, upper in (
            (self.initial_population, "initial_population", 1, None),
            (self.initial_interaction, "initial_interaction", 0.0, 1.0),
            (self.initial_high_allele_frequency, "initial_high_allele_frequency", 0.0, 1.0),
        ):
            if values and len(values) != n:
                raise ValueError(f"{label} must be empty or match patch count")
            if upper is None:
                if values and any(value < lower for value in values):
                    raise ValueError(f"{label} has values below {lower}")
            elif values and any(value < lower or value > upper for value in values):
                raise ValueError(f"{label} must lie in [{lower}, {upper}]")
        if self.density_capacity <= 0.0 or self.area_reference <= 0.0:
            raise ValueError("density_capacity and area_reference must be positive")
        if self.interaction_feedback <= 0.0:
            raise ValueError("interaction_feedback must be positive")
        if not 0.0 <= self.interaction_memory_weight <= 1.0:
            raise ValueError("interaction_memory_weight must lie in [0, 1]")
        if self.high_peak_width <= 0.0:
            raise ValueError("high_peak_width must be positive")
        if self.trait_grid_size < 3:
            raise ValueError("trait_grid_size must be at least 3")
        if not 0.0 <= self.high_trait_cutoff <= 1.0:
            raise ValueError("high_trait_cutoff must lie in [0, 1]")
        if self.realised_high_trait_threshold < 0.0:
            raise ValueError("realised_high_trait_threshold must be nonnegative")
        if self.realised_bin_abundance_threshold < 1 or self.realised_high_trait_abundance_threshold < 1:
            raise ValueError("realised abundance thresholds must be positive integers")
        if self.trait_occupancy_mode not in {"deterministic_viability_selection", "finite_trait_bin_recruitment"}:
            raise ValueError("unknown trait_occupancy_mode")
        if self.trait_occupancy_model != "viability_selection_local_recruitment":
            raise ValueError("unknown trait_occupancy_model")
        if self.trait_selection_floor <= 0.0:
            raise ValueError("trait_selection_floor must be positive")
        if self.genotype_trait_recruitment not in {"resident_trait_only", "two_kernel_recruitment"}:
            raise ValueError("unknown genotype_trait_recruitment")
        if not 0.0 <= self.inheritance_weight <= 1.0:
            raise ValueError("inheritance_weight must lie in [0, 1]")
        if not 0.0 <= self.low_trait_kernel_center <= 1.0 or not 0.0 <= self.high_trait_kernel_center <= 1.0:
            raise ValueError("trait kernel centers must lie in [0, 1]")
        if self.trait_kernel_width <= 0.0:
            raise ValueError("trait_kernel_width must be positive")
        if self.q_feedback_alpha is not None and not 0.0 <= self.q_feedback_alpha <= 1.0:
            raise ValueError("q_feedback_alpha must lie in [0, 1]")
        if self.q_feedback_beta_trait < 0.0:
            raise ValueError("q_feedback_beta_trait must be nonnegative")
        if self.q_feedback_gamma_allele is not None and self.q_feedback_gamma_allele < 0.0:
            raise ValueError("q_feedback_gamma_allele must be nonnegative")
        if self.initial_trait_distribution:
            if len(self.initial_trait_distribution) != n:
                raise ValueError("initial_trait_distribution must be empty or match patch count")
            for distribution in self.initial_trait_distribution:
                if len(distribution) != self.trait_grid_size:
                    raise ValueError("each initial_trait_distribution row must match trait_grid_size")
                if any(value < 0.0 for value in distribution):
                    raise ValueError("initial_trait_distribution cannot contain negative values")
                if sum(distribution) <= 0.0:
                    raise ValueError("each initial_trait_distribution row must have positive mass")
        if self.initial_trait_abundance:
            if len(self.initial_trait_abundance) != n:
                raise ValueError("initial_trait_abundance must be empty or match patch count")
            for abundance in self.initial_trait_abundance:
                if len(abundance) != self.trait_grid_size:
                    raise ValueError("each initial_trait_abundance row must match trait_grid_size")
                if any(value < 0 for value in abundance):
                    raise ValueError("initial_trait_abundance cannot contain negative values")
                if sum(abundance) <= 0:
                    raise ValueError("each initial_trait_abundance row must have positive abundance")
        if not 0.0 < self.effective_fraction <= 1.0:
            raise ValueError("effective_fraction must lie in (0, 1]")
        if not 0.0 <= self.skew_penalty < 1.0:
            raise ValueError("skew_penalty must lie in [0, 1)")
        if not 0.0 <= self.migration_rate <= 1.0:
            raise ValueError("migration_rate must lie in [0, 1]")


@dataclass(frozen=True)
class TraitSpaceSummary:
    viable_components: int
    high_trait_component_present: bool
    high_trait_margin: float


@dataclass(frozen=True)
class TraitOccupancySummary:
    """Realised resident trait distribution summary for one patch.

    This is a simulation state, not a theorem. It is separate from the
    high-allele frequency p_j and from potential viability Omega_tau(q).
    """

    distribution: tuple[float, ...]
    abundance: tuple[int, ...]
    total_abundance: int
    high_trait_abundance: int
    high_trait_mass: float
    realised_high_trait_occupied: bool
    realised_components: int
    potential_high_trait_component_present: bool


@dataclass(frozen=True)
class SimulationSnapshot:
    generation: int
    interaction: tuple[float, ...]
    population: tuple[int, ...]
    effective_size: tuple[float, ...]
    high_allele_frequency: tuple[float, ...]
    trait_space: tuple[TraitSpaceSummary, ...]
    trait_occupancy: tuple[TraitOccupancySummary, ...]
    h_alpha: float
    h_gamma: float
    fst: float | None


@dataclass(frozen=True)
class SimulationResult:
    parameters: DynamicsParameters
    snapshots: tuple[SimulationSnapshot, ...]


@dataclass(frozen=True)
class FirstPassageEvent:
    name: str
    occurred: bool
    time: int | None
    censored: bool
    threshold: float | int | None
    aggregation_rule: str


def sigmoid(value: float) -> float:
    if value >= 0.0:
        inverse = exp(-value)
        return 1.0 / (1.0 + inverse)
    inverse = exp(value)
    return inverse / (1.0 + inverse)


def trait_fitness(z: float, interaction: float, parameters: DynamicsParameters) -> float:
    """Declared continuous trait-performance surface W(z;q)."""
    if not 0.0 <= z <= 1.0 or not 0.0 <= interaction <= 1.0:
        raise ValueError("z and interaction must lie in [0, 1]")
    low_route = parameters.low_base - parameters.low_cost * z * z
    peak = parameters.high_base + parameters.high_interaction_benefit * interaction
    high_route = peak * exp(-((z - 1.0) / parameters.high_peak_width) ** 2)
    return low_route + high_route


def trait_grid(parameters: DynamicsParameters) -> tuple[float, ...]:
    """Return the declared trait grid z_k in [0, 1]."""
    return tuple(index / (parameters.trait_grid_size - 1) for index in range(parameters.trait_grid_size))


def trait_space_summary(interaction: float, parameters: DynamicsParameters) -> TraitSpaceSummary:
    """Summarise viable-grid topology and high-investment component presence."""
    grid = trait_grid(parameters)
    viable = tuple(trait_fitness(z, interaction, parameters) >= parameters.viability_threshold for z in grid)
    components = 0
    in_component = False
    high_present = False
    for z, is_viable in zip(grid, viable):
        if is_viable and not in_component:
            components += 1
            in_component = True
        if not is_viable:
            in_component = False
        if is_viable and z >= parameters.high_trait_cutoff:
            high_present = True
    margin = max(
        trait_fitness(z, interaction, parameters) - parameters.viability_threshold
        for z in grid
        if z >= parameters.high_trait_cutoff
    )
    return TraitSpaceSummary(components, high_present, margin)


def trait_occupancy_summary(
    distribution: Sequence[float],
    potential: TraitSpaceSummary,
    parameters: DynamicsParameters,
    abundance: Sequence[int] | None = None,
) -> TraitOccupancySummary:
    """Summarise realised trait-bin occupancy independently of potential viability."""
    normalised = _normalise_distribution(distribution)
    grid = trait_grid(parameters)
    if abundance is None:
        counts = _abundance_from_distribution(normalised, 0)
    else:
        counts = tuple(int(value) for value in abundance)
    if parameters.trait_occupancy_mode == "finite_trait_bin_recruitment":
        occupied = tuple(value >= parameters.realised_bin_abundance_threshold for value in counts)
    else:
        occupied = tuple(value > 0.0 for value in normalised)
    components = 0
    in_component = False
    for is_occupied in occupied:
        if is_occupied and not in_component:
            components += 1
            in_component = True
        if not is_occupied:
            in_component = False
    high_abundance = sum(value for z, value in zip(grid, counts) if z >= parameters.high_trait_cutoff)
    total_abundance = sum(counts)
    high_mass = sum(value for z, value in zip(grid, normalised) if z >= parameters.high_trait_cutoff)
    if parameters.trait_occupancy_mode == "finite_trait_bin_recruitment":
        realised_high = high_abundance >= parameters.realised_high_trait_abundance_threshold
    else:
        realised_high = high_mass > parameters.realised_high_trait_threshold
    return TraitOccupancySummary(
        distribution=normalised,
        abundance=counts,
        total_abundance=total_abundance,
        high_trait_abundance=high_abundance,
        high_trait_mass=high_mass,
        realised_high_trait_occupied=realised_high,
        realised_components=components,
        potential_high_trait_component_present=potential.high_trait_component_present,
    )


def _heterozygosity(p: float) -> float:
    return 2.0 * p * (1.0 - p)


def _diversity(frequencies: Sequence[float], weights: Sequence[float]) -> tuple[float, float, float | None]:
    total = sum(weights)
    normalised = tuple(weight / total for weight in weights)
    p_bar = sum(weight * p for weight, p in zip(normalised, frequencies))
    h_alpha = sum(weight * _heterozygosity(p) for weight, p in zip(normalised, frequencies))
    h_gamma = _heterozygosity(p_bar)
    fst = None if h_gamma <= 0.0 else 1.0 - h_alpha / h_gamma
    return h_alpha, h_gamma, fst


def _binomial(rng: Random, trials: int, probability: float) -> int:
    probability = min(1.0, max(0.0, probability))
    return sum(rng.random() < probability for _ in range(trials))


def _multinomial(rng: Random, trials: int, probabilities: Sequence[float]) -> tuple[int, ...]:
    if trials < 0:
        raise ValueError("trials must be nonnegative")
    normalised = _normalise_distribution(probabilities)
    counts = [0 for _ in normalised]
    cumulative: list[float] = []
    running = 0.0
    for probability in normalised:
        running += probability
        cumulative.append(running)
    cumulative[-1] = 1.0
    for _ in range(trials):
        draw = rng.random()
        for index, boundary in enumerate(cumulative):
            if draw <= boundary:
                counts[index] += 1
                break
    return tuple(counts)


def _normalise_distribution(distribution: Sequence[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in distribution)
    total = sum(values)
    if total <= 0.0:
        raise ValueError("trait distribution must have positive mass")
    return tuple(value / total for value in values)


def _abundance_from_distribution(distribution: Sequence[float], total: int) -> tuple[int, ...]:
    if total <= 0:
        return tuple(0 for _ in distribution)
    normalised = _normalise_distribution(distribution)
    raw = [value * total for value in normalised]
    counts = [int(value) for value in raw]
    remainder = total - sum(counts)
    order = sorted(range(len(raw)), key=lambda index: raw[index] - counts[index], reverse=True)
    for index in order[:remainder]:
        counts[index] += 1
    return tuple(counts)


def _default_trait_distribution(parameters: DynamicsParameters) -> tuple[float, ...]:
    grid = trait_grid(parameters)
    weights = tuple(max(parameters.trait_selection_floor, 1.0 - z) for z in grid)
    return _normalise_distribution(weights)


def _kernel_distribution(center: float, parameters: DynamicsParameters) -> tuple[float, ...]:
    grid = trait_grid(parameters)
    high_kernel = center >= parameters.high_trait_cutoff
    weights = []
    for z in grid:
        in_declared_region = z >= parameters.high_trait_cutoff if high_kernel else z < parameters.high_trait_cutoff
        if in_declared_region:
            weights.append(max(parameters.trait_selection_floor, exp(-((z - center) / parameters.trait_kernel_width) ** 2)))
        else:
            weights.append(0.0)
    return _normalise_distribution(weights)


def recruit_trait_distribution(
    resident_distribution: Sequence[float],
    high_allele_frequency: float,
    parameters: DynamicsParameters,
) -> tuple[float, ...]:
    """Return the declared pre-selection recruit trait distribution.

    ``two_kernel_recruitment`` is an allele-linked recruitment closure, not a
    Mendelian inheritance model and not mutation across trait bins.
    """
    resident = _normalise_distribution(resident_distribution)
    if parameters.genotype_trait_recruitment == "resident_trait_only":
        return resident
    low_kernel = _kernel_distribution(parameters.low_trait_kernel_center, parameters)
    high_kernel = _kernel_distribution(parameters.high_trait_kernel_center, parameters)
    allele_kernel = tuple(
        (1.0 - high_allele_frequency) * low + high_allele_frequency * high
        for low, high in zip(low_kernel, high_kernel)
    )
    return _normalise_distribution(
        (1.0 - parameters.inheritance_weight) * kernel + parameters.inheritance_weight * mass
        for kernel, mass in zip(allele_kernel, resident)
    )


def _initial_values(
    parameters: DynamicsParameters,
) -> tuple[
    tuple[int, ...],
    tuple[float, ...],
    tuple[float, ...],
    tuple[tuple[float, ...], ...],
    tuple[tuple[int, ...], ...],
]:
    n = len(parameters.patch_areas)
    population = parameters.initial_population or tuple(
        max(1, round(parameters.density_capacity * area * 0.6)) for area in parameters.patch_areas
    )
    interaction = parameters.initial_interaction or tuple(0.5 for _ in range(n))
    frequency = parameters.initial_high_allele_frequency or tuple(0.5 for _ in range(n))
    if parameters.initial_trait_distribution:
        trait_distribution = tuple(_normalise_distribution(row) for row in parameters.initial_trait_distribution)
    else:
        trait_distribution = tuple(_default_trait_distribution(parameters) for _ in range(n))
    if parameters.initial_trait_abundance:
        trait_abundance = tuple(tuple(int(value) for value in row) for row in parameters.initial_trait_abundance)
        if not parameters.initial_population:
            population = tuple(sum(row) for row in trait_abundance)
        trait_distribution = tuple(_normalise_distribution(row) for row in trait_abundance)
    else:
        trait_abundance = tuple(
            _abundance_from_distribution(distribution, count)
            for distribution, count in zip(trait_distribution, population)
        )
    return tuple(population), tuple(interaction), tuple(frequency), trait_distribution, trait_abundance


def _effective_size(population: int, interaction: float, parameters: DynamicsParameters) -> float:
    # Interaction may increase demographic support but can also increase skew.
    return max(1.0, parameters.effective_fraction * population * (1.0 - parameters.skew_penalty * interaction))


def _snapshot(
    generation: int,
    population: tuple[int, ...],
    interaction: tuple[float, ...],
    frequency: tuple[float, ...],
    trait_distribution: tuple[tuple[float, ...], ...],
    trait_abundance: tuple[tuple[int, ...], ...],
    parameters: DynamicsParameters,
) -> SimulationSnapshot:
    effective = tuple(_effective_size(n, q, parameters) for n, q in zip(population, interaction))
    trait_spaces = tuple(trait_space_summary(q, parameters) for q in interaction)
    trait_occupancies = tuple(
        trait_occupancy_summary(mu, potential, parameters, abundance)
        for mu, potential, abundance in zip(trait_distribution, trait_spaces, trait_abundance)
    )
    h_alpha, h_gamma, fst = _diversity(frequency, tuple(float(n) for n in population))
    return SimulationSnapshot(
        generation,
        interaction,
        population,
        effective,
        frequency,
        trait_spaces,
        trait_occupancies,
        h_alpha,
        h_gamma,
        fst,
    )


def update_trait_distribution(
    distribution: Sequence[float],
    interaction: float,
    parameters: DynamicsParameters,
    high_allele_frequency: float = 0.0,
) -> tuple[float, ...]:
    """Update resident trait-bin frequencies by the declared simple model.

    Model: viability selection plus local recruitment. Trait-bin mass at z_k is
    multiplied by max(trait_selection_floor, W(z_k; q_t)) and then normalised.
    There is no mutation or dispersal across trait bins in this first extension.
    """
    grid = trait_grid(parameters)
    normalised = recruit_trait_distribution(distribution, high_allele_frequency, parameters)
    weighted = tuple(
        mass * max(parameters.trait_selection_floor, trait_fitness(z, interaction, parameters))
        for z, mass in zip(grid, normalised)
    )
    return _normalise_distribution(weighted)


def update_trait_abundance(
    abundance: Sequence[int],
    interaction: float,
    high_allele_frequency: float,
    next_population: int,
    parameters: DynamicsParameters,
    rng: Random,
) -> tuple[int, ...]:
    """Update realised trait-bin abundance under the finite recruitment closure."""
    if next_population < 1:
        raise ValueError("next_population must be positive")
    resident = _normalise_distribution(abundance)
    recruit = recruit_trait_distribution(resident, high_allele_frequency, parameters)
    grid = trait_grid(parameters)
    selected = _normalise_distribution(
        mass * max(parameters.trait_selection_floor, trait_fitness(z, interaction, parameters))
        for z, mass in zip(grid, recruit)
    )
    return _multinomial(rng, next_population, selected)


def _feedback_weights(parameters: DynamicsParameters) -> tuple[float, float, float]:
    if parameters.q_feedback_alpha is None and parameters.q_feedback_gamma_allele is None:
        return parameters.interaction_memory_weight, parameters.q_feedback_beta_trait, 1.0 - parameters.interaction_memory_weight
    alpha = parameters.interaction_memory_weight if parameters.q_feedback_alpha is None else parameters.q_feedback_alpha
    gamma = 0.0 if parameters.q_feedback_gamma_allele is None else parameters.q_feedback_gamma_allele
    return alpha, parameters.q_feedback_beta_trait, gamma


def interaction_support_signal(
    interaction: float,
    realised_high_trait_mass: float,
    high_allele_frequency: float,
    parameters: DynamicsParameters,
) -> float:
    alpha, beta_trait, gamma = _feedback_weights(parameters)
    return alpha * interaction + beta_trait * realised_high_trait_mass + gamma * high_allele_frequency


def simulate(parameters: DynamicsParameters) -> SimulationResult:
    """Run the declared finite-population multi-patch life cycle."""
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
                * (
                    (area / parameters.area_reference)
                    * dens
                    * signal
                    - parameters.interaction_barrier
                )
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
        total_weight = sum(weights)
        selected_mean = sum(weight * p for weight, p in zip(weights, selected)) / total_weight
        migrated = tuple(
            (1.0 - parameters.migration_rate) * p + parameters.migration_rate * selected_mean
            for p in selected
        )

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
        for n, q, p in zip(next_population, q_next, migrated):
            n_eff = _effective_size(n, q, parameters)
            gene_copies = max(2, round(2.0 * n_eff))
            next_frequency.append(_binomial(rng, gene_copies, p) / gene_copies)

        population = tuple(next_population)
        interaction = q_next
        frequency = tuple(next_frequency)
        trait_distribution = next_trait_distribution
        trait_abundance = next_trait_abundance
        snapshots.append(
            _snapshot(generation, population, interaction, frequency, trait_distribution, trait_abundance, parameters)
        )

    return SimulationResult(parameters, tuple(snapshots))


def first_high_trait_absence(result: SimulationResult) -> int | None:
    """Backward-compatible alias for first_potential_high_trait_absence()."""
    return first_potential_high_trait_absence(result)


def tau_trait_potential(result: SimulationResult) -> int | None:
    """Return tau_trait_potential: first loss of potential high-trait viability."""
    return first_potential_high_trait_absence(result)


def first_potential_high_trait_absence(result: SimulationResult) -> int | None:
    """Return first generation with all-patch loss of potential high-trait viability."""
    for snapshot in result.snapshots:
        if not any(summary.high_trait_component_present for summary in snapshot.trait_space):
            return snapshot.generation
    return None


def first_realised_high_trait_absence(result: SimulationResult) -> int | None:
    """Return first generation with all-patch loss of realised high-trait occupancy."""
    for snapshot in result.snapshots:
        if not any(summary.realised_high_trait_occupied for summary in snapshot.trait_occupancy):
            return snapshot.generation
    return None


def tau_trait_realised(result: SimulationResult) -> int | None:
    """Return tau_trait_realised: first loss of realised high-trait occupancy."""
    return first_realised_high_trait_absence(result)


def first_allele_loss(result: SimulationResult, allele_threshold: float = 0.0) -> int | None:
    """Return first generation where all patches are at or below allele_threshold."""
    if not 0.0 <= allele_threshold <= 1.0:
        raise ValueError("allele_threshold must lie in [0, 1]")
    for snapshot in result.snapshots:
        if all(p <= allele_threshold for p in snapshot.high_allele_frequency):
            return snapshot.generation
    return None


def tau_allele_loss(result: SimulationResult, allele_threshold: float = 0.0) -> int | None:
    """Return tau_allele_loss under the all_patch_loss aggregation rule."""
    return first_allele_loss(result, allele_threshold)


def first_alpha_warning(result: SimulationResult, warning_threshold: float) -> int | None:
    """Backward-compatible alias for first_h_alpha_warning()."""
    return first_h_alpha_warning(result, warning_threshold)


def first_h_alpha_warning(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return first generation where census-weighted H_alpha crosses a boundary."""
    if not 0.0 <= warning_threshold <= 1.0:
        raise ValueError("warning_threshold must lie in [0, 1]")
    for snapshot in result.snapshots:
        if snapshot.h_alpha <= warning_threshold:
            return snapshot.generation
    return None


def tau_H_alpha(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return tau_H_alpha: first H_alpha warning-threshold crossing."""
    return first_h_alpha_warning(result, warning_threshold)


def first_h_gamma_warning(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return first generation where census-weighted H_gamma crosses a boundary."""
    if not 0.0 <= warning_threshold <= 1.0:
        raise ValueError("warning_threshold must lie in [0, 1]")
    for snapshot in result.snapshots:
        if snapshot.h_gamma <= warning_threshold:
            return snapshot.generation
    return None


def tau_H_gamma(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return tau_H_gamma: first H_gamma warning-threshold crossing."""
    return first_h_gamma_warning(result, warning_threshold)


def first_fst_warning(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return first generation where F_ST reaches or exceeds a boundary."""
    if not 0.0 <= warning_threshold <= 1.0:
        raise ValueError("warning_threshold must lie in [0, 1]")
    for snapshot in result.snapshots:
        if snapshot.fst is not None and snapshot.fst >= warning_threshold:
            return snapshot.generation
    return None


def tau_FST(result: SimulationResult, warning_threshold: float) -> int | None:
    """Return tau_FST: first F_ST warning-threshold crossing."""
    return first_fst_warning(result, warning_threshold)


def first_passage_event(
    name: str,
    time: int | None,
    *,
    threshold: float | int | None,
    aggregation_rule: str,
) -> FirstPassageEvent:
    """Wrap a first-passage time with explicit censoring metadata."""
    return FirstPassageEvent(
        name=name,
        occurred=time is not None,
        time=time,
        censored=time is None,
        threshold=threshold,
        aggregation_rule=aggregation_rule,
    )


def first_passage_events(
    result: SimulationResult,
    *,
    h_alpha_threshold: float,
    h_gamma_threshold: float,
    fst_threshold: float,
    allele_threshold: float = 0.0,
) -> tuple[FirstPassageEvent, ...]:
    """Return predeclared first-passage events with explicit aggregation rules."""
    return (
        first_passage_event(
            "tau_trait_potential",
            tau_trait_potential(result),
            threshold=result.parameters.viability_threshold,
            aggregation_rule="all_patch_loss",
        ),
        first_passage_event(
            "tau_trait_realised",
            tau_trait_realised(result),
            threshold=result.parameters.realised_high_trait_abundance_threshold
            if result.parameters.trait_occupancy_mode == "finite_trait_bin_recruitment"
            else result.parameters.realised_high_trait_threshold,
            aggregation_rule="all_patch_loss",
        ),
        first_passage_event(
            "tau_allele_loss",
            tau_allele_loss(result, allele_threshold),
            threshold=allele_threshold,
            aggregation_rule="all_patch_loss",
        ),
        first_passage_event(
            "tau_H_alpha",
            tau_H_alpha(result, h_alpha_threshold),
            threshold=h_alpha_threshold,
            aggregation_rule="metapopulation_weighted_loss",
        ),
        first_passage_event(
            "tau_H_gamma",
            tau_H_gamma(result, h_gamma_threshold),
            threshold=h_gamma_threshold,
            aggregation_rule="metapopulation_weighted_loss",
        ),
        first_passage_event(
            "tau_FST",
            tau_FST(result, fst_threshold),
            threshold=fst_threshold,
            aggregation_rule="metapopulation_weighted_loss",
        ),
    )
