"""Conditional network allele-floor theorem for the finite-bin simulator.

The simulator applies local selection, census-weighted migration, then finite
allele sampling. If every patch starts above a common allele floor and each
patch's interaction and population remain above declared lower bounds, local
selection maps the common floor to a new selected floor. The migration step is a
convex combination of values all above that selected floor, so migration cannot
reduce it. Only finite sampling can break the common floor.

This is a conditional finite-horizon theorem. It does not prove the required
interaction/population region, nor does it prove trait-bin persistence. Its
interaction lower-envelope uses the explicit theorem conditions
``high_interaction_benefit >= 0`` and ``selection_strength >= 0``.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, trait_fitness


@dataclass(frozen=True)
class NetworkAlleleFloorOneStep:
    """One-step common allele-floor envelope and sampling risk."""

    allele_floor: float
    interaction_lower_bound: float
    population_lower_bound: int
    selected_floor: float
    migrated_floor: float
    gene_copy_lower_bound: int
    per_patch_sampling_failure_upper_bound: float
    any_patch_sampling_failure_upper_bound: float
    deterministic_floor_preserved_before_sampling: bool


@dataclass(frozen=True)
class NetworkAlleleFloorHorizon:
    """Finite-horizon conditional lower probability for a common allele floor."""

    patches: int
    horizon: int
    one_step: NetworkAlleleFloorOneStep
    horizon_failure_upper_bound: float
    horizon_retention_probability_lower_bound: float
    certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _round_lower(value: float) -> int:
    """Safe lower bound for Python's half-even round."""
    return ceil(value - 0.5)


def _require_monotone_high_selection(parameters: DynamicsParameters) -> None:
    """Check the conditions that make q_min a high-fitness lower envelope."""
    if parameters.high_interaction_benefit < 0.0:
        raise ValueError("network theorem requires non-negative high_interaction_benefit")
    if parameters.selection_strength < 0.0:
        raise ValueError("network theorem requires non-negative selection_strength")


def selected_allele_floor(
    allele_floor: float,
    interaction_lower_bound: float,
    parameters: DynamicsParameters,
) -> float:
    """Lower-bound the simulator's selected high-allele frequency.

    The simulator uses a one-allele viability update. Under the stated
    non-negative interaction-benefit and selection-strength conditions, selected
    frequency is increasing in current frequency and interaction. Therefore
    lower bounds on both inputs give the displayed lower envelope.
    """
    _require_monotone_high_selection(parameters)
    p = _probability(allele_floor, "allele_floor")
    q = _probability(interaction_lower_bound, "interaction_lower_bound")
    margin = trait_fitness(1.0, q, parameters) - parameters.viability_threshold
    high_fitness = max(1e-12, 1.0 + parameters.selection_strength * margin)
    return p * high_fitness / (p * high_fitness + (1.0 - p))


def common_floor_migration_bound(
    selected_floor: float, migration_rate: float) -> float:
    """Return the exact common-floor bound after census-weighted migration.

    When every patch has selected frequency at least ``selected_floor``, their
    census-weighted mean is at least that floor. Hence the update

        (1-m) p_j + m mean(p)

    is also at least it for every migration rate in [0, 1].
    """
    selected_floor = _probability(selected_floor, "selected_floor")
    _probability(migration_rate, "migration_rate")
    return selected_floor


def _gene_copy_lower_bound(
    population_lower_bound: int,
    parameters: DynamicsParameters,
) -> int:
    if population_lower_bound < 1:
        raise ValueError("population_lower_bound must be at least one")
    # The simulator's effective size is N*effective_fraction*(1-skew*q).
    # For q in [0,1], q=1 is a safe lower envelope when skew_penalty >= 0.
    if parameters.skew_penalty < 0.0:
        raise ValueError("network theorem requires non-negative skew_penalty")
    effective_lower = max(
        1.0,
        parameters.effective_fraction
        * population_lower_bound
        * (1.0 - parameters.skew_penalty),
    )
    return max(2, _round_lower(2.0 * effective_lower))


def _chernoff_fraction_failure(
    trials_lower: int,
    success_probability_lower: float,
    target_floor: float,
) -> float:
    if trials_lower < 1:
        raise ValueError("trials_lower must be positive")
    pi = _probability(success_probability_lower, "success_probability_lower")
    threshold = _probability(target_floor, "target_floor")
    if threshold >= pi:
        return 1.0
    mu = trials_lower * pi
    delta = 1.0 - threshold / pi
    return min(1.0, exp(-mu * delta * delta / 2.0))


def network_allele_floor_one_step(
    parameters: DynamicsParameters,
    patches: int,
    allele_floor: float,
    interaction_lower_bound: float,
    population_lower_bound: int,
) -> NetworkAlleleFloorOneStep:
    """Certify a one-step common allele floor conditional on q/N lower bounds.

    The statement applies when every patch begins with frequency at least
    ``allele_floor``, has interaction at least ``interaction_lower_bound``, and
    its next sampling population is at least ``population_lower_bound``. It is
    exact for the selection/migration algebra of the declared simulator and
    probabilistic only at the finite sampling step.
    """
    if patches < 1:
        raise ValueError("patches must be at least one")
    floor = _probability(allele_floor, "allele_floor")
    selected = selected_allele_floor(floor, interaction_lower_bound, parameters)
    migrated = common_floor_migration_bound(selected, parameters.migration_rate)
    copies = _gene_copy_lower_bound(population_lower_bound, parameters)
    per_patch_failure = _chernoff_fraction_failure(copies, migrated, floor)
    any_failure = min(1.0, patches * per_patch_failure)
    return NetworkAlleleFloorOneStep(
        allele_floor=floor,
        interaction_lower_bound=float(interaction_lower_bound),
        population_lower_bound=int(population_lower_bound),
        selected_floor=selected,
        migrated_floor=migrated,
        gene_copy_lower_bound=copies,
        per_patch_sampling_failure_upper_bound=per_patch_failure,
        any_patch_sampling_failure_upper_bound=any_failure,
        deterministic_floor_preserved_before_sampling=migrated >= floor,
    )


def network_allele_floor_horizon(
    parameters: DynamicsParameters,
    patches: int,
    allele_floor: float,
    interaction_lower_bound: float,
    population_lower_bound: int,
    horizon: int,
) -> NetworkAlleleFloorHorizon:
    """Union-bound the conditional common-floor failure risk through a horizon."""
    if horizon < 1:
        raise ValueError("horizon must be at least one")
    one_step = network_allele_floor_one_step(
        parameters,
        patches,
        allele_floor,
        interaction_lower_bound,
        population_lower_bound,
    )
    failure = min(1.0, horizon * one_step.any_patch_sampling_failure_upper_bound)
    retention = max(0.0, 1.0 - failure) if one_step.deterministic_floor_preserved_before_sampling else 0.0
    return NetworkAlleleFloorHorizon(
        patches=patches,
        horizon=horizon,
        one_step=one_step,
        horizon_failure_upper_bound=failure,
        horizon_retention_probability_lower_bound=retention,
        certified=one_step.deterministic_floor_preserved_before_sampling and retention > 0.0,
    )
