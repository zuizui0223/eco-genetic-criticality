"""Audit how a full-simulator trajectory departs from the canonical H1 theorem.

The canonical H1 theorem applies to a one-state map.  The full simulator adds
density, trait/allele feedback, finite inheritance, and migration.  This module
measures the resulting update residual rather than silently treating every
simulation result as a canonical bifurcation result.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from math import isfinite

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    SimulationResult,
    interaction_support_signal,
    sigmoid,
)


@dataclass(frozen=True)
class H1TheoremBoundaryAudit:
    """Patchwise canonical-update discrepancy and named departure conditions."""

    patch_count: int
    generations: int
    maximum_canonical_update_residual: float
    mean_canonical_update_residual: float
    maximum_density_deviation_from_one: float
    maximum_support_deviation_from_interaction: float
    interaction_trait_feedback_enabled: bool
    interaction_allele_feedback_enabled: bool
    migration_enabled: bool
    patchwise_canonical_update_certified: bool
    single_patch_canonical_theorem_limit_certified: bool
    departure_labels: tuple[str, ...]

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def _validate_tolerance(tolerance: float) -> float:
    tolerance = float(tolerance)
    if not isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("tolerance must be finite and non-negative")
    return tolerance


def audit_h1_theorem_boundary(
    result: SimulationResult,
    *,
    tolerance: float = 1e-12,
) -> H1TheoremBoundaryAudit:
    """Compare every simulated q update against its patchwise canonical update.

    For a patch with current interaction q, canonical H1 predicts

    ``sigmoid(kappa * ((area / area_reference) * q - barrier))``.

    The full simulator instead inserts its realised density and support signal.
    The returned residual is therefore zero precisely when those additions reduce
    to density one and support equal to q, up to numeric tolerance.  The
    stronger single-patch theorem-limit flag additionally requires no migration.
    """
    tolerance = _validate_tolerance(tolerance)
    parameters: DynamicsParameters = result.parameters
    if len(result.snapshots) < 2:
        raise ValueError("result must contain at least one transition")

    residuals: list[float] = []
    density_deviations: list[float] = []
    support_deviations: list[float] = []
    for current, following in zip(result.snapshots[:-1], result.snapshots[1:]):
        for patch, area in enumerate(parameters.patch_areas):
            carrying = parameters.density_capacity * area
            density = min(1.0, current.population[patch] / carrying)
            occupancy = current.trait_occupancy[patch]
            support = interaction_support_signal(
                current.interaction[patch],
                occupancy.high_trait_mass,
                current.high_allele_frequency[patch],
                parameters,
            )
            canonical_next = sigmoid(
                parameters.interaction_feedback
                * ((area / parameters.area_reference) * current.interaction[patch] - parameters.interaction_barrier)
            )
            residuals.append(abs(following.interaction[patch] - canonical_next))
            density_deviations.append(abs(density - 1.0))
            support_deviations.append(abs(support - current.interaction[patch]))

    max_residual = max(residuals)
    labels: list[str] = []
    if max(density_deviations) > tolerance:
        labels.append("density_not_one")
    if parameters.q_feedback_beta_trait != 0.0:
        labels.append("trait_feedback_enabled")
    _, _, allele_weight = _feedback_weights_for_audit(parameters)
    if allele_weight != 0.0:
        labels.append("allele_feedback_enabled")
    if max(support_deviations) > tolerance:
        labels.append("support_not_equal_interaction")
    if parameters.migration_rate != 0.0:
        labels.append("migration_enabled")
    if len(parameters.patch_areas) != 1:
        labels.append("multiple_patches")

    patchwise = max_residual <= tolerance
    single_patch = patchwise and len(parameters.patch_areas) == 1 and parameters.migration_rate == 0.0
    return H1TheoremBoundaryAudit(
        patch_count=len(parameters.patch_areas),
        generations=parameters.generations,
        maximum_canonical_update_residual=max_residual,
        mean_canonical_update_residual=sum(residuals) / len(residuals),
        maximum_density_deviation_from_one=max(density_deviations),
        maximum_support_deviation_from_interaction=max(support_deviations),
        interaction_trait_feedback_enabled=parameters.q_feedback_beta_trait != 0.0,
        interaction_allele_feedback_enabled=allele_weight != 0.0,
        migration_enabled=parameters.migration_rate != 0.0,
        patchwise_canonical_update_certified=patchwise,
        single_patch_canonical_theorem_limit_certified=single_patch,
        departure_labels=tuple(labels),
    )


def _feedback_weights_for_audit(parameters: DynamicsParameters) -> tuple[float, float, float]:
    """Mirror the full simulator's declared q-support coefficients."""
    if parameters.q_feedback_alpha is None and parameters.q_feedback_gamma_allele is None:
        return (
            parameters.interaction_memory_weight,
            parameters.q_feedback_beta_trait,
            1.0 - parameters.interaction_memory_weight,
        )
    alpha = parameters.interaction_memory_weight if parameters.q_feedback_alpha is None else parameters.q_feedback_alpha
    gamma = 0.0 if parameters.q_feedback_gamma_allele is None else parameters.q_feedback_gamma_allele
    return alpha, parameters.q_feedback_beta_trait, gamma
