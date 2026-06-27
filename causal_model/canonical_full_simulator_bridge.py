"""Exact embedding of the canonical H1 interaction map in the full simulator.

The full multipatch simulator is not generally the canonical H1 map.  This
module documents a strict one-patch parameter limit in which its interaction
trajectory is identical to

    q[t+1] = sigmoid(kappa * ((A/A_ref) q[t] - theta)).

The embedding fixes density at one, removes trait/allele feedback from the
interaction-support signal, and keeps census abundance exactly at carrying
capacity.  Trait and allele states may still be simulated, but they cannot
alter q in this declared limit.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, isfinite

from causal_model.canonical_h1_bifurcation import iterate_canonical_map
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate


@dataclass(frozen=True)
class CanonicalFullSimulatorBridgeCertificate:
    """Numerical equality check between canonical and full-simulator q trajectories."""

    area: float
    area_reference: float
    feedback_strength: float
    barrier: float
    initial_interaction: float
    generations: int
    carrying_population: int
    canonical_interaction: tuple[float, ...]
    full_simulator_interaction: tuple[float, ...]
    maximum_absolute_error: float
    exact_embedding_certified: bool


def canonical_full_simulator_parameters(
    *,
    area: float,
    area_reference: float,
    feedback_strength: float,
    barrier: float,
    initial_interaction: float,
    generations: int,
    carrying_population: int = 100,
    random_seed: int = 1,
) -> DynamicsParameters:
    """Return the declared one-patch full-simulator embedding of canonical H1.

    The choices are algebraic, not fitted:

    - `density_capacity=carrying_population/area` and initial census equal to
      `carrying_population` make density exactly one;
    - baseline growth one and all other growth terms zero retain that census;
    - `q_feedback_alpha=1`, trait feedback zero, and allele feedback zero make
      the support signal exactly current interaction.
    """
    for name, value in (
        ("area", area),
        ("area_reference", area_reference),
        ("feedback_strength", feedback_strength),
        ("barrier", barrier),
        ("initial_interaction", initial_interaction),
    ):
        if not isfinite(value):
            raise ValueError(f"{name} must be finite")
    if area <= 0.0 or area_reference <= 0.0 or feedback_strength <= 0.0:
        raise ValueError("area, area_reference, and feedback_strength must be positive")
    if not 0.0 <= initial_interaction <= 1.0:
        raise ValueError("initial_interaction must lie in [0, 1]")
    if not isinstance(generations, int) or generations < 1:
        raise ValueError("generations must be a positive integer")
    if not isinstance(carrying_population, int) or carrying_population < 1:
        raise ValueError("carrying_population must be a positive integer")

    return DynamicsParameters(
        patch_areas=(float(area),),
        generations=generations,
        initial_population=(carrying_population,),
        initial_interaction=(float(initial_interaction),),
        initial_high_allele_frequency=(0.5,),
        density_capacity=carrying_population / float(area),
        area_reference=float(area_reference),
        interaction_feedback=float(feedback_strength),
        interaction_barrier=float(barrier),
        interaction_memory_weight=1.0,
        q_feedback_alpha=1.0,
        q_feedback_beta_trait=0.0,
        q_feedback_gamma_allele=0.0,
        baseline_growth=1.0,
        interaction_growth=0.0,
        high_allele_growth=0.0,
        migration_rate=0.0,
        random_seed=random_seed,
    )


def canonical_full_simulator_bridge_certificate(
    *,
    area: float,
    area_reference: float = 1.0,
    feedback_strength: float,
    barrier: float,
    initial_interaction: float,
    generations: int = 50,
    carrying_population: int = 100,
    random_seed: int = 1,
    tolerance: float = 1e-12,
) -> CanonicalFullSimulatorBridgeCertificate:
    """Run the embedding and certify trajectory equality within numeric tolerance."""
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    parameters = canonical_full_simulator_parameters(
        area=area,
        area_reference=area_reference,
        feedback_strength=feedback_strength,
        barrier=barrier,
        initial_interaction=initial_interaction,
        generations=generations,
        carrying_population=carrying_population,
        random_seed=random_seed,
    )
    canonical = iterate_canonical_map(
        initial_interaction,
        feedback_strength=feedback_strength,
        area=area,
        area_reference=area_reference,
        barrier=barrier,
        iterations=generations,
    ).values
    full = simulate(parameters)
    trajectory = tuple(snapshot.interaction[0] for snapshot in full.snapshots)
    maximum_error = max(abs(left - right) for left, right in zip(canonical, trajectory))
    return CanonicalFullSimulatorBridgeCertificate(
        area=area,
        area_reference=area_reference,
        feedback_strength=feedback_strength,
        barrier=barrier,
        initial_interaction=initial_interaction,
        generations=generations,
        carrying_population=carrying_population,
        canonical_interaction=canonical,
        full_simulator_interaction=trajectory,
        maximum_absolute_error=maximum_error,
        exact_embedding_certified=maximum_error <= tolerance,
    )
