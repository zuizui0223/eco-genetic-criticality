"""Exact one-step recruitment factorisation for the colonization life cycle.

This module is a deliberately narrow bridge from the abstract channel theorems
to the colonization backend. It does **not** factorise the backend's multi-step,
stochastic invasion growth rate ``lambda``. Instead it exposes the expected
number of juvenile recruits retained at the end of one time step per initial
adult, conditional on a declared local context.

The life-cycle arithmetic mirrors ``_colonization_step``:

    adult survives
    -> adult conceives
    -> offspring either disperses and settles through a corridor,
       or remains local and settles if local room is available
    -> occupied patch survives the end-of-step extinction draw.

For a trait z, the expected retained juvenile recruitment is exactly

    W_recruit(z) = F_local(z) * E_settlement(z),

where

    F_local = survival_probability * conception_probability
    E_settlement = [1 - extinction_rate]
                   * {d(z) * connectivity * expected_target_room
                      + [1-d(z)] * local_room}.

``d(z)`` is the same ``benefit_shape`` dispersal-investment probability used by
the colonization IBM. This is an exact expectation for the specified one-step
context, not an approximation to long-run invasion lambda.
"""
from __future__ import annotations

from dataclasses import dataclass

from causal_model.colonization_metapopulation_abm import (
    ColonizationParameters,
    ColonizationRegime,
)
from causal_model.spatial_metapopulation_abm import _clip, benefit_shape


@dataclass(frozen=True)
class ColonizationRecruitmentContext:
    """Deterministic inputs for one parent's expected juvenile recruitment.

    ``mate_success`` and ``expected_target_room`` are context summaries rather
    than random draws. ``expected_target_room`` is the mean settlement room across
    a uniformly sampled reachable target, matching the target-choice rule in the
    colonization simulator.
    """

    trait: float
    age: int
    local_density: float
    resource: float
    mate_success: float
    local_room: float
    expected_target_room: float
    corridor_available: bool = True

    def __post_init__(self) -> None:
        if not 0.0 <= self.trait <= 1.0:
            raise ValueError("trait must lie in [0, 1]")
        if self.age < 0:
            raise ValueError("age must be non-negative")
        for name in (
            "local_density",
            "resource",
            "mate_success",
            "local_room",
            "expected_target_room",
        ):
            value = getattr(self, name)
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")


@dataclass(frozen=True)
class ColonizationRecruitmentFactors:
    """Declared factors for expected one-step retained juvenile recruitment."""

    trait: float
    survival_probability: float
    conception_probability: float
    dispersal_probability: float
    patch_persistence_probability: float
    local_reproductive_factor: float
    settlement_factor: float
    expected_juvenile_recruitment: float

    @property
    def factorisation_residual(self) -> float:
        return self.expected_juvenile_recruitment - (
            self.local_reproductive_factor * self.settlement_factor
        )

    @property
    def is_strictly_positive(self) -> bool:
        return self.local_reproductive_factor > 0.0 and self.settlement_factor > 0.0


def one_step_recruitment_factors(
    context: ColonizationRecruitmentContext,
    params: ColonizationParameters,
    regime: ColonizationRegime,
) -> ColonizationRecruitmentFactors:
    """Return the exact end-of-step recruitment factorisation for one adult.

    The formula follows the stochastic-event order in ``_colonization_step``.
    A dispersing offspring has no fallback to local settlement when the corridor
    attempt fails, so the expected settlement factor is a weighted sum of distinct
    dispersal and local branches. The global patch-extinction draw is independent
    of trait and location in the current IBM, and therefore multiplies the
    settlement factor by ``1-extinction_rate`` in expectation. At
    ``age >= max_age`` the simulator removes the adult before conception; the
    matching survival probability is zero.
    """
    age_term = (context.age / max(params.max_age, 1)) ** 2
    survival_probability = (
        0.0
        if context.age >= params.max_age
        else _clip(params.base_survival * (1.0 - 0.6 * age_term))
    )
    conception_probability = _clip(
        params.fecundity
        + 0.20 * context.mate_success
        + 0.30 * context.resource
        + regime.repro_baseline
        - params.dispersal_cost * context.trait
    )
    dispersal_probability = benefit_shape(context.trait, params.benefit_saturation)
    connectivity = _clip(regime.connectivity_present) if context.corridor_available else 0.0
    patch_persistence_probability = 1.0 - _clip(params.extinction_rate)
    pre_extinction_settlement = (
        dispersal_probability * connectivity * context.expected_target_room
        + (1.0 - dispersal_probability) * context.local_room
    )
    settlement_factor = patch_persistence_probability * pre_extinction_settlement
    local_reproductive_factor = survival_probability * conception_probability
    expected_juvenile_recruitment = local_reproductive_factor * settlement_factor
    return ColonizationRecruitmentFactors(
        trait=context.trait,
        survival_probability=survival_probability,
        conception_probability=conception_probability,
        dispersal_probability=dispersal_probability,
        patch_persistence_probability=patch_persistence_probability,
        local_reproductive_factor=local_reproductive_factor,
        settlement_factor=settlement_factor,
        expected_juvenile_recruitment=expected_juvenile_recruitment,
    )


def require_theorem_interior(factors: ColonizationRecruitmentFactors) -> ColonizationRecruitmentFactors:
    """Require the strict positivity assumed by N1--N4.

    Boundary states with zero conception or zero settlement remain valid biology
    and valid IBM states, but division-based channel-identifiability theorems do
    not apply there.
    """
    if not factors.is_strictly_positive:
        raise ValueError(
            "N1--N4 require positive F and E; this one-step context lies on a boundary"
        )
    return factors
