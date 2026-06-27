"""Exact bifurcation bookkeeping for the canonical H1 interaction map.

The canonical reduction is

    q[t+1] = sigmoid(kappa * ((A / A_ref) * q[t] - theta)).

It is deliberately narrower than the finite-bin multipatch simulator.  Here,
``q`` is the only feedback state, density is held at one, and no realised-trait
or allele terms enter the interaction update.  Within that declared map the
module gives an analytic bistability certificate, numerical fixed points, and
continuation traces for hysteresis checks.

Nothing here establishes bistability in every positive-feedback system or in
the full coupled simulator.  It establishes the canonical H1 mechanism under
its explicit reduction.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log, nextafter, sqrt
from typing import Literal, Sequence

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    trait_space_summary,
)


Stability = Literal["stable", "unstable", "neutral"]
BifurcationRegime = Literal[
    "monostable",
    "lower_saddle_node",
    "bistable",
    "upper_saddle_node",
]


@dataclass(frozen=True)
class CanonicalFixedPoint:
    """One fixed point of the declared canonical interaction map."""

    interaction: float
    local_multiplier: float
    stability: Stability


@dataclass(frozen=True)
class CanonicalBifurcationCertificate:
    """Analytic and numerical status of the canonical map at one barrier."""

    feedback_strength: float
    area: float
    area_reference: float
    barrier: float
    gain: float
    turning_points: tuple[float, float] | None
    bistable_barrier_interval: tuple[float, float] | None
    regime: BifurcationRegime
    fixed_points: tuple[CanonicalFixedPoint, ...]
    strict_bistability_certified: bool


@dataclass(frozen=True)
class HighTraitBranch:
    """High-trait viability margin evaluated on one stable interaction branch."""

    interaction: float
    high_trait_margin: float
    high_trait_component_present: bool


@dataclass(frozen=True)
class CanonicalH1Certificate:
    """Specified-system H1 certificate: bistability plus a trait-margin sign change."""

    bifurcation: CanonicalBifurcationCertificate
    low_stable_branch: HighTraitBranch | None
    high_stable_branch: HighTraitBranch | None
    branch_dependent_high_trait_mode: bool


@dataclass(frozen=True)
class CanonicalOrbit:
    """A deterministic orbit from one declared initial interaction state."""

    initial_interaction: float
    values: tuple[float, ...]

    @property
    def terminal_interaction(self) -> float:
        return self.values[-1]


@dataclass(frozen=True)
class BarrierContinuationTrace:
    """Terminal interactions along a path of barriers with state carry-over."""

    barriers: tuple[float, ...]
    terminal_interactions: tuple[float, ...]
    initial_interaction: float
    iterations_per_barrier: int


def _validate_map_parameters(
    feedback_strength: float,
    area: float,
    area_reference: float,
    barrier: float | None = None,
) -> None:
    for name, value in (
        ("feedback_strength", feedback_strength),
        ("area", area),
        ("area_reference", area_reference),
    ):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if barrier is not None and not isfinite(barrier):
        raise ValueError("barrier must be finite")


def _validate_interaction(value: float, name: str = "interaction") -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def _sigmoid(value: float) -> float:
    if value >= 0.0:
        inverse = exp(-value)
        return 1.0 / (1.0 + inverse)
    inverse = exp(value)
    return inverse / (1.0 + inverse)


def _logit(value: float) -> float:
    if not 0.0 < value < 1.0:
        raise ValueError("logit is defined here only for values strictly between zero and one")
    return log(value / (1.0 - value))


def canonical_gain(feedback_strength: float, area: float, area_reference: float = 1.0) -> float:
    """Return ``K = kappa * A / A_ref`` for the canonical map."""
    _validate_map_parameters(feedback_strength, area, area_reference)
    return feedback_strength * area / area_reference


def canonical_interaction_update(
    interaction: float,
    *,
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    barrier: float,
) -> float:
    """Apply one update of ``q_next = sigmoid(kappa * ((A/A_ref)q - theta))``."""
    q = _validate_interaction(interaction)
    _validate_map_parameters(feedback_strength, area, area_reference, barrier)
    return _sigmoid(feedback_strength * ((area / area_reference) * q - barrier))


def canonical_turning_points(
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    *,
    tolerance: float = 1e-12,
) -> tuple[float, float] | None:
    """Return the two turning points of the fixed-point equation when ``K > 4``.

    Fixed points satisfy ``logit(q) - Kq + kappa*theta = 0``.  Its derivative
    can vanish only for ``K > 4``; then the turning points are
    ``(1 +- sqrt(1 - 4/K)) / 2``.
    """
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    gain = canonical_gain(feedback_strength, area, area_reference)
    if gain <= 4.0 + tolerance:
        return None
    root = sqrt(1.0 - 4.0 / gain)
    return ((1.0 - root) / 2.0, (1.0 + root) / 2.0)


def canonical_bistable_barrier_interval(
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    *,
    tolerance: float = 1e-12,
) -> tuple[float, float] | None:
    """Return the open barrier interval with exactly three canonical fixed points.

    For ``K = kappa*A/A_ref > 4``, let ``q_- < q_+`` be the turning points.
    Strict bistability holds exactly for

        (K*q_- - logit(q_-))/kappa < theta < (K*q_+ - logit(q_+))/kappa.

    The endpoints are saddle-node boundaries and are not labelled bistable.
    """
    turning = canonical_turning_points(
        feedback_strength,
        area,
        area_reference,
        tolerance=tolerance,
    )
    if turning is None:
        return None
    q_low, q_high = turning
    gain = canonical_gain(feedback_strength, area, area_reference)
    lower = (gain * q_low - _logit(q_low)) / feedback_strength
    upper = (gain * q_high - _logit(q_high)) / feedback_strength
    return (lower, upper)


def _fixed_point_equation(q: float, gain: float, feedback_strength: float, barrier: float) -> float:
    return _logit(q) - gain * q + feedback_strength * barrier


def _bisect_root(function, left: float, right: float, *, tolerance: float) -> float:
    left_value = function(left)
    right_value = function(right)
    if abs(left_value) <= tolerance:
        return left
    if abs(right_value) <= tolerance:
        return right
    if left_value * right_value > 0.0:
        raise ValueError("root bracket must have opposite endpoint signs")
    for _ in range(256):
        middle = (left + right) / 2.0
        middle_value = function(middle)
        if abs(middle_value) <= tolerance or right - left <= tolerance:
            return middle
        if (left_value < 0.0) == (middle_value < 0.0):
            left, left_value = middle, middle_value
        else:
            right, right_value = middle, middle_value
    return (left + right) / 2.0


def _classify_stability(multiplier: float, tolerance: float) -> Stability:
    if abs(multiplier - 1.0) <= tolerance:
        return "neutral"
    return "stable" if multiplier < 1.0 else "unstable"


def _fixed_point(
    interaction: float,
    gain: float,
    *,
    tolerance: float,
) -> CanonicalFixedPoint:
    multiplier = gain * interaction * (1.0 - interaction)
    return CanonicalFixedPoint(
        interaction=interaction,
        local_multiplier=multiplier,
        stability=_classify_stability(multiplier, tolerance),
    )


def _fixed_points_for_regime(
    gain: float,
    feedback_strength: float,
    barrier: float,
    turning: tuple[float, float] | None,
    interval: tuple[float, float] | None,
    regime: BifurcationRegime,
    *,
    tolerance: float,
) -> tuple[CanonicalFixedPoint, ...]:
    edge_low = nextafter(0.0, 1.0)
    edge_high = nextafter(1.0, 0.0)
    equation = lambda q: _fixed_point_equation(q, gain, feedback_strength, barrier)

    if turning is None or interval is None:
        return (_fixed_point(_bisect_root(equation, edge_low, edge_high, tolerance=tolerance), gain, tolerance=tolerance),)

    q_low, q_high = turning
    if regime == "bistable":
        roots = (
            _bisect_root(equation, edge_low, q_low, tolerance=tolerance),
            _bisect_root(equation, q_low, q_high, tolerance=tolerance),
            _bisect_root(equation, q_high, edge_high, tolerance=tolerance),
        )
    elif regime == "lower_saddle_node":
        roots = (q_low, _bisect_root(equation, q_high, edge_high, tolerance=tolerance))
    elif regime == "upper_saddle_node":
        roots = (_bisect_root(equation, edge_low, q_low, tolerance=tolerance), q_high)
    elif barrier < interval[0]:
        roots = (_bisect_root(equation, q_high, edge_high, tolerance=tolerance),)
    else:
        roots = (_bisect_root(equation, edge_low, q_low, tolerance=tolerance),)
    return tuple(_fixed_point(root, gain, tolerance=tolerance) for root in roots)


def canonical_bifurcation_certificate(
    *,
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    barrier: float,
    tolerance: float = 1e-10,
) -> CanonicalBifurcationCertificate:
    """Certify the canonical fixed-point regime at one declared barrier.

    ``strict_bistability_certified`` means that the map has three fixed points:
    stable-low, unstable-middle, and stable-high.  It is an exact conclusion for
    this canonical map, up to the displayed numerical root tolerance.
    """
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive")
    _validate_map_parameters(feedback_strength, area, area_reference, barrier)
    gain = canonical_gain(feedback_strength, area, area_reference)
    turning = canonical_turning_points(
        feedback_strength,
        area,
        area_reference,
        tolerance=tolerance,
    )
    interval = canonical_bistable_barrier_interval(
        feedback_strength,
        area,
        area_reference,
        tolerance=tolerance,
    )

    if interval is None:
        regime: BifurcationRegime = "monostable"
    else:
        lower, upper = interval
        if abs(barrier - lower) <= tolerance:
            regime = "lower_saddle_node"
        elif abs(barrier - upper) <= tolerance:
            regime = "upper_saddle_node"
        elif lower < barrier < upper:
            regime = "bistable"
        else:
            regime = "monostable"

    fixed_points = _fixed_points_for_regime(
        gain,
        feedback_strength,
        barrier,
        turning,
        interval,
        regime,
        tolerance=tolerance,
    )
    return CanonicalBifurcationCertificate(
        feedback_strength=feedback_strength,
        area=area,
        area_reference=area_reference,
        barrier=barrier,
        gain=gain,
        turning_points=turning,
        bistable_barrier_interval=interval,
        regime=regime,
        fixed_points=fixed_points,
        strict_bistability_certified=regime == "bistable",
    )


def canonical_h1_certificate(
    *,
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    barrier: float,
    trait_parameters: DynamicsParameters,
    tolerance: float = 1e-10,
) -> CanonicalH1Certificate:
    """Combine canonical bistability with the declared high-trait viability margin.

    The result certifies the H1 mechanism for this specified canonical map only
    when two stable branches exist and the high-trait margin is negative on the
    low branch and positive on the high branch.
    """
    bifurcation = canonical_bifurcation_certificate(
        feedback_strength=feedback_strength,
        area=area,
        area_reference=area_reference,
        barrier=barrier,
        tolerance=tolerance,
    )
    stable = tuple(point for point in bifurcation.fixed_points if point.stability == "stable")
    if len(stable) != 2:
        return CanonicalH1Certificate(bifurcation, None, None, False)

    low, high = stable
    low_summary = trait_space_summary(low.interaction, trait_parameters)
    high_summary = trait_space_summary(high.interaction, trait_parameters)
    low_branch = HighTraitBranch(
        interaction=low.interaction,
        high_trait_margin=low_summary.high_trait_margin,
        high_trait_component_present=low_summary.high_trait_component_present,
    )
    high_branch = HighTraitBranch(
        interaction=high.interaction,
        high_trait_margin=high_summary.high_trait_margin,
        high_trait_component_present=high_summary.high_trait_component_present,
    )
    return CanonicalH1Certificate(
        bifurcation=bifurcation,
        low_stable_branch=low_branch,
        high_stable_branch=high_branch,
        branch_dependent_high_trait_mode=(
            low_branch.high_trait_margin < 0.0 < high_branch.high_trait_margin
        ),
    )


def iterate_canonical_map(
    initial_interaction: float,
    *,
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    barrier: float,
    iterations: int = 200,
) -> CanonicalOrbit:
    """Iterate the canonical map from one initial state."""
    if iterations < 0:
        raise ValueError("iterations must be non-negative")
    current = _validate_interaction(initial_interaction, "initial_interaction")
    _validate_map_parameters(feedback_strength, area, area_reference, barrier)
    values = [current]
    for _ in range(iterations):
        current = canonical_interaction_update(
            current,
            feedback_strength=feedback_strength,
            area=area,
            area_reference=area_reference,
            barrier=barrier,
        )
        values.append(current)
    return CanonicalOrbit(initial_interaction=values[0], values=tuple(values))


def follow_barrier_path(
    barriers: Sequence[float],
    *,
    initial_interaction: float,
    feedback_strength: float,
    area: float,
    area_reference: float = 1.0,
    iterations_per_barrier: int = 200,
) -> BarrierContinuationTrace:
    """Continue the map along a barrier path, retaining the prior terminal state.

    An increasing and a decreasing path through the same bistable interval can
    end on different branches.  That path dependence is the numerical hysteresis
    check associated with the analytic saddle-node certificate.
    """
    if not barriers:
        raise ValueError("barriers must be nonempty")
    if iterations_per_barrier < 1:
        raise ValueError("iterations_per_barrier must be at least one")
    _validate_map_parameters(feedback_strength, area, area_reference)
    current = _validate_interaction(initial_interaction, "initial_interaction")
    terminal: list[float] = []
    declared_barriers = tuple(float(barrier) for barrier in barriers)
    for barrier in declared_barriers:
        _validate_map_parameters(feedback_strength, area, area_reference, barrier)
        orbit = iterate_canonical_map(
            current,
            feedback_strength=feedback_strength,
            area=area,
            area_reference=area_reference,
            barrier=barrier,
            iterations=iterations_per_barrier,
        )
        current = orbit.terminal_interaction
        terminal.append(current)
    return BarrierContinuationTrace(
        barriers=declared_barriers,
        terminal_interactions=tuple(terminal),
        initial_interaction=_validate_interaction(initial_interaction, "initial_interaction"),
        iterations_per_barrier=iterations_per_barrier,
    )
