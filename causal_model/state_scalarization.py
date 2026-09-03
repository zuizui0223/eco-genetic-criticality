"""Finite criterion for when multiple ecological targets admit one monotone scalar state."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Hashable, Mapping, Sequence


@dataclass(frozen=True)
class ScalarizationAudit:
    target_names: tuple[str, ...]
    chain_under_product_order: bool
    scalar_state: dict[Hashable, int] | None
    crossing_pair: tuple[Hashable, Hashable] | None


def _vector(values: Mapping[str, float], target_names: Sequence[str]) -> tuple[float, ...]:
    return tuple(float(values[name]) for name in target_names)


def weakly_leq(left: Sequence[float], right: Sequence[float]) -> bool:
    return all(a <= b for a, b in zip(left, right))


def crossing(left: Sequence[float], right: Sequence[float]) -> bool:
    """Return whether each vector is better on at least one coordinate."""
    left_better = any(a > b for a, b in zip(left, right))
    right_better = any(a < b for a, b in zip(left, right))
    return left_better and right_better


def common_monotone_scalar_audit(
    *,
    target_values: Mapping[Hashable, Mapping[str, float]],
    target_names: Sequence[str],
) -> ScalarizationAudit:
    """Audit existence of a sufficient scalar with a common 'higher is no worse' order.

    A scalar state h is sought such that every target T_j factors through h and
    every reconstructed target is nondecreasing in h.  For a finite system this
    exists iff the distinct target vectors form a chain under coordinatewise
    product order.
    """
    names = tuple(target_names)
    if not names:
        raise ValueError("at least one target is required")
    if not target_values:
        raise ValueError("at least one state is required")

    vectors = {state: _vector(values, names) for state, values in target_values.items()}
    states = tuple(vectors)
    for i, left_state in enumerate(states):
        for right_state in states[i + 1 :]:
            left, right = vectors[left_state], vectors[right_state]
            if not (weakly_leq(left, right) or weakly_leq(right, left)):
                return ScalarizationAudit(names, False, None, (left_state, right_state))

    unique_vectors = sorted(set(vectors.values()), key=lambda v: (sum(v), v))
    # Under a product-order chain, sum is strictly increasing between distinct
    # comparable vectors and therefore gives a valid chain order.
    rank = {vector: i for i, vector in enumerate(unique_vectors)}
    scalar = {state: rank[vector] for state, vector in vectors.items()}
    return ScalarizationAudit(names, True, scalar, None)


def scalarization_is_valid(
    *,
    target_values: Mapping[Hashable, Mapping[str, float]],
    target_names: Sequence[str],
    scalar_state: Mapping[Hashable, int],
) -> bool:
    """Check factorization plus common monotonicity for one proposed scalar state."""
    names = tuple(target_names)
    states = tuple(target_values)
    if set(states) != set(scalar_state):
        return False

    # Equal scalar labels must imply equal target vectors (factorization).
    for i, a in enumerate(states):
        va = _vector(target_values[a], names)
        for b in states[i + 1 :]:
            vb = _vector(target_values[b], names)
            ha, hb = scalar_state[a], scalar_state[b]
            if ha == hb and va != vb:
                return False
            if ha < hb and not weakly_leq(va, vb):
                return False
            if hb < ha and not weakly_leq(vb, va):
                return False
    return True
