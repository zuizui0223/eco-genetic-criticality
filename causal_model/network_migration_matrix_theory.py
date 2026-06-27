"""Exact migration-matrix bounds for the H3 gene-flow component.

The migration update is written with an explicit destination-by-source matrix

    p_next[i] = sum_j M[i, j] p[j],

where every row of ``M`` is non-negative and sums to one.  This covers
complete-graph mixing, asymmetric source-sink networks, stepping-stone chains,
and distance-derived kernels.  The results below concern allele-frequency
mixing only; they do not claim demographic rescue, trait-bin recolonisation, or
patch persistence without additional life-cycle equations.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Sequence


@dataclass(frozen=True)
class MigrationMatrixCertificate:
    """Validation and common-floor result for one declared migration matrix."""

    matrix: tuple[tuple[float, ...], ...]
    source_floor: float
    destination_lower_bounds: tuple[float, ...]
    common_floor_preserved: bool


@dataclass(frozen=True)
class FocalRescueCertificate:
    """A sufficient lower-bound certificate for one focal destination patch."""

    destination_index: int
    source_lower_bounds: tuple[float, ...]
    destination_lower_bound: float
    target_floor: float
    rescue_certified: bool


def _probability(value: float, name: str) -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return value


def validate_destination_by_source_matrix(
    matrix: Sequence[Sequence[float]],
    *,
    tolerance: float = 1e-12,
) -> tuple[tuple[float, ...], ...]:
    """Validate a finite row-stochastic destination-by-source matrix.

    ``matrix[i][j]`` is the contribution weight from source ``j`` to destination
    ``i``.  Each destination row must be a convex combination of source values.
    """
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    rows = tuple(tuple(float(value) for value in row) for row in matrix)
    if not rows:
        raise ValueError("migration matrix must be nonempty")
    width = len(rows)
    if any(len(row) != width for row in rows):
        raise ValueError("migration matrix must be square")
    for row in rows:
        if any(not isfinite(value) or value < -tolerance for value in row):
            raise ValueError("migration weights must be finite and non-negative")
        if abs(sum(row) - 1.0) > tolerance:
            raise ValueError("each destination row must sum to one")
    return tuple(tuple(max(0.0, value) for value in row) for row in rows)


def mix_allele_frequencies(
    source_frequencies: Sequence[float],
    matrix: Sequence[Sequence[float]],
) -> tuple[float, ...]:
    """Apply one exact convex-combination allele migration update."""
    checked_matrix = validate_destination_by_source_matrix(matrix)
    sources = tuple(_probability(value, "source frequency") for value in source_frequencies)
    if len(sources) != len(checked_matrix):
        raise ValueError("source_frequencies length must match migration matrix size")
    return tuple(
        sum(weight * frequency for weight, frequency in zip(row, sources))
        for row in checked_matrix
    )


def complete_graph_mixing_matrix(
    source_weights: Sequence[float],
    migration_rate: float,
) -> tuple[tuple[float, ...], ...]:
    """Return the matrix for the repository's current global-mean mixing update.

    With normalised source weights ``w_j`` and migration rate ``m``, the current
    update is ``(1-m)p_i + m*sum_j w_j p_j``.  Its explicit matrix is
    ``M[i,j]=(1-m)I[i=j]+m*w_j``.
    """
    weights = tuple(float(value) for value in source_weights)
    if not weights or any(not isfinite(value) or value < 0.0 for value in weights):
        raise ValueError("source_weights must be nonempty, finite, and non-negative")
    total = sum(weights)
    if total <= 0.0:
        raise ValueError("source_weights must have positive total")
    migration_rate = _probability(migration_rate, "migration_rate")
    normalised = tuple(value / total for value in weights)
    return tuple(
        tuple(
            (1.0 - migration_rate if destination == source else 0.0)
            + migration_rate * normalised[source]
            for source in range(len(normalised))
        )
        for destination in range(len(normalised))
    )


def common_floor_certificate(
    matrix: Sequence[Sequence[float]],
    source_floor: float,
) -> MigrationMatrixCertificate:
    """Certify preservation of a common allele floor under any valid network.

    If every source frequency is at least ``source_floor``, every destination is
    also at least that floor because each destination is a convex combination of
    source values.  This generalises the complete-graph global-mean result.
    """
    checked_matrix = validate_destination_by_source_matrix(matrix)
    floor = _probability(source_floor, "source_floor")
    lower_bounds = tuple(sum(weight * floor for weight in row) for row in checked_matrix)
    return MigrationMatrixCertificate(
        matrix=checked_matrix,
        source_floor=floor,
        destination_lower_bounds=lower_bounds,
        common_floor_preserved=all(bound >= floor - 1e-12 for bound in lower_bounds),
    )


def destination_lower_bound(
    matrix: Sequence[Sequence[float]],
    source_lower_bounds: Sequence[float],
    destination_index: int,
) -> float:
    """Return the sharp convex-combination lower bound for one destination."""
    checked_matrix = validate_destination_by_source_matrix(matrix)
    bounds = tuple(_probability(value, "source lower bound") for value in source_lower_bounds)
    if len(bounds) != len(checked_matrix):
        raise ValueError("source_lower_bounds length must match migration matrix size")
    if not 0 <= destination_index < len(checked_matrix):
        raise ValueError("destination_index is out of range")
    return sum(weight * bound for weight, bound in zip(checked_matrix[destination_index], bounds))


def focal_rescue_certificate(
    matrix: Sequence[Sequence[float]],
    source_lower_bounds: Sequence[float],
    destination_index: int,
    target_floor: float,
) -> FocalRescueCertificate:
    """Certify migration rescue when the destination lower bound reaches a target.

    This is a sufficient condition.  It requires declared lower bounds for all
    source patches; it does not assume that a source-sink network automatically
    rescues a patch merely because migration is nonzero.
    """
    bounds = tuple(_probability(value, "source lower bound") for value in source_lower_bounds)
    target = _probability(target_floor, "target_floor")
    lower = destination_lower_bound(matrix, bounds, destination_index)
    return FocalRescueCertificate(
        destination_index=destination_index,
        source_lower_bounds=bounds,
        destination_lower_bound=lower,
        target_floor=target,
        rescue_certified=lower >= target,
    )
