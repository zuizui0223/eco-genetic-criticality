"""General mathematical principles for eco-genetic criticality.

This module deliberately separates theorem-level assumptions from later special
models such as logistic interaction feedback or Wright--Fisher transmission.

The theorem layer has five parts.

G0 -- Finite transmission variance erodes expected heterozygosity.
    For a post-selection allele frequency p* and next-generation random frequency
    P', if E[P'|p*]=p* and Var(P'|p*)=v>0, then

        E[H(P')|p*] = H(p*) - 2v < H(p*),

    where H(p)=2p(1-p). No binomial transmission assumption is required.

P0 -- Global feedback contraction bound.
    For q = g(kappa(Aq-theta)), if sup_x |g'(x)| <= M and kappa*A*M < 1,
    the update is a contraction and therefore has a unique fixed point. The
    reverse inequality only removes this universal uniqueness guarantee; it does
    not itself prove bistability without further shape assumptions.

P1 -- Trait-mode lifting.
    If two stable interaction states q_L<q_H exist and a high-trait mode has
    viability margin negative at q_L and positive at q_H, then its existence is
    branch/history dependent.

P2 -- Partition non-additivity.
    If a collective interaction mechanism requires a patchwise threshold A>A_c,
    total area alone is insufficient: a partition with max_j A_j<=A_c cannot
    express that mechanism even when sum_j A_j>A_c.

G1 -- Conditional eco-genetic coupling.
    If effective reproductive size Psi(A,q,xi) increases with q and transmission
    variance decreases with effective size, then lower interaction branches have
    greater expected heterozygosity erosion. This is conditional on those
    biological monotonicities; they are not universal facts.

The proofs are in ``docs/general_eco_genetic_principles.md``. The functions here
are small algebraic helpers and regression checks, not replacements for proofs.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isclose
from typing import Sequence


@dataclass(frozen=True)
class TransmissionMomentResult:
    """Exact heterozygosity identity from the first two transmission moments."""

    post_selection_frequency: float
    next_frequency_mean: float
    next_frequency_variance: float
    heterozygosity_after_selection: float
    expected_heterozygosity_after_transmission: float
    expected_heterozygosity_loss: float
    unbiased: bool
    finite_variance: bool


@dataclass(frozen=True)
class FeedbackBound:
    """Global contraction certificate for an interaction fixed-point update."""

    patch_size: float
    feedback_strength: float
    max_response_slope: float
    global_lipschitz_bound: float
    uniqueness_certified: bool
    bistability_not_ruled_out: bool


@dataclass(frozen=True)
class TraitModeBranchResult:
    """Branch-specific high-trait viability margins."""

    low_interaction_state: float
    high_interaction_state: float
    low_branch_margin: float
    high_branch_margin: float
    branch_dependent_mode: bool


@dataclass(frozen=True)
class PartitionThresholdResult:
    """Patchwise-capacity conclusion under an explicitly supplied threshold."""

    total_area: float
    patch_sizes: tuple[float, ...]
    critical_patch_size: float
    total_area_exceeds_threshold: bool
    any_patch_exceeds_threshold: bool
    collective_mechanism_possible_by_threshold: bool


@dataclass(frozen=True)
class EcoGeneticOrdering:
    """Conditional ordering of transmission variance and diversity erosion across branches."""

    low_interaction_state: float
    high_interaction_state: float
    low_effective_size: float
    high_effective_size: float
    low_transmission_variance: float
    high_transmission_variance: float
    low_expected_heterozygosity_loss: float
    high_expected_heterozygosity_loss: float
    ordering_supported: bool


def heterozygosity(allele_frequency: float) -> float:
    """Return gene diversity H(p)=2p(1-p) for a biallelic locus."""
    if not 0.0 <= allele_frequency <= 1.0:
        raise ValueError("allele_frequency must lie in [0, 1]")
    return 2.0 * allele_frequency * (1.0 - allele_frequency)


def transmission_moment_identity(
    post_selection_frequency: float,
    next_frequency_mean: float,
    next_frequency_variance: float,
    *,
    tolerance: float = 1e-12,
) -> TransmissionMomentResult:
    """Apply E[H(P')]=H(E[P'])-2 Var(P') exactly.

    The identity follows from H(x)=2x-2x^2. ``unbiased`` records whether the
    provided transmission mean equals p* within tolerance; the algebraic identity
    itself does not need that equality. A strictly positive variance only implies
    erosion relative to H(p*) when transmission is unbiased.
    """
    for name, value in (
        ("post_selection_frequency", post_selection_frequency),
        ("next_frequency_mean", next_frequency_mean),
    ):
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1]")
    if next_frequency_variance < 0.0:
        raise ValueError("next_frequency_variance must be non-negative")
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    h_after_selection = heterozygosity(post_selection_frequency)
    expected_h = heterozygosity(next_frequency_mean) - 2.0 * next_frequency_variance
    if expected_h < -tolerance:
        raise ValueError("moments imply impossible negative expected heterozygosity")
    expected_h = max(0.0, expected_h)
    loss = h_after_selection - expected_h
    unbiased = isclose(next_frequency_mean, post_selection_frequency, abs_tol=tolerance)
    return TransmissionMomentResult(
        post_selection_frequency=post_selection_frequency,
        next_frequency_mean=next_frequency_mean,
        next_frequency_variance=next_frequency_variance,
        heterozygosity_after_selection=h_after_selection,
        expected_heterozygosity_after_transmission=expected_h,
        expected_heterozygosity_loss=loss,
        unbiased=unbiased,
        finite_variance=next_frequency_variance > tolerance,
    )


def finite_transmission_strictly_erodes(result: TransmissionMomentResult, *, tolerance: float = 1e-12) -> bool:
    """Return G0's strict conclusion under unbiased positive-variance transmission."""
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    return (
        result.unbiased
        and result.finite_variance
        and result.expected_heterozygosity_after_transmission
        < result.heterozygosity_after_selection - tolerance
    )


def feedback_contraction_bound(
    patch_size: float,
    feedback_strength: float,
    max_response_slope: float,
) -> FeedbackBound:
    """Return the P0 global contraction certificate.

    For q -> g(kappa(Aq-theta)), a global derivative bound M for g gives a
    Lipschitz bound kappa*A*M. Strictly below one, Banach's theorem certifies a
    unique fixed point for every barrier theta. At or above one, uniqueness is not
    certified by this bound; multistability remains model-shape dependent.
    """
    if patch_size <= 0.0:
        raise ValueError("patch_size must be positive")
    if feedback_strength <= 0.0:
        raise ValueError("feedback_strength must be positive")
    if max_response_slope < 0.0:
        raise ValueError("max_response_slope must be non-negative")
    lipschitz = patch_size * feedback_strength * max_response_slope
    return FeedbackBound(
        patch_size=patch_size,
        feedback_strength=feedback_strength,
        max_response_slope=max_response_slope,
        global_lipschitz_bound=lipschitz,
        uniqueness_certified=lipschitz < 1.0,
        bistability_not_ruled_out=lipschitz >= 1.0,
    )


def trait_mode_branch_result(
    low_interaction_state: float,
    high_interaction_state: float,
    low_branch_margin: float,
    high_branch_margin: float,
) -> TraitModeBranchResult:
    """Apply P1 from branch-specific viability margins.

    A margin is max_{z in Z_H} W(z;q)-tau for a declared high-trait region. The
    theorem requires negative margin on the low branch and positive margin on the
    high branch; exact zero is a boundary and is deliberately not called a robust
    branch difference.
    """
    if not 0.0 <= low_interaction_state < high_interaction_state <= 1.0:
        raise ValueError("interaction states must satisfy 0 <= q_low < q_high <= 1")
    return TraitModeBranchResult(
        low_interaction_state=low_interaction_state,
        high_interaction_state=high_interaction_state,
        low_branch_margin=low_branch_margin,
        high_branch_margin=high_branch_margin,
        branch_dependent_mode=low_branch_margin < 0.0 < high_branch_margin,
    )


def partition_threshold_result(
    patch_sizes: Sequence[float],
    critical_patch_size: float,
) -> PartitionThresholdResult:
    """Apply P2 under a supplied patchwise collective-capacity threshold."""
    sizes = tuple(float(size) for size in patch_sizes)
    if not sizes or any(size <= 0.0 for size in sizes):
        raise ValueError("patch_sizes must be nonempty and positive")
    if critical_patch_size <= 0.0:
        raise ValueError("critical_patch_size must be positive")
    total = sum(sizes)
    any_exceeds = any(size > critical_patch_size for size in sizes)
    return PartitionThresholdResult(
        total_area=total,
        patch_sizes=sizes,
        critical_patch_size=critical_patch_size,
        total_area_exceeds_threshold=total > critical_patch_size,
        any_patch_exceeds_threshold=any_exceeds,
        collective_mechanism_possible_by_threshold=any_exceeds,
    )


def eco_genetic_ordering(
    low_interaction_state: float,
    high_interaction_state: float,
    low_effective_size: float,
    high_effective_size: float,
    low_transmission_variance: float,
    high_transmission_variance: float,
) -> EcoGeneticOrdering:
    """Record G1's conditional branch ordering.

    The caller must derive effective sizes and transmission variances from a
    declared life cycle. This helper refuses to infer them from patch area alone.
    Under the supplied ordering, expected heterozygosity loss is exactly twice the
    transmission variance by G0 for unbiased transmission.
    """
    if not 0.0 <= low_interaction_state < high_interaction_state <= 1.0:
        raise ValueError("interaction states must satisfy 0 <= q_low < q_high <= 1")
    if low_effective_size <= 0.0 or high_effective_size <= 0.0:
        raise ValueError("effective sizes must be positive")
    if low_transmission_variance < 0.0 or high_transmission_variance < 0.0:
        raise ValueError("transmission variances must be non-negative")
    low_loss = 2.0 * low_transmission_variance
    high_loss = 2.0 * high_transmission_variance
    supported = (
        low_effective_size < high_effective_size
        and low_transmission_variance > high_transmission_variance
        and low_loss > high_loss
    )
    return EcoGeneticOrdering(
        low_interaction_state=low_interaction_state,
        high_interaction_state=high_interaction_state,
        low_effective_size=low_effective_size,
        high_effective_size=high_effective_size,
        low_transmission_variance=low_transmission_variance,
        high_transmission_variance=high_transmission_variance,
        low_expected_heterozygosity_loss=low_loss,
        high_expected_heterozygosity_loss=high_loss,
        ordering_supported=supported,
    )
