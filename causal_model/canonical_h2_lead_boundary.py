"""Closed-form first-passage boundary for the canonical deterministic H2 map."""
from __future__ import annotations

from dataclasses import dataclass
from math import ceil, isfinite, log


@dataclass(frozen=True)
class CanonicalH2LeadBoundary:
    """Analytic threshold-crossing times and the exact discrete lead relation."""

    effective_population_size: float
    trait_retention: float
    heterozygosity_retention: float
    continuous_genetic_crossing_time: float
    continuous_trait_crossing_time: float | None
    discrete_genetic_warning_time: int
    discrete_trait_loss_time: int | None
    strict_expected_genetic_lead_certified: bool


def _probability(value: float, name: str, *, open_lower: bool = False) -> float:
    value = float(value)
    if not isfinite(value) or not 0.0 <= value <= 1.0 or (open_lower and value == 0.0):
        interval = "(0, 1]" if open_lower else "[0, 1]"
        raise ValueError(f"{name} must lie in {interval}")
    return value


def canonical_h2_lead_boundary(
    *,
    initial_heterozygosity: float,
    heterozygosity_threshold: float,
    initial_trait_abundance: float,
    trait_threshold: float,
    effective_population_size: float,
    trait_retention: float,
) -> CanonicalH2LeadBoundary:
    """Solve canonical H2 first passages without iterating generations.

    For ``H[t]=H0*(1-1/(2N_e))**t`` and ``T[t]=T0*r_T**t``, the first integer
    crossings are ceilings of their log-ratio crossing times.  If ``r_T=1``,
    the trait event is absent and remains represented as ``None``.
    """
    h0 = _probability(initial_heterozygosity, "initial_heterozygosity", open_lower=True)
    h_star = _probability(heterozygosity_threshold, "heterozygosity_threshold")
    if h_star >= h0:
        raise ValueError("heterozygosity_threshold must be below initial_heterozygosity")
    if not isfinite(initial_trait_abundance) or initial_trait_abundance <= 0.0:
        raise ValueError("initial_trait_abundance must be finite and positive")
    if not isfinite(trait_threshold) or not 0.0 < trait_threshold < initial_trait_abundance:
        raise ValueError("trait_threshold must be positive and below initial_trait_abundance")
    if not isfinite(effective_population_size) or effective_population_size < 1.0:
        raise ValueError("effective_population_size must be finite and at least one")
    r_t = _probability(trait_retention, "trait_retention", open_lower=True)
    alpha = 1.0 - 1.0 / (2.0 * effective_population_size)
    genetic_continuous = log(h_star / h0) / log(alpha)
    genetic_discrete = max(0, ceil(genetic_continuous))
    if r_t == 1.0:
        trait_continuous = None
        trait_discrete = None
    else:
        trait_continuous = log(trait_threshold / initial_trait_abundance) / log(r_t)
        trait_discrete = max(0, ceil(trait_continuous))
    return CanonicalH2LeadBoundary(
        effective_population_size=float(effective_population_size),
        trait_retention=r_t,
        heterozygosity_retention=alpha,
        continuous_genetic_crossing_time=genetic_continuous,
        continuous_trait_crossing_time=trait_continuous,
        discrete_genetic_warning_time=genetic_discrete,
        discrete_trait_loss_time=trait_discrete,
        strict_expected_genetic_lead_certified=(
            trait_discrete is not None and genetic_discrete < trait_discrete
        ),
    )
