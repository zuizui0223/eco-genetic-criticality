"""A transparent sufficient-condition chain linking canonical H1, H2 and H3.

This module does not convert separate assumptions into a universal ecological
law.  Instead it records precisely when four declared components can be
composed:

1. H1 has two stable interaction branches with a high-trait margin switch;
2. effective size is an increasing affine function of interaction;
3. the H2 expected-drift map has a strict genetic lead on the low branch; and
4. H3 fragmentation removes local support, unless external arrivals satisfy a
   separately declared rescue condition.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from causal_model.canonical_h1_bifurcation import CanonicalH1Certificate
from causal_model.canonical_h2_h3_certificates import (
    CanonicalH2LeadCertificate,
    CanonicalH3FragmentationCertificate,
    CanonicalH3RescueCertificate,
)


@dataclass(frozen=True)
class CoupledEcoGeneticCertificate:
    """Composed theorem status for the explicitly declared H1→H2→H3 skeleton."""

    h1_branch_switch_certified: bool
    h3_fragmentation_removes_support_certified: bool
    low_branch_effective_population_size: float | None
    high_branch_effective_population_size: float | None
    lower_effective_size_on_low_branch_certified: bool
    h2_expected_genetic_lead_certified: bool
    external_demographic_rescue_certified: bool
    external_high_trait_rescue_certified: bool
    fragmentation_to_genetic_warning_chain_certified: bool
    rescue_interrupts_fragmentation_chain: bool


def affine_effective_population_size(
    interaction: float,
    *,
    baseline_effective_population_size: float,
    interaction_slope: float,
) -> float:
    """Evaluate the declared closure ``N_e(q)=N_0 + beta*q``.

    A positive slope is required in the coupled theorem: it is the explicit
    ecological assumption that stronger interaction raises effective size.
    """
    for name, value in (
        ("interaction", interaction),
        ("baseline_effective_population_size", baseline_effective_population_size),
        ("interaction_slope", interaction_slope),
    ):
        if not isfinite(value):
            raise ValueError(f"{name} must be finite")
    if not 0.0 <= interaction <= 1.0:
        raise ValueError("interaction must lie in [0, 1]")
    if baseline_effective_population_size <= 0.0:
        raise ValueError("baseline_effective_population_size must be positive")
    if interaction_slope < 0.0:
        raise ValueError("interaction_slope must be non-negative")
    return baseline_effective_population_size + interaction_slope * interaction


def coupled_eco_genetic_certificate(
    *,
    h1: CanonicalH1Certificate,
    h2_low_branch: CanonicalH2LeadCertificate,
    h3_fragmentation: CanonicalH3FragmentationCertificate,
    baseline_effective_population_size: float,
    interaction_slope: float,
    h3_rescue: CanonicalH3RescueCertificate | None = None,
    tolerance: float = 1e-12,
) -> CoupledEcoGeneticCertificate:
    """Compose canonical H1, H2 and H3 conclusions under explicit closures.

    The chain is certified only if H1 has a branch-dependent high-trait mode,
    equal fragmentation removes local support, the low H1 branch has lower
    declared effective size than the high branch, and the supplied H2 map has a
    strict expected genetic lead.  Rescue does not invalidate the theorem; it
    identifies the stated mechanism by which the fragmentation chain is
    interrupted in a recipient patch.
    """
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    h1_ok = h1.branch_dependent_high_trait_mode
    low_ne: float | None = None
    high_ne: float | None = None
    lower_ne = False
    if h1_ok:
        assert h1.low_stable_branch is not None
        assert h1.high_stable_branch is not None
        low_ne = affine_effective_population_size(
            h1.low_stable_branch.interaction,
            baseline_effective_population_size=baseline_effective_population_size,
            interaction_slope=interaction_slope,
        )
        high_ne = affine_effective_population_size(
            h1.high_stable_branch.interaction,
            baseline_effective_population_size=baseline_effective_population_size,
            interaction_slope=interaction_slope,
        )
        lower_ne = low_ne + tolerance < high_ne

    demographic_rescue = bool(h3_rescue and h3_rescue.rescue_certified)
    high_trait_rescue = bool(h3_rescue and h3_rescue.high_trait_rescue_certified)
    chain = (
        h1_ok
        and h3_fragmentation.fragmentation_removes_local_support_certified
        and lower_ne
        and h2_low_branch.strict_expected_genetic_lead_certified
    )
    return CoupledEcoGeneticCertificate(
        h1_branch_switch_certified=h1_ok,
        h3_fragmentation_removes_support_certified=h3_fragmentation.fragmentation_removes_local_support_certified,
        low_branch_effective_population_size=low_ne,
        high_branch_effective_population_size=high_ne,
        lower_effective_size_on_low_branch_certified=lower_ne,
        h2_expected_genetic_lead_certified=h2_low_branch.strict_expected_genetic_lead_certified,
        external_demographic_rescue_certified=demographic_rescue,
        external_high_trait_rescue_certified=high_trait_rescue,
        fragmentation_to_genetic_warning_chain_certified=chain,
        rescue_interrupts_fragmentation_chain=chain and (demographic_rescue or high_trait_rescue),
    )
