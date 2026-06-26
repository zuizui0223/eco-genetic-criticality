"""High-probability H-alpha lead bounds for finite trait recruitment.

This module bridges the deterministic L3 certificate to a finite-bin stochastic
closure. A nonzero realised trait abundance cannot be guaranteed pathwise under
multinomial recruitment: any positive-probability trait bin can be missed in a
finite draw. Therefore this module proves lower bounds on the probability of a
lead event rather than claiming an all-random-realizations guarantee.

At a declared time t, suppose:

    E[H_t] <= h_0 lambda_bar^t

and realised high-trait recruitment in each generation s=1,...,t has a
conditional Binomial lower-tail failure bound. Markov's inequality bounds the
chance that H_t is still above its warning threshold, and a union bound controls
trait-loss risk through t. Their intersection gives a conservative lower bound
on:

    P(tau_H <= t < tau_trait_realised).

The functions do not infer these assumptions from a simulator. The caller must
establish the diversity and recruitment bounds for the declared closure.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import exp


@dataclass(frozen=True)
class BinomialTraitPersistenceBound:
    """One-generation lower-tail bound for finite high-trait recruitment."""

    cohort_size_lower_bound: int
    high_trait_recruit_probability_lower_bound: float
    occupancy_threshold: int
    expected_high_trait_abundance_lower_bound: float
    per_generation_failure_upper_bound: float


@dataclass(frozen=True)
class ProbabilisticLeadCertificate:
    """Conservative lower bound on a genetic lead event at a declared time."""

    time: int
    expected_diversity_upper_bound: float
    diversity_warning_failure_upper_bound: float
    trait_persistence_failure_upper_bound: float
    lead_probability_lower_bound: float


def binomial_trait_persistence_bound(
    cohort_size_lower_bound: int,
    high_trait_recruit_probability_lower_bound: float,
    occupancy_threshold: int,
) -> BinomialTraitPersistenceBound:
    """Bound one-generation realised trait-loss risk with Chernoff's inequality.

    Assume conditional on every allowed history that the next realised high-trait
    abundance dominates a Binomial(n_min, p_min) variable, where ``n_min`` is the
    cohort-size lower bound and ``p_min`` the high-trait recruitment-probability
    lower bound. For occupancy threshold a < mu=n_min*p_min:

        P(N_H <= a) <= exp[-mu (1-a/mu)^2 / 2].

    If a >= mu, this function returns the trivial upper bound one rather than
    pretending that a lower-tail certificate exists.
    """
    if cohort_size_lower_bound < 1:
        raise ValueError("cohort_size_lower_bound must be at least one")
    if not 0.0 < high_trait_recruit_probability_lower_bound <= 1.0:
        raise ValueError("high_trait_recruit_probability_lower_bound must lie in (0, 1]")
    if occupancy_threshold < 0:
        raise ValueError("occupancy_threshold must be non-negative")

    mu = cohort_size_lower_bound * high_trait_recruit_probability_lower_bound
    if occupancy_threshold >= mu:
        failure = 1.0
    else:
        delta = 1.0 - occupancy_threshold / mu
        failure = exp(-mu * delta * delta / 2.0)
    return BinomialTraitPersistenceBound(
        cohort_size_lower_bound=cohort_size_lower_bound,
        high_trait_recruit_probability_lower_bound=high_trait_recruit_probability_lower_bound,
        occupancy_threshold=occupancy_threshold,
        expected_high_trait_abundance_lower_bound=mu,
        per_generation_failure_upper_bound=min(1.0, failure),
    )


def trait_persistence_union_bound(
    per_generation_failure_upper_bound: float,
    time: int,
) -> float:
    """Upper-bound any realised trait-loss event through generations 1,...,time."""
    if not 0.0 <= per_generation_failure_upper_bound <= 1.0:
        raise ValueError("per_generation_failure_upper_bound must lie in [0, 1]")
    if time < 0:
        raise ValueError("time must be non-negative")
    return min(1.0, time * per_generation_failure_upper_bound)


def probabilistic_h_alpha_lead_certificate(
    initial_diversity: float,
    warning_threshold: float,
    multiplier_upper_bound: float,
    time: int,
    trait_persistence_failure_upper_bound: float,
) -> ProbabilisticLeadCertificate:
    """Certify a conservative lower probability for tau_H <= t < tau_trait.

    Assumptions:

        E[H_t] <= h_0 lambda_bar^t,  0 < lambda_bar < 1,

    and ``trait_persistence_failure_upper_bound`` upper-bounds the probability
    that realised high-trait occupancy has been lost at or before time ``t``.

    Markov's inequality gives:

        P(H_t > h_warn) <= E[H_t] / h_warn.

    Since H_t <= h_warn implies tau_H <= t, and combining the two failure events
    with a union bound yields:

        P(tau_H <= t < tau_trait) >=
        max(0, 1 - h_0 lambda_bar^t / h_warn - epsilon_trait).

    The result is a probability certificate, not a pathwise guarantee.
    """
    if not 0.0 < initial_diversity <= 1.0:
        raise ValueError("initial_diversity must lie in (0, 1]")
    if not 0.0 < warning_threshold <= 1.0:
        raise ValueError("warning_threshold must lie in (0, 1]")
    if not 0.0 < multiplier_upper_bound < 1.0:
        raise ValueError("multiplier_upper_bound must lie in (0, 1)")
    if time < 0:
        raise ValueError("time must be non-negative")
    if not 0.0 <= trait_persistence_failure_upper_bound <= 1.0:
        raise ValueError("trait_persistence_failure_upper_bound must lie in [0, 1]")

    expected_diversity = initial_diversity * multiplier_upper_bound**time
    diversity_failure = min(1.0, expected_diversity / warning_threshold)
    lead_lower = max(0.0, 1.0 - diversity_failure - trait_persistence_failure_upper_bound)
    return ProbabilisticLeadCertificate(
        time=time,
        expected_diversity_upper_bound=expected_diversity,
        diversity_warning_failure_upper_bound=diversity_failure,
        trait_persistence_failure_upper_bound=trait_persistence_failure_upper_bound,
        lead_probability_lower_bound=lead_lower,
    )


def finite_bin_h_alpha_lead_certificate(
    initial_diversity: float,
    warning_threshold: float,
    multiplier_upper_bound: float,
    time: int,
    cohort_size_lower_bound: int,
    high_trait_recruit_probability_lower_bound: float,
    occupancy_threshold: int,
) -> tuple[BinomialTraitPersistenceBound, ProbabilisticLeadCertificate]:
    """Combine a finite-bin Chernoff bound with the H-alpha lead certificate."""
    trait_bound = binomial_trait_persistence_bound(
        cohort_size_lower_bound,
        high_trait_recruit_probability_lower_bound,
        occupancy_threshold,
    )
    certificate = probabilistic_h_alpha_lead_certificate(
        initial_diversity=initial_diversity,
        warning_threshold=warning_threshold,
        multiplier_upper_bound=multiplier_upper_bound,
        time=time,
        trait_persistence_failure_upper_bound=trait_persistence_union_bound(
            trait_bound.per_generation_failure_upper_bound,
            time,
        ),
    )
    return trait_bound, certificate
