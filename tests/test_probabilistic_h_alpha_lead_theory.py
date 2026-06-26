import math

import pytest

from causal_model.probabilistic_h_alpha_lead_theory import (
    binomial_trait_persistence_bound,
    finite_bin_h_alpha_lead_certificate,
    probabilistic_h_alpha_lead_certificate,
    trait_persistence_union_bound,
)


def test_binomial_trait_persistence_bound_uses_chernoff_lower_tail() -> None:
    bound = binomial_trait_persistence_bound(
        cohort_size_lower_bound=100,
        high_trait_recruit_probability_lower_bound=0.5,
        occupancy_threshold=20,
    )

    assert bound.expected_high_trait_abundance_lower_bound == pytest.approx(50.0)
    assert bound.per_generation_failure_upper_bound == pytest.approx(math.exp(-9.0))


def test_binomial_trait_persistence_returns_trivial_bound_without_margin() -> None:
    bound = binomial_trait_persistence_bound(10, 0.2, occupancy_threshold=2)
    assert bound.expected_high_trait_abundance_lower_bound == pytest.approx(2.0)
    assert bound.per_generation_failure_upper_bound == 1.0


def test_union_bound_accumulates_finite_trait_loss_risk() -> None:
    assert trait_persistence_union_bound(0.1, 3) == pytest.approx(0.3)
    assert trait_persistence_union_bound(0.6, 2) == 1.0


def test_probabilistic_lead_certificate_combines_markov_and_trait_risk() -> None:
    certificate = probabilistic_h_alpha_lead_certificate(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.5,
        time=2,
        trait_persistence_failure_upper_bound=0.1,
    )

    assert certificate.expected_diversity_upper_bound == pytest.approx(0.2)
    assert certificate.diversity_warning_failure_upper_bound == pytest.approx(0.4)
    assert certificate.lead_probability_lower_bound == pytest.approx(0.5)


def test_finite_bin_certificate_reports_both_assumption_layers() -> None:
    trait, certificate = finite_bin_h_alpha_lead_certificate(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.5,
        time=2,
        cohort_size_lower_bound=100,
        high_trait_recruit_probability_lower_bound=0.5,
        occupancy_threshold=20,
    )

    assert trait.per_generation_failure_upper_bound < 0.001
    assert certificate.trait_persistence_failure_upper_bound < 0.001
    assert certificate.lead_probability_lower_bound > 0.59


def test_invalid_inputs_are_rejected() -> None:
    with pytest.raises(ValueError):
        binomial_trait_persistence_bound(0, 0.5, 1)
    with pytest.raises(ValueError):
        trait_persistence_union_bound(-0.1, 1)
    with pytest.raises(ValueError):
        probabilistic_h_alpha_lead_certificate(0.8, 0.5, 1.0, 2, 0.1)
