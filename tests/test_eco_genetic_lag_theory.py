import pytest

from causal_model.eco_genetic_lag_theory import (
    assess_genetic_lag,
    certify_trait_persistence_bound,
    conditional_lead_certificate,
    cumulative_multiplier,
    diversity_trajectory,
    exact_lead_condition,
    first_warning_time,
    uniform_upper_multiplier_bound,
)


def test_genetic_lead_is_exact_first_passage_before_trait_collapse():
    assessment = assess_genetic_lag(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multipliers=(0.9, 0.9, 0.9, 0.9),
        trait_collapse_time=4,
    )
    assert assessment.diversity_trajectory == pytest.approx((0.8, 0.72, 0.648, 0.5832, 0.52488))
    assert assessment.warning_time is None
    assert not assessment.genetic_lead


def test_no_lead_and_lead_exist_in_same_recursion_family():
    no_lead = assess_genetic_lag(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multipliers=(0.99, 0.99, 0.99, 0.99),
        trait_collapse_time=4,
    )
    lead = assess_genetic_lag(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multipliers=(0.7, 0.7, 0.95, 0.95),
        trait_collapse_time=4,
    )
    assert no_lead.warning_time is None
    assert not no_lead.genetic_lead
    assert lead.warning_time == 2
    assert lead.genetic_lead
    assert exact_lead_condition(0.8, 0.5, (0.7, 0.7, 0.95, 0.95), 4)


def test_warning_at_collapse_is_not_a_genetic_lead():
    assessment = assess_genetic_lag(
        initial_diversity=0.8,
        warning_threshold=0.4,
        multipliers=(0.9, 0.9, 0.9, 0.9),
        trait_collapse_time=4,
    )
    assert assessment.warning_time is None
    assert not assessment.genetic_lead


def test_cumulative_multiplier_uses_empty_product_at_time_zero():
    multipliers = (0.8, 0.9, 0.7)
    assert cumulative_multiplier(multipliers, 0) == pytest.approx(1.0)
    assert cumulative_multiplier(multipliers, 2) == pytest.approx(0.72)


def test_uniform_bound_gives_sufficient_not_necessary_lead_certificate():
    guaranteed = uniform_upper_multiplier_bound(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.7,
        trait_collapse_time=4,
    )
    no_guarantee = uniform_upper_multiplier_bound(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.95,
        trait_collapse_time=4,
    )
    assert guaranteed.latest_guaranteed_warning_time == 2
    assert guaranteed.lead_guaranteed
    assert no_guarantee.latest_guaranteed_warning_time > 4
    assert not no_guarantee.lead_guaranteed


def test_trait_persistence_bound_returns_first_possible_loss_time():
    bound = certify_trait_persistence_bound((5.0, 3.0, 0.1, 0.0, 0.0))

    assert bound.certified_trait_loss_lower_bound == 3
    assert bound.extinction_threshold == 0.0


def test_conditional_lead_certificate_combines_decay_and_persistence_bounds():
    certificate = conditional_lead_certificate(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.7,
        trait_lower_bounds=(4.0, 3.0, 2.0, 1.0, 0.5),
    )

    assert certificate.diversity_bound.latest_guaranteed_warning_time == 2
    assert certificate.trait_persistence_bound.certified_trait_loss_lower_bound == 5
    assert certificate.lead_guaranteed


def test_conditional_lead_certificate_is_not_triggered_without_trait_margin():
    certificate = conditional_lead_certificate(
        initial_diversity=0.8,
        warning_threshold=0.5,
        multiplier_upper_bound=0.7,
        trait_lower_bounds=(4.0, 1.0, 0.0, 0.0),
    )

    assert certificate.diversity_bound.latest_guaranteed_warning_time == 2
    assert certificate.trait_persistence_bound.certified_trait_loss_lower_bound == 2
    assert not certificate.lead_guaranteed


def test_first_warning_time_and_inputs():
    assert first_warning_time((0.8, 0.6, 0.5, 0.4), 0.5) == 2
    assert diversity_trajectory(0.6, (0.9, 0.9)) == pytest.approx((0.6, 0.54, 0.486))
    with pytest.raises(ValueError):
        assess_genetic_lag(0.8, 0.5, (0.9,), trait_collapse_time=2)
    with pytest.raises(ValueError):
        uniform_upper_multiplier_bound(0.8, 0.5, 1.0, 3)
    with pytest.raises(ValueError):
        certify_trait_persistence_bound(())
