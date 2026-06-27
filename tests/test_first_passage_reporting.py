import pytest

from causal_model.first_passage_reporting import compare_first_passage_times


def test_comparison_separates_conditional_and_unconditional_lead_frequencies():
    comparison = compare_first_passage_times(
        warning_times=(1, 5, None, 4),
        reference_times=(3, 4, 2, None),
    )

    assert comparison.replicate_count == 4
    assert comparison.valid_pair_count == 2
    assert comparison.censored_pair_count == 2
    assert comparison.valid_pair_probability == pytest.approx(0.5)
    assert comparison.time_differences == (-2, 1)
    assert comparison.lead_count == 1
    assert comparison.conditional_lead_probability == pytest.approx(0.5)
    assert comparison.unconditional_observed_lead_fraction == pytest.approx(0.25)
    assert comparison.median_time_difference == pytest.approx(-0.5)


def test_no_valid_event_pair_is_not_reported_as_zero_conditional_lead_probability():
    comparison = compare_first_passage_times(
        warning_times=(None, 4),
        reference_times=(3, None),
    )

    assert comparison.valid_pair_count == 0
    assert comparison.censored_pair_count == 2
    assert comparison.conditional_lead_probability is None
    assert comparison.unconditional_observed_lead_fraction == 0.0
    assert comparison.time_differences == ()


def test_comparison_rejects_mismatched_or_empty_vectors():
    with pytest.raises(ValueError):
        compare_first_passage_times((1,), (1, 2))
    with pytest.raises(ValueError):
        compare_first_passage_times((), ())
