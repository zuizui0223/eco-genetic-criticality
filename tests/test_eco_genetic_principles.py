from pathlib import Path

import pytest

from causal_model.eco_genetic_principles import (
    eco_genetic_ordering,
    feedback_contraction_bound,
    finite_transmission_strictly_erodes,
    partition_threshold_result,
    trait_mode_branch_result,
    transmission_moment_identity,
)


def test_g0_identity_requires_only_first_two_moments():
    result = transmission_moment_identity(
        post_selection_frequency=0.4,
        next_frequency_mean=0.4,
        next_frequency_variance=0.03,
    )
    assert result.heterozygosity_after_selection == pytest.approx(0.48)
    assert result.expected_heterozygosity_after_transmission == pytest.approx(0.42)
    assert result.expected_heterozygosity_loss == pytest.approx(0.06)
    assert result.unbiased
    assert result.finite_variance
    assert finite_transmission_strictly_erodes(result)


def test_g0_does_not_overclaim_under_biased_transmission():
    result = transmission_moment_identity(
        post_selection_frequency=0.4,
        next_frequency_mean=0.5,
        next_frequency_variance=0.01,
    )
    assert not result.unbiased
    assert not finite_transmission_strictly_erodes(result)


def test_zero_transmission_variance_does_not_create_drift_erosion():
    result = transmission_moment_identity(0.4, 0.4, 0.0)
    assert result.expected_heterozygosity_loss == pytest.approx(0.0)
    assert not finite_transmission_strictly_erodes(result)


def test_p0_certifies_uniqueness_only_under_strict_global_contraction():
    below = feedback_contraction_bound(1.0, 2.0, 0.4)
    boundary = feedback_contraction_bound(1.0, 2.0, 0.5)
    above = feedback_contraction_bound(1.0, 2.0, 0.6)

    assert below.global_lipschitz_bound == pytest.approx(0.8)
    assert below.uniqueness_certified
    assert not below.bistability_not_ruled_out
    assert boundary.global_lipschitz_bound == pytest.approx(1.0)
    assert not boundary.uniqueness_certified
    assert boundary.bistability_not_ruled_out
    assert above.global_lipschitz_bound == pytest.approx(1.2)
    assert not above.uniqueness_certified
    assert above.bistability_not_ruled_out


def test_trait_mode_lifting_requires_strict_branch_margin_sign_change():
    present = trait_mode_branch_result(0.2, 0.8, -0.1, 0.05)
    boundary = trait_mode_branch_result(0.2, 0.8, 0.0, 0.05)
    assert present.branch_dependent_mode
    assert not boundary.branch_dependent_mode


def test_partition_nonadditivity_does_not_follow_from_total_area_alone():
    result = partition_threshold_result((1.5, 1.5, 1.5), critical_patch_size=2.0)
    assert result.total_area == pytest.approx(4.5)
    assert result.total_area_exceeds_threshold
    assert not result.any_patch_exceeds_threshold
    assert not result.collective_mechanism_possible_by_threshold


def test_conditional_eco_genetic_coupling_requires_explicit_ordering():
    supported = eco_genetic_ordering(
        0.2, 0.8,
        low_effective_size=10.0,
        high_effective_size=30.0,
        low_transmission_variance=0.015,
        high_transmission_variance=0.005,
    )
    reversed_case = eco_genetic_ordering(
        0.2, 0.8,
        low_effective_size=10.0,
        high_effective_size=30.0,
        low_transmission_variance=0.005,
        high_transmission_variance=0.015,
    )
    assert supported.ordering_supported
    assert supported.low_expected_heterozygosity_loss == pytest.approx(0.03)
    assert supported.high_expected_heterozygosity_loss == pytest.approx(0.01)
    assert not reversed_case.ordering_supported


def test_hypothesis_program_keeps_central_hypotheses_outside_theorem_claims():
    text = Path("docs/eco_genetic_hypothesis_program.md").read_text(encoding="utf-8")
    assert "H_critical" in text
    assert "H_genetic_lag" in text
    assert "H_fragmentation" in text
    assert "simulation result" in text
    assert "not a theorem" in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        transmission_moment_identity(0.4, 0.4, -0.1)
    with pytest.raises(ValueError):
        feedback_contraction_bound(0.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        trait_mode_branch_result(0.8, 0.2, -1.0, 1.0)
    with pytest.raises(ValueError):
        partition_threshold_result((), 1.0)
