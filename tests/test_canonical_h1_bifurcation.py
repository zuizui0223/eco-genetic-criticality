import pytest

from causal_model.canonical_h1_bifurcation import (
    canonical_bifurcation_certificate,
    canonical_bistable_barrier_interval,
    canonical_h1_certificate,
    follow_barrier_path,
    iterate_canonical_map,
)
from causal_model.multipatch_criticality_dynamics import DynamicsParameters


@pytest.fixture
def canonical_kwargs():
    return {
        "feedback_strength": 8.0,
        "area": 1.0,
        "area_reference": 1.0,
        "barrier": 0.5,
    }


def test_gain_at_or_below_four_has_no_strict_bistable_interval():
    assert canonical_bistable_barrier_interval(4.0, 1.0) is None
    certificate = canonical_bifurcation_certificate(
        feedback_strength=4.0,
        area=1.0,
        barrier=0.5,
    )
    assert certificate.regime == "monostable"
    assert not certificate.strict_bistability_certified
    assert len(certificate.fixed_points) == 1
    assert certificate.fixed_points[0].interaction == pytest.approx(0.5, abs=1e-8)


def test_exact_barrier_interval_certifies_three_fixed_points(canonical_kwargs):
    interval = canonical_bistable_barrier_interval(
        canonical_kwargs["feedback_strength"],
        canonical_kwargs["area"],
        canonical_kwargs["area_reference"],
    )
    assert interval is not None
    assert interval[0] == pytest.approx(0.3667900062)
    assert interval[1] == pytest.approx(0.6332099938)

    certificate = canonical_bifurcation_certificate(**canonical_kwargs)
    assert certificate.regime == "bistable"
    assert certificate.strict_bistability_certified
    assert len(certificate.fixed_points) == 3
    assert [point.stability for point in certificate.fixed_points] == [
        "stable",
        "unstable",
        "stable",
    ]
    assert certificate.fixed_points[0].interaction < 0.1
    assert certificate.fixed_points[1].interaction == pytest.approx(0.5, abs=1e-8)
    assert certificate.fixed_points[2].interaction > 0.9


def test_h1_certificate_requires_bistability_and_branch_margin_sign_change(canonical_kwargs):
    trait_parameters = DynamicsParameters(patch_areas=(1.0,))
    certificate = canonical_h1_certificate(
        **canonical_kwargs,
        trait_parameters=trait_parameters,
    )
    assert certificate.bifurcation.strict_bistability_certified
    assert certificate.low_stable_branch is not None
    assert certificate.high_stable_branch is not None
    assert certificate.low_stable_branch.high_trait_margin < 0.0
    assert certificate.high_stable_branch.high_trait_margin > 0.0
    assert certificate.branch_dependent_high_trait_mode

    no_high_trait_gain = canonical_h1_certificate(
        **canonical_kwargs,
        trait_parameters=DynamicsParameters(
            patch_areas=(1.0,),
            high_interaction_benefit=0.0,
        ),
    )
    assert not no_high_trait_gain.branch_dependent_high_trait_mode


def test_initial_conditions_converge_to_opposite_stable_branches(canonical_kwargs):
    low = iterate_canonical_map(0.01, **canonical_kwargs, iterations=100)
    high = iterate_canonical_map(0.99, **canonical_kwargs, iterations=100)
    assert low.terminal_interaction < 0.1
    assert high.terminal_interaction > 0.9


def test_upward_and_downward_barrier_paths_show_hysteresis():
    upward_barriers = (0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.8)
    descending_barriers = tuple(reversed(upward_barriers))
    common = {
        "feedback_strength": 8.0,
        "area": 1.0,
        "area_reference": 1.0,
        "iterations_per_barrier": 150,
    }
    upward = follow_barrier_path(
        upward_barriers,
        initial_interaction=0.999,
        **common,
    )
    downward = follow_barrier_path(
        descending_barriers,
        initial_interaction=0.001,
        **common,
    )
    up_at_middle = upward.terminal_interactions[upward.barriers.index(0.5)]
    down_at_middle = downward.terminal_interactions[downward.barriers.index(0.5)]

    assert up_at_middle > 0.9
    assert down_at_middle < 0.1
    assert upward.terminal_interactions[-1] < 0.1
    assert downward.terminal_interactions[-1] > 0.9


def test_invalid_arguments_are_rejected():
    with pytest.raises(ValueError):
        canonical_bifurcation_certificate(
            feedback_strength=0.0,
            area=1.0,
            barrier=0.5,
        )
    with pytest.raises(ValueError):
        follow_barrier_path(
            (),
            initial_interaction=0.5,
            feedback_strength=8.0,
            area=1.0,
        )
