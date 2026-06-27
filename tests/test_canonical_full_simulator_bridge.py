import pytest

from causal_model.canonical_full_simulator_bridge import (
    canonical_full_simulator_bridge_certificate,
    canonical_full_simulator_parameters,
)
from causal_model.multipatch_criticality_dynamics import simulate


def test_full_simulator_reproduces_canonical_h1_trajectory_exactly_in_declared_limit():
    certificate = canonical_full_simulator_bridge_certificate(
        area=1.0,
        area_reference=1.0,
        feedback_strength=8.0,
        barrier=0.5,
        initial_interaction=0.01,
        generations=40,
    )

    assert certificate.exact_embedding_certified
    assert certificate.maximum_absolute_error == pytest.approx(0.0, abs=1e-15)
    assert certificate.canonical_interaction == pytest.approx(certificate.full_simulator_interaction, abs=1e-15)
    assert certificate.full_simulator_interaction[-1] < 0.1


def test_full_simulator_embedding_reaches_the_other_canonical_branch_from_high_initial_state():
    certificate = canonical_full_simulator_bridge_certificate(
        area=1.0,
        area_reference=1.0,
        feedback_strength=8.0,
        barrier=0.5,
        initial_interaction=0.99,
        generations=40,
        carrying_population=37,
    )

    assert certificate.exact_embedding_certified
    assert certificate.full_simulator_interaction[-1] > 0.9


def test_embedding_keeps_census_at_declared_carrying_population():
    parameters = canonical_full_simulator_parameters(
        area=2.5,
        area_reference=1.0,
        feedback_strength=7.0,
        barrier=0.4,
        initial_interaction=0.3,
        generations=12,
        carrying_population=73,
    )
    result = simulate(parameters)

    assert all(snapshot.population == (73,) for snapshot in result.snapshots)


def test_invalid_embedding_inputs_are_rejected():
    with pytest.raises(ValueError):
        canonical_full_simulator_parameters(
            area=0.0,
            area_reference=1.0,
            feedback_strength=8.0,
            barrier=0.5,
            initial_interaction=0.1,
            generations=10,
        )
