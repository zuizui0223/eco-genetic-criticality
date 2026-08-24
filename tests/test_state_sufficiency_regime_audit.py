from dataclasses import replace

import pytest

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate


def _base() -> DynamicsParameters:
    return DynamicsParameters(
        patch_areas=(1.0, 1.0),
        generations=1,
        initial_population=(20, 20),
        initial_interaction=(0.8, 0.2),
        initial_high_allele_frequency=(0.8, 0.2),
        random_seed=20260824,
    )


def _weighted_mean(values, weights):
    total = sum(weights)
    return sum(v * w for v, w in zip(values, weights)) / total


def test_exact_present_state_is_future_sufficient_under_declared_closure():
    """With identical full state, forcing parameters and RNG seed, history adds no input."""
    parameters = _base()
    first = simulate(parameters)
    second = simulate(parameters)
    assert first == second


def test_coarse_aggregate_state_is_not_future_sufficient_when_patch_alignment_differs():
    """Same marginals and aggregate genetic summaries can hide different joint patch states."""
    aligned = _base()
    anti_aligned = replace(aligned, initial_high_allele_frequency=(0.2, 0.8))

    aligned_result = simulate(aligned)
    anti_result = simulate(anti_aligned)
    a0 = aligned_result.snapshots[0]
    b0 = anti_result.snapshots[0]

    # Same census, same q and p marginals, same weighted means, and same standard
    # aggregate diversity summaries. Only the patchwise q-p alignment differs.
    assert a0.population == b0.population == (20, 20)
    assert sorted(a0.interaction) == sorted(b0.interaction) == [0.2, 0.8]
    assert sorted(a0.high_allele_frequency) == sorted(b0.high_allele_frequency) == [0.2, 0.8]
    assert _weighted_mean(a0.interaction, a0.population) == pytest.approx(
        _weighted_mean(b0.interaction, b0.population)
    )
    assert _weighted_mean(a0.high_allele_frequency, a0.population) == pytest.approx(
        _weighted_mean(b0.high_allele_frequency, b0.population)
    )
    assert a0.h_alpha == pytest.approx(b0.h_alpha)
    assert a0.h_gamma == pytest.approx(b0.h_gamma)
    assert a0.fst == pytest.approx(b0.fst)
    assert [x.high_trait_mass for x in a0.trait_occupancy] == pytest.approx(
        [x.high_trait_mass for x in b0.trait_occupancy]
    )

    # The declared interaction update is patchwise. Positive q-p alignment
    # therefore produces a different next interaction field than anti-alignment,
    # despite identical coarse summaries at t=0.
    a1 = aligned_result.snapshots[1]
    b1 = anti_result.snapshots[1]
    assert a1.interaction != pytest.approx(b1.interaction)
    assert sorted(a1.interaction) != pytest.approx(sorted(b1.interaction))
