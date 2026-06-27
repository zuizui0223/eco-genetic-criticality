import pytest

from causal_model.network_h3_experiments import simulate_h3_ensemble
from causal_model.network_h3_lifecycle import NetworkLifecycleParameters, PatchState


def test_ensemble_reports_recolonisation_and_preserves_per_replicate_event_times():
    parameters = NetworkLifecycleParameters(
        capacities=(40, 40),
        source_to_destination_kernel=((0.0, 1.0), (0.0, 1.0)),
        generations=2,
        adult_survival_probability=1.0,
        emigration_probability=1.0,
        recruitment_per_adult=0.0,
        persistence_threshold=1,
        colonisation_threshold=1,
        random_seed=3,
    )
    _, summary = simulate_h3_ensemble(
        parameters,
        (PatchState(10, 10, 20), PatchState(0, 0, 0)),
        replicates=3,
    )

    assert len(summary.replicates) == 3
    assert summary.recolonisation_probability == pytest.approx(1.0)
    assert summary.metapopulation_extinction_probability == pytest.approx(0.0)
    assert all(replicate.recolonisation_time == 1 for replicate in summary.replicates)
    assert all(replicate.metapopulation_extinction_time is None for replicate in summary.replicates)


def test_ensemble_rejects_zero_replicates():
    parameters = NetworkLifecycleParameters(
        capacities=(10,),
        source_to_destination_kernel=((1.0,),),
    )
    with pytest.raises(ValueError):
        simulate_h3_ensemble(parameters, (PatchState(1, 0, 0),), replicates=0)
