import pytest

from causal_model.network_h3_lifecycle import (
    NetworkLifecycleParameters,
    PatchState,
    first_metapopulation_extinction,
    first_realised_high_trait_loss,
    first_recolonisation,
    isolated_kernel,
    simulate_network_lifecycle,
)


def _parameters(kernel, **overrides):
    values = {
        "capacities": (120, 120),
        "source_to_destination_kernel": kernel,
        "generations": 1,
        "adult_survival_probability": 1.0,
        "emigration_probability": 1.0,
        "recruitment_per_adult": 0.0,
        "high_trait_recruitment_multiplier": 1.0,
        "persistence_threshold": 1,
        "colonisation_threshold": 1,
        "random_seed": 17,
    }
    values.update(overrides)
    return NetworkLifecycleParameters(**values)


def test_isolation_prevents_recolonisation_of_an_empty_patch():
    result = simulate_network_lifecycle(
        _parameters(isolated_kernel(2)),
        (PatchState(20, 20, 40), PatchState(0, 0, 0)),
    )
    destination = result.snapshots[1]

    assert destination.states[0] == PatchState(20, 20, 40)
    assert destination.states[1] == PatchState(0, 0, 0)
    assert destination.transitions is not None
    assert destination.transitions[1].status == "empty"
    assert first_recolonisation(result) is None


def test_directed_individual_dispersal_recolonises_demography_trait_and_allele_state():
    directed = ((0.0, 1.0), (0.0, 1.0))
    result = simulate_network_lifecycle(
        _parameters(directed),
        (PatchState(20, 20, 40), PatchState(0, 0, 0)),
    )
    destination = result.snapshots[1]

    assert destination.states[0] == PatchState(0, 0, 0)
    assert destination.states[1] == PatchState(20, 20, 40)
    assert destination.transitions is not None
    assert destination.transitions[1].status == "recolonised"
    assert first_recolonisation(result) == 1


def test_incoming_individuals_can_rescue_a_patch_that_falls_below_persistence_threshold():
    source_to_destination = ((0.75, 0.25), (0.0, 1.0))
    result = simulate_network_lifecycle(
        _parameters(
            source_to_destination,
            persistence_threshold=5,
            colonisation_threshold=5,
            random_seed=14,
        ),
        (PatchState(100, 100, 200), PatchState(1, 0, 0)),
    )
    destination = result.snapshots[1]

    assert destination.transitions is not None
    assert destination.transitions[1].status == "rescued"
    assert destination.transitions[1].inbound_migrants >= 5
    assert destination.states[1].population >= 5
    assert destination.states[1].high_trait_abundance > 0
    assert destination.states[1].high_allele_copies > 0


def test_subthreshold_arrival_does_not_create_a_spurious_colony():
    directed = ((0.0, 1.0), (0.0, 1.0))
    result = simulate_network_lifecycle(
        _parameters(directed, colonisation_threshold=2),
        (PatchState(1, 1, 2), PatchState(0, 0, 0)),
    )
    destination = result.snapshots[1]

    assert destination.states[1] == PatchState(0, 0, 0)
    assert destination.transitions is not None
    assert destination.transitions[1].status == "empty"


def test_local_extinction_and_trait_loss_are_recorded_without_migration():
    result = simulate_network_lifecycle(
        _parameters(
            isolated_kernel(2),
            adult_survival_probability=0.0,
            emigration_probability=0.0,
            generations=2,
        ),
        (PatchState(5, 5, 10), PatchState(4, 0, 2)),
    )

    assert first_metapopulation_extinction(result) == 1
    assert first_realised_high_trait_loss(result) == 1
    assert result.snapshots[1].h_alpha == 0.0
    assert result.snapshots[1].h_gamma == 0.0
    assert result.snapshots[1].fst is None


def test_recruitment_retains_finite_trait_and_allele_counts_within_capacity():
    parameters = NetworkLifecycleParameters(
        capacities=(30,),
        source_to_destination_kernel=((1.0,),),
        generations=2,
        adult_survival_probability=1.0,
        emigration_probability=0.0,
        recruitment_per_adult=1.0,
        high_trait_recruitment_multiplier=2.0,
        persistence_threshold=1,
        colonisation_threshold=1,
        random_seed=9,
    )
    result = simulate_network_lifecycle(parameters, (PatchState(10, 10, 20),))

    for snapshot in result.snapshots:
        state = snapshot.states[0]
        assert 0 <= state.high_trait_abundance <= state.population <= 30
        assert 0 <= state.high_allele_copies <= 2 * state.population
    assert result.snapshots[-1].states[0].population > 10


def test_invalid_kernel_and_initial_population_above_capacity_are_rejected():
    with pytest.raises(ValueError):
        NetworkLifecycleParameters(
            capacities=(10, 10),
            source_to_destination_kernel=((0.7, 0.2), (0.3, 0.7)),
        )
    parameters = _parameters(isolated_kernel(2), capacities=(10, 10))
    with pytest.raises(ValueError):
        simulate_network_lifecycle(parameters, (PatchState(11, 0, 0), PatchState(0, 0, 0)))
