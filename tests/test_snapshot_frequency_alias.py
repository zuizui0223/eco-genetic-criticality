from causal_model.multipatch_criticality_dynamics import SimulationSnapshot


def test_generic_frequency_alias_matches_high_allele_frequency():
    snapshot = SimulationSnapshot(
        generation=0,
        interaction=(0.5,),
        population=(10,),
        effective_size=(5.0,),
        high_allele_frequency=(0.25,),
        trait_space=(),
        trait_occupancy=(),
        h_alpha=0.0,
        h_gamma=0.0,
        fst=None,
    )
    assert snapshot.allele_frequency == snapshot.high_allele_frequency
