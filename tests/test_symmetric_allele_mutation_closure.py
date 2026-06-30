from dataclasses import replace

from causal_model.finite_h1_mutation_window_audit import run_finite_h1_mutation_window_audit
from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate
from causal_model.multipatch_criticality_experiments import quick_profile
from causal_model.symmetric_allele_mutation_closure import (
    apply_symmetric_allele_mutation,
    simulate_with_symmetric_allele_mutation,
)


def test_zero_mutation_delegates_to_identical_legacy_trajectory():
    parameters = DynamicsParameters(
        patch_areas=(1.0,),
        generations=4,
        initial_population=(20,),
        initial_interaction=(0.6,),
        initial_high_allele_frequency=(0.7,),
        random_seed=101,
    )
    assert simulate_with_symmetric_allele_mutation(parameters, mutation_rate=0.0) == simulate(parameters)


def test_symmetric_mutation_map_moves_fixed_states_inward_and_keeps_half_fixed():
    assert apply_symmetric_allele_mutation(1.0, 0.1) == 0.9
    assert apply_symmetric_allele_mutation(0.0, 0.1) == 0.1
    assert apply_symmetric_allele_mutation(0.5, 0.1) == 0.5


def test_tiny_mutation_window_audit_retains_every_seed_replicate():
    spec = replace(
        quick_profile(),
        total_area=4.0,
        patch_count=4,
        generations=1,
        replicates=1,
        master_seed=17,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.5,),
    )
    cells = run_finite_h1_mutation_window_audit(
        spec,
        mutation_rates=(0.0, 0.1),
        master_seeds=(101, 202),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        hold_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )
    assert len(cells) == 2
    assert {cell.mutation_rate for cell in cells} == {0.0, 0.1}
    assert all(len(cell.replicates) == 2 for cell in cells)
    assert all(cell.summary["denominators"]["total_seed_replicates"] == 2 for cell in cells)
