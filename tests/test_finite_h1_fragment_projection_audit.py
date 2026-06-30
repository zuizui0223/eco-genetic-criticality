from dataclasses import replace

from causal_model.finite_h1_fragment_projection_audit import (
    FullState,
    _contingency_allocate,
    _positive_largest_remainder,
    project_full_state,
    run_finite_h1_fragment_projection_audit,
)
from causal_model.multipatch_criticality_dynamics import DynamicsParameters
from causal_model.multipatch_criticality_experiments import quick_profile


def test_positive_largest_remainder_preserves_total_and_patch_positivity():
    allocation = _positive_largest_remainder(10, (1.0, 1.0, 1.0, 1.0))
    assert sum(allocation) == 10
    assert all(value >= 1 for value in allocation)


def test_contingency_allocation_preserves_every_row_and_trait_bin_margin():
    matrix = _contingency_allocate((3, 4, 3), (2, 5, 3))
    assert tuple(sum(row) for row in matrix) == (3, 4, 3)
    assert tuple(sum(matrix[row][column] for row in range(3)) for column in range(3)) == (2, 5, 3)


def test_projected_initial_snapshot_preserves_declared_absolute_and_intensive_states():
    source = FullState(
        patch_areas=(4.0,),
        population=(10,),
        interaction=(0.8,),
        high_allele_frequency=(0.7,),
        trait_abundance=((2, 3, 5),),
    )
    target = DynamicsParameters(
        patch_areas=(1.0, 1.0, 1.0, 1.0),
        generations=1,
        trait_grid_size=3,
        density_capacity=40.0,
    )
    parameters, invariants = project_full_state(source, target)
    assert sum(parameters.initial_population) == 10
    assert tuple(sum(row) for row in parameters.initial_trait_abundance) == parameters.initial_population
    assert invariants.projection_supported is True
    assert invariants.target_total_population == 10
    assert invariants.target_trait_bin_totals == (2, 3, 5)
    assert invariants.target_weighted_interaction == 0.8
    assert invariants.target_weighted_high_allele_frequency == 0.7


def test_tiny_projection_audit_retains_declared_seed_replicates():
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
    cells = run_finite_h1_fragment_projection_audit(
        spec,
        master_seeds=(101, 202),
        endpoint_padding_fraction=0.5,
        stage_generations=1,
        hold_generations=1,
        nested_barrier_points=(5, 9),
        interaction_separation_threshold=0.0,
        maximum_normalized_bracket_width=1.0,
    )
    assert len(cells) == 1
    assert len(cells[0].replicates) == 2
    assert cells[0].summary["denominators"]["total_seed_replicates"] == 2
