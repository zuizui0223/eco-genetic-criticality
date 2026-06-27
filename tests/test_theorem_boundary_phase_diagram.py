from dataclasses import replace

from causal_model.canonical_full_simulator_bridge import canonical_full_simulator_parameters
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    quick_profile,
    scenario_equal_migrating,
    scenario_one_large,
)
from causal_model.theorem_boundary_phase_diagram import (
    run_theorem_boundary_phase_diagram,
    write_theorem_boundary_phase_artifacts,
)


def test_boundary_phase_diagram_certifies_the_canonical_one_patch_limit():
    base = canonical_full_simulator_parameters(
        area=1.0,
        area_reference=1.0,
        feedback_strength=8.0,
        barrier=0.5,
        initial_interaction=0.01,
        generations=4,
        carrying_population=50,
    )
    spec = ExperimentSpec(
        experiment_id="canonical_limit_test",
        total_area=1.0,
        patch_count=1,
        generations=4,
        replicates=2,
        area_reference_values=(1.0,),
        interaction_feedback_values=(8.0,),
        interaction_barrier_values=(0.5,),
        migration_rate=0.0,
        base_parameters=base,
    )
    cells = run_theorem_boundary_phase_diagram(spec, scenarios=(scenario_one_large(spec),))

    assert len(cells) == 1
    cell = cells[0]
    assert all(rep.audit.single_patch_canonical_theorem_limit_certified for rep in cell.replicates)
    assert cell.summary["scope"]["single_patch_canonical_theorem_limit_probability"] == 1.0
    assert cell.summary["scope"]["maximum_canonical_update_residual"]["maximum"] == 0.0


def test_boundary_phase_diagram_records_named_departures_and_writes_artifacts(tmp_path):
    spec = replace(
        quick_profile(),
        total_area=2.0,
        patch_count=2,
        generations=3,
        replicates=2,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.4,),
        migration_rate=0.2,
    )
    cells = run_theorem_boundary_phase_diagram(spec, scenarios=(scenario_equal_migrating(spec),))
    cell = cells[0]
    departure = cell.summary["scope"]["departure_probabilities"]

    assert departure["migration_enabled"] == 1.0
    assert departure["multiple_patches"] == 1.0
    assert cell.summary["scope"]["single_patch_canonical_theorem_limit_probability"] == 0.0

    csv_path = tmp_path / "theorem_boundary.csv"
    json_path = tmp_path / "theorem_boundary.json"
    write_theorem_boundary_phase_artifacts(cells, csv_path=csv_path, json_path=json_path)

    assert csv_path.exists()
    assert json_path.exists()
    assert "scope.maximum_canonical_update_residual.mean" in csv_path.read_text(encoding="utf-8")
    assert "h1_theorem_boundary" in json_path.read_text(encoding="utf-8")
