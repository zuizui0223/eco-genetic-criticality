from dataclasses import replace
from pathlib import Path

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate
from causal_model.multipatch_criticality_experiments import (
    PROFILE_FULL,
    PROFILE_QUICK,
    PROFILE_STANDARD,
    derived_seed,
    full_profile,
    parameter_grid,
    quick_profile,
    results_to_csv_rows,
    run_parameter_grid,
    scenario_equal_isolated,
    scenario_equal_migrating,
    scenario_one_large,
    standard_profile,
    summarise_replicate,
)
from examples.multipatch_phase_diagram_pilot import pilot_spec, render_markdown
from examples.multipatch_h_alpha_lead_phase_boundary import (
    LEAD_PROBABILITY_KEY,
    heatmap_matrix,
    render_markdown as render_h_alpha_boundary_markdown,
)


def _low_trait_distribution(size: int) -> tuple[float, ...]:
    return tuple(1.0 if index == 0 else 0.0 for index in range(size))


def _high_trait_distribution(size: int) -> tuple[float, ...]:
    return tuple(1.0 if index == size - 1 else 0.0 for index in range(size))


def test_scenario_constructors_preserve_total_area() -> None:
    spec = quick_profile()
    one_large = scenario_one_large(spec)
    isolated = scenario_equal_isolated(spec)
    migrating = scenario_equal_migrating(spec)
    assert one_large.total_area == spec.total_area
    assert sum(isolated.patch_areas) == spec.total_area
    assert sum(migrating.patch_areas) == spec.total_area


def test_isolated_and_migrating_differ_only_in_migration_settings() -> None:
    spec = replace(quick_profile(), migration_rate=0.25)
    isolated = scenario_equal_isolated(spec)
    migrating = scenario_equal_migrating(spec)
    assert isolated.patch_areas == migrating.patch_areas
    assert isolated.migration_rate == 0.0
    assert migrating.migration_rate == 0.25


def test_seed_schedule_is_reproducible_and_shared_across_scenarios() -> None:
    spec = quick_profile()
    cell = parameter_grid(spec)[1]
    assert derived_seed(spec.master_seed, cell.cell_index, 2) == derived_seed(spec.master_seed, cell.cell_index, 2)
    assert derived_seed(spec.master_seed, cell.cell_index, 2) != derived_seed(spec.master_seed, cell.cell_index, 1)


def test_experiment_result_keeps_censored_first_passages_and_metadata_explicit() -> None:
    spec = replace(quick_profile(), replicates=2)
    result = run_parameter_grid(spec, scenarios=(scenario_one_large(spec),))[0]
    summary = result.summary
    assert "censored_event_counts" in summary
    assert "tau_trait_realised" in summary["first_passage"]
    assert "tau_allele_loss" in summary["first_passage"]
    assert "aggregation_rule" in summary["first_passage"]["tau_allele_loss"]
    assert "valid_event_pair_counts" in summary
    assert all("events" in replicate.as_dict() for replicate in result.replicates)


def test_potential_and_realised_trait_outcomes_can_differ() -> None:
    size = 21
    params = DynamicsParameters(
        patch_areas=(1.0,),
        generations=1,
        initial_interaction=(1.0,),
        initial_trait_distribution=(_low_trait_distribution(size),),
        trait_grid_size=size,
        high_interaction_benefit=2.0,
        viability_threshold=1.0,
        realised_high_trait_threshold=1e-6,
        random_seed=3,
    )
    summary = summarise_replicate(
        simulate(params),
        replicate_index=0,
        seed=3,
        h_alpha_warning_threshold=0.0,
        h_gamma_warning_threshold=0.0,
        fst_warning_threshold=1.0,
    )
    assert summary.potential_high_trait_viable is True
    assert summary.realised_high_trait_persists is False


def test_allele_persistence_and_realised_trait_occupancy_can_differ() -> None:
    size = 21
    params = DynamicsParameters(
        patch_areas=(1.0,),
        generations=2,
        initial_high_allele_frequency=(0.5,),
        initial_trait_distribution=(_low_trait_distribution(size),),
        trait_grid_size=size,
        random_seed=5,
    )
    summary = summarise_replicate(
        simulate(params),
        replicate_index=0,
        seed=5,
        h_alpha_warning_threshold=0.0,
        h_gamma_warning_threshold=0.0,
        fst_warning_threshold=1.0,
    )
    assert summary.tau_trait_realised == 0
    assert summary.h_alpha > 0.0
    assert summary.tau_H_alpha is None
    assert summary.tau_allele_loss is None


def test_standard_profile_uses_finite_two_kernel_coupled_closure() -> None:
    spec = standard_profile()
    parameters = spec.base_parameters
    assert spec.profile == PROFILE_STANDARD
    assert parameters.trait_occupancy_mode == "finite_trait_bin_recruitment"
    assert parameters.genotype_trait_recruitment == "two_kernel_recruitment"
    assert parameters.q_feedback_beta_trait > 0.0
    assert parameters.q_feedback_gamma_allele is not None
    assert parameters.q_feedback_gamma_allele > 0.0


def test_finite_profile_reports_abundance_independently_of_mass() -> None:
    spec = replace(standard_profile(), generations=2, replicates=1, area_reference_values=(1.0,), interaction_feedback_values=(3.0,), interaction_barrier_values=(0.5,))
    result = run_parameter_grid(spec, scenarios=(scenario_one_large(spec),))[0]
    replicate = result.replicates[0]
    assert len(replicate.final_high_trait_abundance_by_patch) == 1
    assert "realised_high_trait_abundance" in result.summary["metrics"]
    row = results_to_csv_rows((result,))[0]
    assert "metrics.realised_high_trait_abundance.mean" in row


def test_h_alpha_h_gamma_and_fst_are_returned_independently() -> None:
    spec = replace(quick_profile(), replicates=2)
    result = run_parameter_grid(spec, scenarios=(scenario_equal_isolated(spec),))[0]
    metrics = result.summary["metrics"]
    assert set(metrics) >= {"H_alpha", "H_gamma", "F_ST"}
    row = results_to_csv_rows((result,))[0]
    assert "metrics.H_alpha.mean" in row
    assert "metrics.H_gamma.mean" in row
    assert "metrics.F_ST.mean" in row


def test_quick_profile_completes_quickly_and_remains_deterministic() -> None:
    spec = quick_profile()
    results = run_parameter_grid(spec)
    assert spec.profile == PROFILE_QUICK
    assert spec.base_parameters.trait_occupancy_mode == "deterministic_viability_selection"
    assert len(results) == len(parameter_grid(spec)) * 3


def test_full_phase_diagram_is_opt_in_not_triggered_by_quick_profile() -> None:
    quick = quick_profile()
    full = full_profile()
    assert quick.profile == PROFILE_QUICK
    assert full.profile == PROFILE_FULL
    assert quick.replicates < full.replicates
    assert len(parameter_grid(quick)) < len(parameter_grid(full))


def test_experiment_layer_does_not_reinterpret_theorem_claims() -> None:
    import causal_model.multipatch_criticality_experiments as experiments
    text = experiments.__doc__
    assert text is not None
    assert "simulation/reporting layer" in text
    assert "does not alter the theorem layer" in text


def test_replicate_summary_reports_patch_distributions() -> None:
    size = 21
    params = DynamicsParameters(
        patch_areas=(2.0, 2.0),
        generations=1,
        initial_trait_distribution=(_high_trait_distribution(size), _low_trait_distribution(size)),
        trait_grid_size=size,
        random_seed=7,
    )
    summary = summarise_replicate(
        simulate(params),
        replicate_index=0,
        seed=7,
        h_alpha_warning_threshold=0.2,
        h_gamma_warning_threshold=0.2,
        fst_warning_threshold=0.2,
    )
    assert len(summary.final_q_by_patch) == 2
    assert len(summary.final_population_by_patch) == 2
    assert len(summary.final_effective_size_by_patch) == 2
    assert len(summary.final_p_by_patch) == 2
    assert len(summary.final_high_trait_abundance_by_patch) == 2


def test_pilot_report_keeps_math_layer_before_ecosystem_projection() -> None:
    spec = replace(
        pilot_spec(),
        generations=1,
        replicates=1,
        area_reference_values=(0.8, 1.0, 1.2),
        interaction_feedback_values=(3.0, 4.5),
        interaction_barrier_values=(0.55, 0.75),
    )
    report = render_markdown(run_parameter_grid(spec))

    assert "model-specific stochastic pilot" in report
    assert "not an empirical ecosystem projection" in report
    assert "Omega_tau^potential != realised N_H > 0 != p > 0 != H_alpha > 0" in report
    assert "not to map them onto Campanula immediately" in report


def test_h_alpha_phase_boundary_uses_simulation_feedback_label() -> None:
    rows = (
        _boundary_row("equal_isolated", 1.2, 4.5, 0.75, 1.0),
        _boundary_row("one_large", 1.2, 4.5, 0.75, 0.0),
    )
    report = render_h_alpha_boundary_markdown(rows, Path("docs/figures/h_alpha_lead_phase_boundary.svg"))

    assert "interaction_feedback" in report
    assert "not the canonical logistic theorem parameter `kappa`" in report
    assert "Campanula mapping is deliberately deferred" in report


def test_h_alpha_phase_boundary_matrix_reads_lead_probability() -> None:
    rows = (
        _boundary_row("equal_isolated", 0.8, 3.0, 0.35, 0.25),
        _boundary_row("equal_isolated", 0.8, 3.0, 0.75, 0.50),
        _boundary_row("equal_isolated", 1.2, 3.0, 0.35, 0.75),
        _boundary_row("equal_isolated", 1.2, 3.0, 0.75, 1.00),
    )
    area_values, theta_values, matrix = heatmap_matrix(rows, scenario_id="equal_isolated", interaction_feedback=3.0)

    assert area_values == (0.8, 1.2)
    assert theta_values == (0.35, 0.75)
    assert matrix == [[0.25, 0.5], [0.75, 1.0]]


def _boundary_row(
    scenario_id: str,
    area_reference: float,
    interaction_feedback: float,
    theta: float,
    h_alpha_lead_probability: float,
) -> dict[str, object]:
    return {
        "scenario_id": scenario_id,
        "area_reference": area_reference,
        "interaction_feedback": interaction_feedback,
        "interaction_barrier": theta,
        LEAD_PROBABILITY_KEY: h_alpha_lead_probability,
    }
