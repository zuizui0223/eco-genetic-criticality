from causal_model.h2_relative_warning_contract import (
    DEFAULT_RELATIVE_DECLINE_FRACTIONS,
    RelativeWarningDefinition,
    compare_relative_warning,
)
from causal_model.h2r_independent_relative_validation import (
    H2RValidationRecord,
    H2RValidationTrajectory,
    _summarise,
)
from causal_model.h2r_validation_domain import SELECTED_VALIDATION_DOMAIN, h2r_validation_domain_manifest


def _record(seed: int, *, trait_loss_time: int | None) -> H2RValidationRecord:
    alpha = (0.50, 0.50, 0.45, 0.40)
    gamma = (0.50, 0.49, 0.45, 0.40)
    comparisons = tuple(
        compare_relative_warning(
            alpha if diversity_id == "H_alpha" else gamma,
            trait_loss_time=trait_loss_time,
            definition=RelativeWarningDefinition(diversity_id, decline),
        )
        for diversity_id in ("H_alpha", "H_gamma")
        for decline in DEFAULT_RELATIVE_DECLINE_FRACTIONS
    )
    outcome = H2RValidationTrajectory(
        trait_loss_time_post_baseline=trait_loss_time,
        h_alpha_series=alpha,
        h_gamma_series=gamma,
        comparisons=comparisons,
    )
    return H2RValidationRecord(
        master_seed=seed,
        replicate_index=0,
        calibration_seed=seed + 10,
        h1_resolution_supported=True,
        h1_full_state_source_prepared=True,
        anchor_barrier=0.5,
        canonical_interval_width=1.0,
        projection_supported=True,
        baseline_realised_high_trait_present=True,
        trajectory_seed=seed + 100,
        barrier_first_generation=0.505,
        barrier_at_hold=0.65,
        outcome=outcome,
    )


def test_validation_domain_exactly_locks_calibration_selected_cell_and_schedule():
    domain = SELECTED_VALIDATION_DOMAIN
    assert (domain.mutation_rate, domain.area_reference, domain.interaction_feedback) == (0.10, 0.8, 6.0)
    assert domain.schedule.schedule_id == "ramp30_hold90_d0p15"
    assert domain.calibration_seed_block_trait_loss_probabilities == (0.50, 0.40, 0.40, 0.50, 0.50)
    manifest = h2r_validation_domain_manifest()
    assert manifest["selected_domain_count"] == 1
    assert manifest["selection_used_warning_outcomes"] is False


def test_validation_summary_retains_censored_trajectories_and_does_not_reselect_schedule():
    records = (_record(1, trait_loss_time=3), _record(2, trait_loss_time=None))
    summary = _summarise(records, (1, 2))
    alpha_summary = next(
        value
        for value in summary["endpoint_summaries"]
        if value["definition"] == {"diversity_id": "H_alpha", "relative_decline_fraction": 0.10, "require_positive_baseline": True}
    )
    assert alpha_summary["trajectory_available_count"] == 2
    assert alpha_summary["valid_pair_count"] == 1
    assert alpha_summary["censored_count"] == 1
    assert alpha_summary["warning_lead_count"] == 1
    assert summary["interpretation"]["selection_repeated"] is False
