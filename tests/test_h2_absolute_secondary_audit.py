from causal_model.h2_absolute_secondary_audit import (
    AbsoluteThresholdDefinition,
    audit_h2a_from_h2r_validation_payload,
    compare_fixed_absolute_warning,
    first_absolute_threshold_time,
)


def test_absolute_threshold_time_includes_existing_generation_zero_rule():
    assert first_absolute_threshold_time((0.19, 0.30), 0.20) == 0
    assert first_absolute_threshold_time((0.30, 0.20, 0.10), 0.20) == 1
    assert first_absolute_threshold_time((0.30, 0.21), 0.20) is None


def test_fixed_absolute_comparison_keeps_lead_tie_lag_and_censoring_separate():
    definition = AbsoluteThresholdDefinition("H_alpha")
    lead = compare_fixed_absolute_warning((0.6, 0.2, 0.1), trait_loss_time=3, definition=definition)
    tie = compare_fixed_absolute_warning((0.6, 0.2, 0.1), trait_loss_time=1, definition=definition)
    lag = compare_fixed_absolute_warning((0.6, 0.2, 0.1), trait_loss_time=0, definition=definition)
    censored = compare_fixed_absolute_warning((0.6, 0.4), trait_loss_time=2, definition=definition)
    assert (lead.warning_leads, tie.warning_ties, lag.warning_lags) == (True, True, True)
    assert censored.censored is True
    assert censored.warning_leads is None


def test_secondary_audit_reports_mixed_ordering_without_new_selection():
    payload = {
        "records": [
            {
                "master_seed": 1,
                "replicate_index": 0,
                "calibration_seed": 11,
                "outcome": {
                    "trait_loss_time_post_baseline": 3,
                    "h_alpha_series": [0.6, 0.2, 0.1, 0.1],
                    "h_gamma_series": [0.6, 0.5, 0.4, 0.1],
                },
            },
            {
                "master_seed": 2,
                "replicate_index": 0,
                "calibration_seed": 22,
                "outcome": {
                    "trait_loss_time_post_baseline": 1,
                    "h_alpha_series": [0.6, 0.5, 0.2],
                    "h_gamma_series": [0.6, 0.5, 0.2],
                },
            },
            {
                "master_seed": 3,
                "replicate_index": 0,
                "calibration_seed": 33,
                "outcome": None,
            },
        ]
    }
    audit = audit_h2a_from_h2r_validation_payload(payload)
    alpha = next(item for item in audit["endpoint_summaries"] if item["definition"]["diversity_id"] == "H_alpha")
    assert audit["denominators"] == {
        "attempted_source_records": 3,
        "trajectory_available_count": 2,
        "trajectory_unavailable_count": 1,
    }
    assert alpha["valid_pair_count"] == 2
    assert alpha["warning_lead_count"] == 1
    assert alpha["warning_lag_count"] == 1
    assert audit["audit"]["new_simulation"] is False
    assert audit["canonical_interpretation"]["h2a_global_truth_value"] == "not_assigned"
