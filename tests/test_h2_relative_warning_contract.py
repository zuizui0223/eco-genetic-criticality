from causal_model.h2_relative_warning_contract import (
    DeteriorationCalibrationCandidate,
    RelativeWarningDefinition,
    compare_relative_warning,
    first_relative_decline_time,
    h2r_protocol_manifest,
    select_trait_loss_only_calibration,
)


def test_relative_warning_excludes_baseline_and_records_strict_lead():
    definition = RelativeWarningDefinition("H_alpha", 0.10)
    threshold, warning_time = first_relative_decline_time((0.50, 0.49, 0.45, 0.40), definition)
    assert threshold == 0.45
    assert warning_time == 2
    comparison = compare_relative_warning((0.50, 0.49, 0.45, 0.40), trait_loss_time=4, definition=definition)
    assert comparison.valid_pair is True
    assert comparison.warning_leads is True
    assert comparison.lead_time_trait_minus_warning == 2


def test_zero_baseline_and_missing_trait_loss_remain_ineligible_or_censored():
    definition = RelativeWarningDefinition("H_gamma", 0.05)
    ineligible = compare_relative_warning((0.0, 0.0, 0.0), trait_loss_time=2, definition=definition)
    assert ineligible.baseline_eligible is False
    assert ineligible.warning_leads is None
    censored = compare_relative_warning((0.6, 0.5, 0.4), trait_loss_time=None, definition=definition)
    assert censored.censored is True
    assert censored.warning_leads is None


def test_calibration_selects_only_from_all_seed_trait_loss_availability_without_warning_inputs():
    too_sparse = DeteriorationCalibrationCandidate(60, 0.15, (0.10, 0.20, 0.25, 0.20, 0.15))
    eligible = DeteriorationCalibrationCandidate(60, 0.30, (0.40, 0.45, 0.50, 0.55, 0.60))
    also_eligible = DeteriorationCalibrationCandidate(120, 0.45, (0.45, 0.50, 0.50, 0.50, 0.55))
    selected = select_trait_loss_only_calibration((too_sparse, eligible, also_eligible))
    assert selected.selected == also_eligible
    assert selected.reason.startswith("selected using trait-loss availability")
    no_selection = select_trait_loss_only_calibration((too_sparse,))
    assert no_selection.selected is None


def test_manifest_keeps_original_h2_unresolved_and_h1_h3_unchanged():
    manifest = h2r_protocol_manifest()
    assert manifest["original_proposition"]["status_after_stationary_mutation_primary_chain"] == "unresolved_due_to_right_censoring_not_false"
    assert manifest["invariants"]["h1_h3_modified"] is False
    assert manifest["calibration"]["selection_uses"] == "trait-loss event availability only"
