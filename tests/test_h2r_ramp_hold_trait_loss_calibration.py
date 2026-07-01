import pytest

from causal_model.h2r_ramp_hold_trait_loss_calibration import (
    RampHoldCalibrationRecord,
    RampHoldSchedule,
    _select_schedule,
    _summarise_schedules,
    ramp_and_hold_barrier_schedule,
)


def _record(
    seed: int,
    replicate_index: int,
    schedule: RampHoldSchedule,
    *,
    trait_loss_time: int | None,
    baseline_present: bool = True,
) -> RampHoldCalibrationRecord:
    return RampHoldCalibrationRecord(
        mutation_rate=0.15,
        area_reference=0.8,
        interaction_feedback=6.0,
        master_seed=seed,
        replicate_index=replicate_index,
        calibration_seed=seed * 100 + replicate_index,
        schedule=schedule,
        h1_resolution_supported=True,
        h1_full_state_source_prepared=True,
        anchor_barrier=0.4,
        canonical_interval_width=0.8,
        projection_supported=True,
        baseline_realised_high_trait_present=baseline_present,
        trajectory_seed=seed * 1000 + replicate_index,
        barrier_first_generation=0.45,
        barrier_at_hold=0.60,
        trait_loss_time_post_baseline=trait_loss_time,
    )


def test_ramp_hold_schedule_reaches_final_barrier_then_stays_there():
    schedule = RampHoldSchedule(ramp_generations=4, hold_generations=3, total_normalized_barrier_increase=0.25)
    barriers = ramp_and_hold_barrier_schedule(
        anchor_barrier=0.4,
        canonical_interval_width=0.8,
        schedule=schedule,
    )
    assert barriers == pytest.approx((0.45, 0.50, 0.55, 0.60, 0.60, 0.60, 0.60))
    assert barriers[0] > 0.4
    assert len(barriers) == schedule.total_generations


def test_ramp_hold_summary_keeps_seed_blocks_and_excludes_baseline_ineligible_records():
    schedule = RampHoldSchedule(30, 90, 0.30)
    records = (
        _record(1, 0, schedule, trait_loss_time=40),
        _record(1, 1, schedule, trait_loss_time=None),
        _record(2, 0, schedule, trait_loss_time=35),
        _record(2, 1, schedule, trait_loss_time=None, baseline_present=False),
    )
    summary = _summarise_schedules(records, (1, 2), (schedule,))[0]
    assert summary.seed_block_baseline_eligible_counts == (2, 1)
    assert summary.seed_block_trait_loss_counts == (1, 1)
    assert summary.seed_block_trait_loss_probabilities == (0.5, 1.0)
    assert summary.pooled_trait_loss_probability == pytest.approx(2 / 3)
    assert "H_alpha" not in summary.as_dict()
    assert "H_gamma" not in summary.as_dict()


def test_selection_uses_trait_loss_only_and_returns_full_ramp_hold_identity():
    schedule = RampHoldSchedule(30, 90, 0.30)
    records = tuple(
        _record(seed, replicate, schedule, trait_loss_time=20 if replicate in {0, 1} else None)
        for seed in (1, 2, 3, 4, 5)
        for replicate in range(5)
    )
    summary = _summarise_schedules(records, (1, 2, 3, 4, 5), (schedule,))[0]
    selection = _select_schedule((summary,))
    assert selection["selection_uses_warning_outcomes"] is False
    assert selection["selected_schedule"] == schedule.as_dict()
