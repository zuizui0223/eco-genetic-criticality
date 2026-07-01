import pytest

from causal_model.h2r_trait_loss_calibration import (
    H2RTraitLossCalibrationRecord,
    _select_schedule,
    _summarise_schedules,
    linear_normalized_barrier_schedule,
)


def _record(
    seed: int,
    replicate_index: int,
    *,
    trait_loss_time: int | None,
    baseline_present: bool = True,
) -> H2RTraitLossCalibrationRecord:
    return H2RTraitLossCalibrationRecord(
        mutation_rate=0.15,
        area_reference=0.8,
        interaction_feedback=6.0,
        master_seed=seed,
        replicate_index=replicate_index,
        calibration_seed=seed + 100 + replicate_index,
        horizon=60,
        total_normalized_barrier_increase=0.30,
        h1_resolution_supported=True,
        h1_full_state_source_prepared=True,
        anchor_barrier=0.5,
        canonical_interval_width=1.0,
        projection_supported=True,
        baseline_realised_high_trait_present=baseline_present,
        trajectory_seed=seed + 1000 + replicate_index,
        barrier_first_generation=0.505,
        barrier_final_generation=0.80,
        trait_loss_time_post_baseline=trait_loss_time,
    )


def test_linear_normalized_schedule_starts_after_anchor_and_ends_at_declared_increase():
    schedule = linear_normalized_barrier_schedule(
        anchor_barrier=0.4,
        canonical_interval_width=0.8,
        total_normalized_barrier_increase=0.25,
        horizon=4,
    )
    assert schedule == pytest.approx((0.45, 0.50, 0.55, 0.60))
    with pytest.raises(ValueError, match="positive"):
        linear_normalized_barrier_schedule(
            anchor_barrier=0.4,
            canonical_interval_width=0.8,
            total_normalized_barrier_increase=0.0,
            horizon=4,
        )


def test_schedule_summary_uses_post_baseline_trait_loss_only_and_keeps_seed_blocks_separate():
    records = (
        _record(1, 0, trait_loss_time=20),
        _record(1, 1, trait_loss_time=None),
        _record(2, 0, trait_loss_time=15),
        _record(2, 1, trait_loss_time=None, baseline_present=False),
    )
    summary = _summarise_schedules(records, (1, 2), (60,), (0.30,))[0]
    assert summary.seed_block_baseline_eligible_counts == (2, 1)
    assert summary.seed_block_trait_loss_counts == (1, 1)
    assert summary.seed_block_trait_loss_probabilities == (0.5, 1.0)
    assert summary.pooled_trait_loss_probability == pytest.approx(2 / 3)
    assert "H_alpha" not in summary.as_dict()
    assert "H_gamma" not in summary.as_dict()


def test_schedule_selection_cannot_depend_on_warning_fields_because_summaries_have_none():
    records = tuple(
        _record(seed, replicate, trait_loss_time=10 if replicate in {0, 1} else None)
        for seed in (1, 2)
        for replicate in range(5)
    )
    summary = _summarise_schedules(records, (1, 2), (60,), (0.30,))[0]
    selection = _select_schedule((summary,))
    assert selection["selection_uses_warning_outcomes"] is False
    assert selection["selected"] is not None
