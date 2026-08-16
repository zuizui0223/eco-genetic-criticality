from __future__ import annotations

from pathlib import Path

import pytest

from causal_model.h3_fragmentation_gradient_sensitivity import (
    CANONICAL_SCIENTIFIC_COMMIT,
    FRAGMENT_PATCH_COUNTS,
    FRESH_GRADIENT_MASTER_SEEDS,
    _outcome_seed,
    _validate_patch_counts,
    equal_isolated_scenario,
    theoretical_fragment_k,
)


ROOT = Path(__file__).resolve().parents[1]


def test_preregistered_gradient_constants_are_frozen() -> None:
    assert FRAGMENT_PATCH_COUNTS == (1, 2, 3, 4, 6, 8, 12, 16)
    assert FRESH_GRADIENT_MASTER_SEEDS == (20260820, 20260821, 20260822, 20260823, 20260824)
    assert CANONICAL_SCIENTIFIC_COMMIT == "dd8ee379d0d3518194c767d16402042525bc00dc"


def test_equal_isolated_gradient_preserves_total_area_and_zero_migration() -> None:
    for patch_count in FRAGMENT_PATCH_COUNTS:
        scenario = equal_isolated_scenario(4.0, patch_count)
        assert scenario.scenario_id == f"equal_isolated_n{patch_count}"
        assert len(scenario.patch_areas) == patch_count
        assert sum(scenario.patch_areas) == pytest.approx(4.0)
        assert len(set(scenario.patch_areas)) == 1
        assert scenario.migration_rate == 0.0


def test_theoretical_k_tracks_patch_area_without_becoming_a_fitted_threshold() -> None:
    assert theoretical_fragment_k(interaction_feedback=6.0, patch_area=1.0, area_reference=1.0) == pytest.approx(6.0)
    assert theoretical_fragment_k(interaction_feedback=6.0, patch_area=0.5, area_reference=1.0) == pytest.approx(3.0)
    with pytest.raises(ValueError, match="area_reference"):
        theoretical_fragment_k(interaction_feedback=6.0, patch_area=1.0, area_reference=0.0)


def test_patch_count_validation_requires_ordered_one_patch_reference() -> None:
    assert _validate_patch_counts(FRAGMENT_PATCH_COUNTS) == FRAGMENT_PATCH_COUNTS
    with pytest.raises(ValueError, match="begin"):
        _validate_patch_counts((2, 4, 8))
    with pytest.raises(ValueError, match="ascending"):
        _validate_patch_counts((1, 4, 2))
    with pytest.raises(ValueError, match="unique"):
        _validate_patch_counts((1, 2, 2))


def test_outcome_seed_is_deterministic_and_patch_count_specific() -> None:
    seed = 123456
    values = [_outcome_seed(seed, patch_count) for patch_count in FRAGMENT_PATCH_COUNTS]
    assert values == [_outcome_seed(seed, patch_count) for patch_count in FRAGMENT_PATCH_COUNTS]
    assert len(values) == len(set(values))
    assert all(0 <= value < 2**31 - 1 for value in values)


def test_protocol_declares_new_sensitivity_without_reopening_parent_ledger() -> None:
    text = (ROOT / "docs/H3_FRAGMENTATION_GRADIENT_SENSITIVITY_PROTOCOL.md").read_text(encoding="utf-8")
    assert "declared before inspecting any new gradient outcome" in text
    assert "does not alter the closed H1/H3 evidence ledger" in text
    assert "20260820, 20260821, 20260822, 20260823, 20260824" in text
    assert "1, 2, 3, 4, 6, 8, 12, 16" in text
    assert "No warning endpoint is selected, tuned, or evaluated" in text
    assert "repeated measures" in text
