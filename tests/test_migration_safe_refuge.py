import pytest

from causal_model.migration_safe_refuge import migration_safe_allele_lower_bound


def test_common_source_lower_bound_is_preserved_by_migration() -> None:
    bound = migration_safe_allele_lower_bound(0.7, 0.4, 0.4, 0.4)
    assert bound.post_migration_lower_bound == pytest.approx(0.4)
    assert bound.retains_target


def test_migration_can_erode_a_focal_refuge_when_sources_are_lower() -> None:
    bound = migration_safe_allele_lower_bound(0.5, 0.8, 0.2, 0.6)
    assert bound.post_migration_lower_bound == pytest.approx(0.5)
    assert not bound.retains_target


def test_migration_can_rescue_a_low_focal_patch() -> None:
    bound = migration_safe_allele_lower_bound(0.5, 0.2, 0.8, 0.4)
    assert bound.post_migration_lower_bound == pytest.approx(0.5)
    assert bound.retains_target


def test_invalid_probability_is_rejected() -> None:
    with pytest.raises(ValueError):
        migration_safe_allele_lower_bound(1.1, 0.2, 0.2, 0.2)
