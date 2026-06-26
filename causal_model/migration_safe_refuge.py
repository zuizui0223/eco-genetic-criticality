"""Exact lower bounds for the simulator's allele-migration step.

The declared simulator uses

    p_j^mig = (1-m) p_j^sel + m sum_i w_i p_i^sel,

where weights are current census weights. This module records the exact convex
combination lower envelope needed to relax the no-migration restriction in a
refuge certificate. It does not assert a lower source bound; callers must
supply one or prove it for the chosen metapopulation region.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MigrationLowerBound:
    migration_rate: float
    focal_selected_lower_bound: float
    metapopulation_selected_mean_lower_bound: float
    post_migration_lower_bound: float
    retains_target: bool


def migration_safe_allele_lower_bound(
    migration_rate: float,
    focal_selected_lower_bound: float,
    metapopulation_selected_mean_lower_bound: float,
    target_lower_bound: float,
) -> MigrationLowerBound:
    """Return the exact convex-combination lower bound after migration.

    If a focal patch has selected frequency at least ``a`` and the census-
    weighted selected metapopulation mean is at least ``b``, then the simulator
    update guarantees

        p_focal^mig >= (1-m) a + m b.

    Thus migration preserves a focal allele threshold iff this lower envelope is
    at least that target. In particular, migration cannot reduce a common lower
    bound shared by every source patch.
    """
    for value, name in (
        (migration_rate, "migration_rate"),
        (focal_selected_lower_bound, "focal_selected_lower_bound"),
        (metapopulation_selected_mean_lower_bound, "metapopulation_selected_mean_lower_bound"),
        (target_lower_bound, "target_lower_bound"),
    ):
        if not 0.0 <= value <= 1.0:
            raise ValueError(f"{name} must lie in [0, 1]")
    lower = (1.0 - migration_rate) * focal_selected_lower_bound + migration_rate * metapopulation_selected_mean_lower_bound
    return MigrationLowerBound(
        migration_rate=float(migration_rate),
        focal_selected_lower_bound=float(focal_selected_lower_bound),
        metapopulation_selected_mean_lower_bound=float(metapopulation_selected_mean_lower_bound),
        post_migration_lower_bound=lower,
        retains_target=lower >= target_lower_bound,
    )
