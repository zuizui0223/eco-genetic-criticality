"""Declared primary-analysis domain from independent mutation-H1 validation.

This module freezes a selection rule *before* any mutation-conditioned H2/H3
trajectory is evaluated.  It records the full 27-cell independent validation
ledger rather than retaining only favourable cells.

A cell is primary-analysis eligible only when every independent master-seed
block achieved a same-replicate joint-support probability of at least 0.75:

    finite H1 full-state hold
    AND polymorphic high branch
    AND H-alpha/H-gamma above their baseline-warning thresholds.

The evidence comes from Actions run 28436777080, with master seeds
20260710--20260714 and 20 replicates per pair per seed.  This is a Type S
selection rule for the declared symmetric-mutation closure.  It is not a claim
that excluded cells are biologically impossible.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable

VALIDATION_RUN_ID = 28436777080
VALIDATION_MASTER_SEEDS = (20260710, 20260711, 20260712, 20260713, 20260714)
VALIDATION_REPLICATES_PER_PAIR_PER_SEED = 20
PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT = 0.75


@dataclass(frozen=True, order=True)
class MutationH1DomainCell:
    """Independent-validation evidence for one mutation/ecological cell."""

    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    pooled_joint_support: float
    seed_block_joint_support: tuple[float, float, float, float, float]

    @property
    def minimum_seed_block_support(self) -> float:
        return min(self.seed_block_joint_support)

    @property
    def primary_analysis_eligible(self) -> bool:
        return self.minimum_seed_block_support >= PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["minimum_seed_block_support"] = self.minimum_seed_block_support
        result["primary_analysis_eligible"] = self.primary_analysis_eligible
        return result


# Complete 27-cell ledger from independent_validation_v1.  Order is mutation
# rate, then A_ref, then kappa, so reports are stable and auditable.
VALIDATED_CELLS: tuple[MutationH1DomainCell, ...] = (
    MutationH1DomainCell(0.10, 0.8, 3.0, 0.57, (0.60, 0.50, 0.60, 0.50, 0.65)),
    MutationH1DomainCell(0.10, 0.8, 4.5, 0.79, (0.95, 0.75, 0.65, 0.80, 0.80)),
    MutationH1DomainCell(0.10, 0.8, 6.0, 0.85, (0.90, 0.90, 0.90, 0.80, 0.75)),
    MutationH1DomainCell(0.10, 1.0, 3.0, 0.57, (0.60, 0.50, 0.65, 0.50, 0.60)),
    MutationH1DomainCell(0.10, 1.0, 4.5, 0.68, (0.80, 0.60, 0.65, 0.70, 0.65)),
    MutationH1DomainCell(0.10, 1.0, 6.0, 0.81, (0.75, 0.90, 0.70, 0.75, 0.95)),
    MutationH1DomainCell(0.10, 1.2, 3.0, 0.34, (0.45, 0.35, 0.35, 0.30, 0.25)),
    MutationH1DomainCell(0.10, 1.2, 4.5, 0.73, (0.60, 0.70, 0.75, 0.85, 0.75)),
    MutationH1DomainCell(0.10, 1.2, 6.0, 0.72, (0.75, 0.90, 0.60, 0.75, 0.60)),
    MutationH1DomainCell(0.15, 0.8, 3.0, 0.77, (0.85, 0.70, 0.80, 0.80, 0.70)),
    MutationH1DomainCell(0.15, 0.8, 4.5, 0.89, (0.85, 1.00, 0.75, 0.90, 0.95)),
    MutationH1DomainCell(0.15, 0.8, 6.0, 0.93, (0.95, 0.90, 0.90, 0.90, 1.00)),
    MutationH1DomainCell(0.15, 1.0, 3.0, 0.70, (0.75, 0.65, 0.50, 0.90, 0.70)),
    MutationH1DomainCell(0.15, 1.0, 4.5, 0.83, (0.80, 0.85, 0.70, 0.95, 0.85)),
    MutationH1DomainCell(0.15, 1.0, 6.0, 0.91, (0.85, 0.75, 0.95, 1.00, 1.00)),
    MutationH1DomainCell(0.15, 1.2, 3.0, 0.11, (0.05, 0.10, 0.10, 0.15, 0.15)),
    MutationH1DomainCell(0.15, 1.2, 4.5, 0.86, (0.90, 0.90, 0.75, 0.80, 0.95)),
    MutationH1DomainCell(0.15, 1.2, 6.0, 0.79, (0.65, 0.95, 0.75, 0.75, 0.85)),
    MutationH1DomainCell(0.20, 0.8, 3.0, 0.87, (0.80, 0.80, 0.85, 1.00, 0.90)),
    MutationH1DomainCell(0.20, 0.8, 4.5, 0.91, (0.85, 0.95, 0.90, 1.00, 0.85)),
    MutationH1DomainCell(0.20, 0.8, 6.0, 0.95, (0.90, 0.95, 0.95, 0.95, 1.00)),
    MutationH1DomainCell(0.20, 1.0, 3.0, 0.80, (0.85, 0.80, 0.60, 0.85, 0.90)),
    MutationH1DomainCell(0.20, 1.0, 4.5, 0.95, (0.90, 1.00, 0.95, 1.00, 0.90)),
    MutationH1DomainCell(0.20, 1.0, 6.0, 0.87, (0.85, 0.75, 0.90, 0.95, 0.90)),
    MutationH1DomainCell(0.20, 1.2, 3.0, 0.00, (0.00, 0.00, 0.00, 0.00, 0.00)),
    MutationH1DomainCell(0.20, 1.2, 4.5, 0.83, (0.90, 0.85, 0.80, 0.75, 0.85)),
    MutationH1DomainCell(0.20, 1.2, 6.0, 0.94, (0.95, 1.00, 0.85, 0.95, 0.95)),
)


def primary_analysis_cells() -> tuple[MutationH1DomainCell, ...]:
    """Return eligible cells under the frozen all-seed-block rule."""
    validate_primary_domain()
    return tuple(cell for cell in VALIDATED_CELLS if cell.primary_analysis_eligible)


def excluded_cells() -> tuple[MutationH1DomainCell, ...]:
    """Return validation cells retained outside the primary H2/H3 analysis."""
    validate_primary_domain()
    return tuple(cell for cell in VALIDATED_CELLS if not cell.primary_analysis_eligible)


def find_validated_cell(
    mutation_rate: float,
    area_reference: float,
    interaction_feedback: float,
) -> MutationH1DomainCell:
    """Retrieve an evidence cell exactly; reject unvalidated configurations."""
    query = (float(mutation_rate), float(area_reference), float(interaction_feedback))
    for cell in VALIDATED_CELLS:
        if (cell.mutation_rate, cell.area_reference, cell.interaction_feedback) == query:
            return cell
    raise ValueError(
        "configuration was not part of independent_validation_v1: "
        f"mutation_rate={query[0]}, area_reference={query[1]}, interaction_feedback={query[2]}"
    )


def is_primary_analysis_cell(
    mutation_rate: float,
    area_reference: float,
    interaction_feedback: float,
) -> bool:
    """Return whether a validated cell is eligible for primary H2/H3 analysis."""
    return find_validated_cell(mutation_rate, area_reference, interaction_feedback).primary_analysis_eligible


def domain_manifest() -> dict[str, object]:
    """Serialize policy, evidence, and both cell sets for later campaign manifests."""
    validate_primary_domain()
    primary = primary_analysis_cells()
    excluded = excluded_cells()
    return {
        "source": {
            "validation_run_id": VALIDATION_RUN_ID,
            "campaign_role": "independent_validation_v1",
            "master_seeds": list(VALIDATION_MASTER_SEEDS),
            "replicates_per_pair_per_seed": VALIDATION_REPLICATES_PER_PAIR_PER_SEED,
        },
        "selection_rule": {
            "predicate": "minimum independent master-seed-block joint support >= 0.75",
            "minimum_seed_block_support": PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT,
            "joint_support": "finite H1 full-state hold AND polymorphic high branch AND H-alpha/H-gamma above warning thresholds in the same replicate",
        },
        "counts": {
            "validated_cells": len(VALIDATED_CELLS),
            "primary_analysis_cells": len(primary),
            "excluded_cells": len(excluded),
        },
        "primary_analysis_cells": [cell.as_dict() for cell in primary],
        "excluded_cells": [cell.as_dict() for cell in excluded],
        "complete_validation_ledger": [cell.as_dict() for cell in VALIDATED_CELLS],
    }


def validate_primary_domain(cells: Iterable[MutationH1DomainCell] = VALIDATED_CELLS) -> None:
    """Guard the frozen ledger against accidental edits or partial replacement."""
    values = tuple(cells)
    if len(values) != 27:
        raise ValueError("independent validation ledger must retain all 27 mutation-rate/parameter cells")
    keys = {(cell.mutation_rate, cell.area_reference, cell.interaction_feedback) for cell in values}
    if len(keys) != len(values):
        raise ValueError("independent validation ledger contains duplicate cells")
    expected_rates = {0.10, 0.15, 0.20}
    if {cell.mutation_rate for cell in values} != expected_rates:
        raise ValueError("independent validation ledger must retain all three declared mutation rates")
    if any(len(cell.seed_block_joint_support) != len(VALIDATION_MASTER_SEEDS) for cell in values):
        raise ValueError("each cell must preserve one support value for every independent master seed")
    if any(not 0.0 <= support <= 1.0 for cell in values for support in cell.seed_block_joint_support):
        raise ValueError("seed-block supports must lie in [0, 1]")
    if len(tuple(cell for cell in values if cell.minimum_seed_block_support >= PRIMARY_MINIMUM_SEED_BLOCK_SUPPORT)) != 12:
        raise ValueError("frozen policy must produce exactly 12 primary-analysis cells")
