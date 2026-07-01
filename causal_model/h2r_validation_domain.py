"""Frozen domain for independent H2-R relative-warning validation.

The ramp-and-hold trait-loss-only calibration (Actions run 28496735824) examined
all 12 mutation-H1 primary cells without calculating any genetic warning value.
Exactly one configuration met the predeclared all-seed-block trait-loss
availability rule.  This module records that selection so the next validation
cannot search across cells or schedules after seeing H-alpha/H-gamma outcomes.

This is scope control, not H2-R evidence.  The selected configuration is the
only domain of the upcoming independent validation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

from causal_model.h2r_ramp_hold_trait_loss_calibration import RampHoldSchedule

RAMP_HOLD_CALIBRATION_RUN_ID = 28496735824
RAMP_HOLD_CALIBRATION_MASTER_SEEDS = (20261010, 20261011, 20261012, 20261013, 20261014)
RAMP_HOLD_CALIBRATION_REPLICATES_PER_CELL_SEED = 5
TRAIT_LOSS_TARGET_INTERVAL = (0.30, 0.70)


@dataclass(frozen=True)
class H2RValidationDomain:
    mutation_rate: float
    area_reference: float
    interaction_feedback: float
    schedule: RampHoldSchedule
    calibration_pooled_trait_loss_probability: float
    calibration_seed_block_trait_loss_probabilities: tuple[float, ...]
    calibration_seed_block_denominators: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.mutation_rate <= 0.0:
            raise ValueError("mutation_rate must be positive")
        if self.area_reference <= 0.0 or self.interaction_feedback <= 0.0:
            raise ValueError("area_reference and interaction_feedback must be positive")
        if len(self.calibration_seed_block_trait_loss_probabilities) != 5:
            raise ValueError("five calibration seed-block probabilities are required")
        if len(self.calibration_seed_block_denominators) != 5:
            raise ValueError("five calibration seed-block denominators are required")
        lower, upper = TRAIT_LOSS_TARGET_INTERVAL
        if any(not lower <= value <= upper for value in self.calibration_seed_block_trait_loss_probabilities):
            raise ValueError("each selected calibration seed-block probability must lie in the target interval")
        if any(value < 1 for value in self.calibration_seed_block_denominators):
            raise ValueError("selected calibration seed-block denominators must be positive")

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["schedule"] = self.schedule.as_dict()
        return result


SELECTED_VALIDATION_DOMAIN = H2RValidationDomain(
    mutation_rate=0.10,
    area_reference=0.8,
    interaction_feedback=6.0,
    schedule=RampHoldSchedule(
        ramp_generations=30,
        hold_generations=90,
        total_normalized_barrier_increase=0.15,
    ),
    calibration_pooled_trait_loss_probability=0.45,
    calibration_seed_block_trait_loss_probabilities=(0.50, 0.40, 0.40, 0.50, 0.50),
    calibration_seed_block_denominators=(2, 5, 5, 4, 4),
)


def h2r_validation_domain_manifest() -> dict[str, object]:
    return {
        "domain_selection_campaign": "h2r_ramp_hold_trait_loss_only_calibration_v2",
        "domain_selection_run_id": RAMP_HOLD_CALIBRATION_RUN_ID,
        "domain_selection_master_seeds": list(RAMP_HOLD_CALIBRATION_MASTER_SEEDS),
        "domain_selection_replicates_per_cell_seed": RAMP_HOLD_CALIBRATION_REPLICATES_PER_CELL_SEED,
        "selection_target_trait_loss_probability_interval": list(TRAIT_LOSS_TARGET_INTERVAL),
        "selected_domain_count": 1,
        "selected_domain": SELECTED_VALIDATION_DOMAIN.as_dict(),
        "selection_used_warning_outcomes": False,
        "excluded_primary_cell_count": 11,
        "validation_scope": "the selected configuration only; do not generalize to the 11 calibration-unselected primary cells",
    }
