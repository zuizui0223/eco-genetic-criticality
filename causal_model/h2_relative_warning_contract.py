"""Predeclared contract for the conditional relative-warning H2 proposition.

H2-A (the existing proposition) compares fixed absolute H-alpha/H-gamma warning
thresholds with realised trait loss.  It is retained unchanged.  The
mutation-conditioned primary chain showed that, under a stationary mutation
closure and a finite 30-generation horizon, most records were right-censored:
there were too few same-replicate warning/loss pairs to assign H2-A a truth
value for that closure.

H2-R is therefore a *separate* dynamic hypothesis, not a re-labelling of H2-A:

    given a predeclared monotone decline in interaction support,
    tau_delta_H_x(r) < tau_trait_realised,

where x is H-alpha or H-gamma and tau_delta_H_x(r) is the first post-baseline
observation at which H_x has declined by at least r relative to its own
baseline.  A finite run may support, contradict, or leave this proposition
censored; no outcome is silently converted to another.

This module defines the measurement and calibration contract only.  It does not
alter H1, H3, the symmetric-mutation closure, or the existing fixed-threshold
H2-A implementation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Iterable, Sequence


H2_ABSOLUTE_ORIGINAL_ID = "H2-A"
H2_RELATIVE_WARNING_ID = "H2-R"
STATIONARY_PRIMARY_CHAIN_RUN_ID = 28456092898

# These values are a family, not a post-hoc chosen single threshold. Every value
# must be reported in a later H2-R validation artifact.
DEFAULT_RELATIVE_DECLINE_FRACTIONS = (0.05, 0.10, 0.20)
DEFAULT_CALIBRATION_HORIZONS = (60, 120)
DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES = (0.15, 0.30, 0.45)
TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL = (0.30, 0.70)


@dataclass(frozen=True)
class RelativeWarningDefinition:
    """A within-trajectory, post-baseline relative-diversity warning rule."""

    diversity_id: str
    relative_decline_fraction: float
    require_positive_baseline: bool = True

    def __post_init__(self) -> None:
        if self.diversity_id not in {"H_alpha", "H_gamma"}:
            raise ValueError("diversity_id must be H_alpha or H_gamma")
        if not 0.0 < self.relative_decline_fraction < 1.0:
            raise ValueError("relative_decline_fraction must lie in (0, 1)")

    def threshold(self, baseline: float) -> float | None:
        """Return the fixed within-run threshold, or None for an ineligible baseline."""
        _validate_heterozygosity(baseline)
        if self.require_positive_baseline and baseline <= 0.0:
            return None
        return (1.0 - self.relative_decline_fraction) * baseline


@dataclass(frozen=True)
class RelativeWarningComparison:
    """One censored first-passage comparison for H2-R."""

    definition: RelativeWarningDefinition
    baseline: float
    warning_threshold: float | None
    warning_time: int | None
    trait_loss_time: int | None
    baseline_eligible: bool
    valid_pair: bool
    censored: bool
    warning_leads: bool | None
    lead_time_trait_minus_warning: int | None

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["definition"] = asdict(self.definition)
        return result


@dataclass(frozen=True)
class DeteriorationCalibrationCandidate:
    """A candidate schedule judged only by trait-loss event availability."""

    horizon: int
    total_normalized_barrier_increase: float
    trait_loss_probability_by_seed_block: tuple[float, ...]

    def __post_init__(self) -> None:
        if self.horizon < 1:
            raise ValueError("horizon must be positive")
        if self.total_normalized_barrier_increase <= 0.0:
            raise ValueError("total_normalized_barrier_increase must be positive")
        if not self.trait_loss_probability_by_seed_block:
            raise ValueError("at least one seed-block probability is required")
        if any(not 0.0 <= value <= 1.0 for value in self.trait_loss_probability_by_seed_block):
            raise ValueError("trait-loss probabilities must lie in [0, 1]")

    @property
    def pooled_trait_loss_probability(self) -> float:
        return sum(self.trait_loss_probability_by_seed_block) / len(self.trait_loss_probability_by_seed_block)

    @property
    def every_seed_block_in_target_interval(self) -> bool:
        lower, upper = TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL
        return all(lower <= value <= upper for value in self.trait_loss_probability_by_seed_block)


@dataclass(frozen=True)
class H2RCalibrationSelection:
    """Outcome of a trait-loss-only schedule calibration for one domain cell."""

    selected: DeteriorationCalibrationCandidate | None
    candidates: tuple[DeteriorationCalibrationCandidate, ...]
    reason: str

    def as_dict(self) -> dict[str, object]:
        return {
            "selected": None if self.selected is None else asdict(self.selected),
            "candidates": [asdict(candidate) for candidate in self.candidates],
            "reason": self.reason,
            "selection_uses_warning_outcomes": False,
        }


def first_relative_decline_time(
    values: Sequence[float],
    definition: RelativeWarningDefinition,
) -> tuple[float | None, int | None]:
    """Return baseline-derived threshold and first *post-baseline* crossing.

    Generation zero is deliberately excluded: a low or zero baseline is not an
    early warning.  A missing crossing is returned as ``None`` for censoring.
    """
    observations = tuple(float(value) for value in values)
    if not observations:
        raise ValueError("values must be nonempty")
    for value in observations:
        _validate_heterozygosity(value)
    threshold = definition.threshold(observations[0])
    if threshold is None:
        return None, None
    for time, value in enumerate(observations[1:], start=1):
        if value <= threshold:
            return threshold, time
    return threshold, None


def compare_relative_warning(
    values: Sequence[float],
    *,
    trait_loss_time: int | None,
    definition: RelativeWarningDefinition,
) -> RelativeWarningComparison:
    """Keep baseline ineligibility and event censoring distinct for H2-R."""
    observations = tuple(float(value) for value in values)
    if not observations:
        raise ValueError("values must be nonempty")
    threshold, warning_time = first_relative_decline_time(observations, definition)
    baseline_eligible = threshold is not None
    valid_pair = baseline_eligible and warning_time is not None and trait_loss_time is not None
    lead = None if not valid_pair else warning_time < trait_loss_time
    return RelativeWarningComparison(
        definition=definition,
        baseline=observations[0],
        warning_threshold=threshold,
        warning_time=warning_time,
        trait_loss_time=trait_loss_time,
        baseline_eligible=baseline_eligible,
        valid_pair=valid_pair,
        censored=not valid_pair,
        warning_leads=lead,
        lead_time_trait_minus_warning=None if not valid_pair else trait_loss_time - warning_time,
    )


def select_trait_loss_only_calibration(
    candidates: Iterable[DeteriorationCalibrationCandidate],
) -> H2RCalibrationSelection:
    """Select a schedule without inspecting H-alpha/H-gamma warning outcomes.

    A candidate is eligible only if *every* calibration seed block has a trait
    loss probability in the declared 0.30--0.70 interval.  Among eligible
    candidates, choose the one nearest the midpoint (0.50), then the shorter
    horizon, then the smaller normalized barrier increase.  If none qualify,
    return an explicit no-selection result instead of widening the rule after
    looking at warnings.
    """
    values = tuple(candidates)
    if not values:
        raise ValueError("candidates must be nonempty")
    eligible = tuple(candidate for candidate in values if candidate.every_seed_block_in_target_interval)
    if not eligible:
        return H2RCalibrationSelection(None, values, "no candidate met the all-seed-block trait-loss availability target")
    selected = min(
        eligible,
        key=lambda candidate: (
            abs(candidate.pooled_trait_loss_probability - 0.50),
            candidate.horizon,
            candidate.total_normalized_barrier_increase,
        ),
    )
    return H2RCalibrationSelection(selected, values, "selected using trait-loss availability only; warning outcomes were not inspected")


def h2r_protocol_manifest() -> dict[str, object]:
    """Serialize the revised proposition and its anti-post-hoc safeguards."""
    return {
        "original_proposition": {
            "id": H2_ABSOLUTE_ORIGINAL_ID,
            "statement": "fixed absolute H-alpha/H-gamma warning thresholds precede realised high-trait loss under declared conditions",
            "status_after_stationary_mutation_primary_chain": "unresolved_due_to_right_censoring_not_false",
            "source_run_id": STATIONARY_PRIMARY_CHAIN_RUN_ID,
        },
        "revised_proposition": {
            "id": H2_RELATIVE_WARNING_ID,
            "statement": "under a predeclared monotone interaction-support deterioration schedule, a post-baseline relative decline in H-alpha or H-gamma precedes realised high-trait loss conditionally on an eligible observed event pair",
            "relative_warning_definition": "H_x(t) <= (1-r) H_x(0), evaluated only for t > 0",
            "relative_decline_fractions": list(DEFAULT_RELATIVE_DECLINE_FRACTIONS),
            "evidence_label": "H: dynamic hypothesis; any finite result is Type S",
        },
        "calibration": {
            "candidate_horizons": list(DEFAULT_CALIBRATION_HORIZONS),
            "candidate_total_normalized_barrier_increases": list(DEFAULT_TOTAL_NORMALIZED_BARRIER_INCREASES),
            "target_trait_loss_probability_interval_per_seed_block": list(TARGET_TRAIT_LOSS_PROBABILITY_INTERVAL),
            "selection_uses": "trait-loss event availability only",
            "selection_forbidden_inputs": ["H-alpha warning lead", "H-gamma warning lead", "lead-time magnitude"],
        },
        "invariants": {
            "h1_h3_modified": False,
            "symmetric_mutation_closure_modified": False,
            "absolute_h2_implementation_modified": False,
            "baseline_crossing_is_not_warning": True,
            "missing_warning_or_trait_loss_is_censored": True,
        },
    }


def _validate_heterozygosity(value: float) -> None:
    if not 0.0 <= float(value) <= 1.0:
        raise ValueError("heterozygosity values must lie in [0, 1]")
