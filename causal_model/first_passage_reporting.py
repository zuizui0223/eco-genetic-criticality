"""Censoring-aware summaries for ordered first-passage events.

A stochastic lead claim such as ``tau_warning < tau_trait`` has two distinct
frequencies:

1. the conditional lead probability among replicates where both event times are
   observed; and
2. the unconditional observed-lead fraction among all replicates.

They are equal only when every replicate supplies both events.  This module
keeps the denominator visible rather than silently treating censored runs as
no-lead outcomes or as observed terminal-time events.
"""
from __future__ import annotations

from dataclasses import dataclass
from statistics import median
from typing import Sequence


@dataclass(frozen=True)
class FirstPassageComparison:
    """Censoring-aware comparison of a warning time against a reference time.

    ``time_differences`` stores ``warning_time - reference_time`` only for
    valid pairs.  A negative difference is therefore a warning lead.
    """

    replicate_count: int
    valid_pair_count: int
    censored_pair_count: int
    valid_pair_probability: float
    lead_count: int
    conditional_lead_probability: float | None
    unconditional_observed_lead_fraction: float
    time_differences: tuple[int, ...]
    median_time_difference: float | None


def compare_first_passage_times(
    warning_times: Sequence[int | None],
    reference_times: Sequence[int | None],
) -> FirstPassageComparison:
    """Compare predeclared event times without imputing censored observations.

    The conditional probability is ``None`` when no replicate observes both
    events.  Returning zero in that case would falsely claim evidence of no
    lead rather than absence of an estimable comparison.
    """
    if len(warning_times) != len(reference_times):
        raise ValueError("warning_times and reference_times must have equal length")
    if not warning_times:
        raise ValueError("at least one replicate is required")

    differences = tuple(
        warning - reference
        for warning, reference in zip(warning_times, reference_times)
        if warning is not None and reference is not None
    )
    replicate_count = len(warning_times)
    valid_pair_count = len(differences)
    lead_count = sum(difference < 0 for difference in differences)
    conditional = None if not differences else lead_count / valid_pair_count
    return FirstPassageComparison(
        replicate_count=replicate_count,
        valid_pair_count=valid_pair_count,
        censored_pair_count=replicate_count - valid_pair_count,
        valid_pair_probability=valid_pair_count / replicate_count,
        lead_count=lead_count,
        conditional_lead_probability=conditional,
        unconditional_observed_lead_fraction=lead_count / replicate_count,
        time_differences=differences,
        median_time_difference=None if not differences else median(differences),
    )


def comparison_as_dict(comparison: FirstPassageComparison) -> dict[str, object]:
    """Return a stable serialisable representation for experiment artifacts."""
    return {
        "replicate_count": comparison.replicate_count,
        "valid_pair_count": comparison.valid_pair_count,
        "censored_pair_count": comparison.censored_pair_count,
        "valid_pair_probability": comparison.valid_pair_probability,
        "lead_count": comparison.lead_count,
        "conditional_lead_probability": comparison.conditional_lead_probability,
        "unconditional_observed_lead_fraction": comparison.unconditional_observed_lead_fraction,
        "time_differences": list(comparison.time_differences),
        "median_time_difference": comparison.median_time_difference,
    }
