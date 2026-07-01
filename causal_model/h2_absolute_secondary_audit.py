"""Secondary audit of the pre-existing fixed-threshold H2-A definition.

This module performs no new simulation.  It reads the raw trajectory series from
an H2-R independent-validation JSON artifact and applies the original H2-A
thresholds without changing the threshold, cell, schedule, or source seeds.

H2-A is a different proposition from H2-R:

    H2-A: first H_alpha <= 0.20 or H_gamma <= 0.20 precedes realised trait loss.
    H2-R: a post-baseline relative decline from each trajectory's own baseline
          precedes realised trait loss.

The audit deliberately reports all observed lead, tie, lag, and censored
records.  It cannot convert mixed ordering into a universal refutation, but
observed lags prevent the fixed threshold from being retained as a robust
canonical early-warning rule in the audited finite closure.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

H2A_FIXED_THRESHOLD = 0.20
H2R_INDEPENDENT_VALIDATION_RUN_ID = 28500796310
H2R_INDEPENDENT_VALIDATION_ARTIFACT_DIGEST = "sha256:7d2bbed84ddf57486896c0ca231fd82f2e0915699e391c155d288f5c9db8a6ff"


@dataclass(frozen=True)
class AbsoluteThresholdDefinition:
    diversity_id: str
    threshold: float = H2A_FIXED_THRESHOLD

    def __post_init__(self) -> None:
        if self.diversity_id not in {"H_alpha", "H_gamma"}:
            raise ValueError("diversity_id must be H_alpha or H_gamma")
        if not 0.0 <= self.threshold <= 1.0:
            raise ValueError("threshold must lie in [0, 1]")


@dataclass(frozen=True)
class AbsoluteThresholdComparison:
    definition: AbsoluteThresholdDefinition
    warning_time: int | None
    trait_loss_time: int | None
    valid_pair: bool
    censored: bool
    warning_leads: bool | None
    warning_ties: bool | None
    warning_lags: bool | None
    lead_time_trait_minus_warning: int | None

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["definition"] = asdict(self.definition)
        return result


def first_absolute_threshold_time(values: Sequence[float], threshold: float) -> int | None:
    """Return the first finite-model time at or below the *pre-existing* threshold."""
    if not 0.0 <= threshold <= 1.0:
        raise ValueError("threshold must lie in [0, 1]")
    if not values:
        raise ValueError("values must be nonempty")
    for time, value in enumerate(values):
        numeric = float(value)
        if not 0.0 <= numeric <= 1.0:
            raise ValueError("diversity values must lie in [0, 1]")
        if numeric <= threshold:
            return time
    return None


def compare_fixed_absolute_warning(
    values: Sequence[float],
    *,
    trait_loss_time: int | None,
    definition: AbsoluteThresholdDefinition,
) -> AbsoluteThresholdComparison:
    """Compare a fixed H2-A crossing with realised trait loss without dropping censoring."""
    warning_time = first_absolute_threshold_time(values, definition.threshold)
    valid = warning_time is not None and trait_loss_time is not None
    lead_time = None if not valid else trait_loss_time - warning_time
    return AbsoluteThresholdComparison(
        definition=definition,
        warning_time=warning_time,
        trait_loss_time=trait_loss_time,
        valid_pair=valid,
        censored=not valid,
        warning_leads=None if not valid else warning_time < trait_loss_time,
        warning_ties=None if not valid else warning_time == trait_loss_time,
        warning_lags=None if not valid else warning_time > trait_loss_time,
        lead_time_trait_minus_warning=lead_time,
    )


def audit_h2a_from_h2r_validation_payload(payload: Mapping[str, Any]) -> dict[str, object]:
    """Reanalyse an H2-R JSON artifact at H2-A's fixed 0.20 thresholds.

    The required payload shape is the JSON output from
    ``h2r_independent_relative_validation``.  Unavailable trajectories are kept
    in the global denominator but naturally cannot contribute an event ordering.
    """
    records = payload.get("records")
    if not isinstance(records, list) or not records:
        raise ValueError("payload must contain a nonempty records list")
    definitions = (
        AbsoluteThresholdDefinition("H_alpha"),
        AbsoluteThresholdDefinition("H_gamma"),
    )
    comparison_rows: list[dict[str, object]] = []
    unavailable_count = 0
    for record in records:
        if not isinstance(record, Mapping):
            raise ValueError("every record must be a mapping")
        outcome = record.get("outcome")
        if outcome is None:
            unavailable_count += 1
            continue
        if not isinstance(outcome, Mapping):
            raise ValueError("record outcome must be a mapping or null")
        trait_time = outcome.get("trait_loss_time_post_baseline")
        if trait_time is not None:
            trait_time = int(trait_time)
        for definition in definitions:
            series_key = "h_alpha_series" if definition.diversity_id == "H_alpha" else "h_gamma_series"
            values = outcome.get(series_key)
            if not isinstance(values, list):
                raise ValueError(f"available trajectory lacks {series_key}")
            comparison = compare_fixed_absolute_warning(
                values,
                trait_loss_time=trait_time,
                definition=definition,
            )
            comparison_rows.append(
                {
                    "master_seed": int(record["master_seed"]),
                    "replicate_index": int(record["replicate_index"]),
                    "calibration_seed": int(record["calibration_seed"]),
                    "trajectory_available": True,
                    **comparison.as_dict(),
                }
            )
    summaries = [_summarise_definition(comparison_rows, definition) for definition in definitions]
    return {
        "audit": {
            "kind": "secondary_fixed_threshold_reanalysis",
            "source_campaign": "h2r_independent_relative_warning_validation_v1",
            "source_run_id": H2R_INDEPENDENT_VALIDATION_RUN_ID,
            "source_artifact_digest": H2R_INDEPENDENT_VALIDATION_ARTIFACT_DIGEST,
            "new_simulation": False,
            "thresholds_preexisted_source_run": True,
            "threshold_selection_after_observation": False,
            "cell_or_schedule_selection_repeated": False,
            "scope": "selected H2-R validation domain only",
            "finite_evidence_label": "Type S",
        },
        "denominators": {
            "attempted_source_records": len(records),
            "trajectory_available_count": len(records) - unavailable_count,
            "trajectory_unavailable_count": unavailable_count,
        },
        "definitions": [asdict(definition) for definition in definitions],
        "endpoint_summaries": summaries,
        "comparison_rows": comparison_rows,
        "canonical_interpretation": {
            "h2a_global_truth_value": "not_assigned",
            "h2a_selected_domain_status": "mixed_ordering_not_retained_as_robust_absolute_warning_rule",
            "reason": "each fixed-threshold metric has observed trait-loss-before-warning lags among valid pairs",
            "h2r_status_unchanged": "separate relative-warning proposition; not altered by this audit",
        },
    }


def write_h2a_secondary_audit_artifacts(
    audit: Mapping[str, Any],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write all pair-level audit records and a self-contained JSON ledger."""
    csv_target, json_target = Path(csv_path), Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = list(audit["comparison_rows"])
    if not rows:
        raise ValueError("audit contains no available trajectory comparison rows")
    flattened = [_flatten_comparison_row(row) for row in rows]
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted({key for row in flattened for key in row}))
        writer.writeheader()
        writer.writerows(flattened)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)


def _summarise_definition(
    rows: Iterable[Mapping[str, object]],
    definition: AbsoluteThresholdDefinition,
) -> dict[str, object]:
    subset = [row for row in rows if row["definition"] == asdict(definition)]
    if not subset:
        raise RuntimeError("missing expected absolute-threshold definition")
    valid = [row for row in subset if row["valid_pair"] is True]
    lead = [row for row in valid if row["warning_leads"] is True]
    tie = [row for row in valid if row["warning_ties"] is True]
    lag = [row for row in valid if row["warning_lags"] is True]
    seed_rows = []
    for seed in sorted({int(row["master_seed"]) for row in subset}):
        by_seed = [row for row in subset if int(row["master_seed"]) == seed]
        valid_seed = [row for row in by_seed if row["valid_pair"] is True]
        seed_rows.append(
            {
                "master_seed": seed,
                "trajectory_available_count": len(by_seed),
                "warning_observed_count": sum(row["warning_time"] is not None for row in by_seed),
                "trait_loss_observed_count": sum(row["trait_loss_time"] is not None for row in by_seed),
                "valid_pair_count": len(valid_seed),
                "warning_lead_count": sum(row["warning_leads"] is True for row in valid_seed),
                "warning_tie_count": sum(row["warning_ties"] is True for row in valid_seed),
                "warning_lag_count": sum(row["warning_lags"] is True for row in valid_seed),
            }
        )
    return {
        "definition": asdict(definition),
        "trajectory_available_count": len(subset),
        "warning_observed_count": sum(row["warning_time"] is not None for row in subset),
        "trait_loss_observed_count": sum(row["trait_loss_time"] is not None for row in subset),
        "valid_pair_count": len(valid),
        "censored_count": sum(row["censored"] is True for row in subset),
        "warning_lead_count": len(lead),
        "warning_tie_count": len(tie),
        "warning_lag_count": len(lag),
        "warning_lead_probability_among_valid_pairs": None if not valid else len(lead) / len(valid),
        "seed_blocks": seed_rows,
    }


def _flatten_comparison_row(row: Mapping[str, object]) -> dict[str, object]:
    result = dict(row)
    definition = result.pop("definition")
    assert isinstance(definition, Mapping)
    result["diversity_id"] = definition["diversity_id"]
    result["threshold"] = definition["threshold"]
    return result
