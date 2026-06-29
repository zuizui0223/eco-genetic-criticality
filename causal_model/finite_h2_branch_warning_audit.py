"""Branch-conditioned finite-model audit for H2 genetic warning order.

The canonical H2 certificate proves a first-passage ordering only for a
constant-effective-size deterministic expectation map.  The finite coupled
simulator has stochastic drift, changing local effective size, trait feedback,
and potentially different long-run H1 branches.  This module tests H2 only
within finite replicates for which the existing finite H1 branch-separation
audit already supports the stronger finite H1 mechanism predicate.

For those preconditioned low-start and high-start paths separately, it compares
H-alpha, H-gamma, and allele-loss first passage against realised high-trait
loss.  Every comparison records valid-pair denominators and censoring; absent
events are never silently converted into a terminal-time observation.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.finite_h1_branch_separation_audit import (
    FiniteH1BranchPair,
    FiniteH1BranchSeparationCell,
    run_finite_h1_branch_separation_audit,
)
from causal_model.multipatch_criticality_dynamics import FirstPassageEvent
from causal_model.multipatch_criticality_experiments import ExperimentSpec, LandscapeScenario, ReplicateSummary

_WARNING_ATTRIBUTES = (
    ("H_alpha", "tau_H_alpha"),
    ("H_gamma", "tau_H_gamma"),
    ("allele_loss", "tau_allele_loss"),
)
_BRANCH_IDS = ("low_start", "high_start")


@dataclass(frozen=True)
class BranchWarningComparison:
    """One warning-versus-realised-trait first-passage comparison."""

    branch_id: str
    warning_id: str
    warning_time: int | None
    trait_loss_time: int | None
    warning_threshold: float | int | None
    warning_aggregation_rule: str | None
    trait_threshold: float | int | None
    trait_aggregation_rule: str | None
    valid_pair: bool
    censored: bool
    warning_leads: bool | None
    lead_time_trait_minus_warning: int | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class FiniteH2BranchWarningReplicate:
    """H2 results conditioned on the H1 finite branch-mechanism predicate."""

    replicate_index: int
    seed: int
    finite_h1_mechanism_supported: bool | None
    low_start: tuple[BranchWarningComparison, ...] | None
    high_start: tuple[BranchWarningComparison, ...] | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "finite_h1_mechanism_supported": self.finite_h1_mechanism_supported,
            "low_start": None if self.low_start is None else [comparison.as_dict() for comparison in self.low_start],
            "high_start": None if self.high_start is None else [comparison.as_dict() for comparison in self.high_start],
        }


@dataclass(frozen=True)
class FiniteH2BranchWarningCell:
    """One H1 finite-branch cell with low/high branch-conditioned H2 results."""

    experiment_id: str
    profile: str
    scenario_id: str
    parameters: object
    h1_precondition: Mapping[str, object]
    replicates: tuple[FiniteH2BranchWarningReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": asdict(self.parameters),
            "h1_precondition": dict(self.h1_precondition),
            "replicate_count": len(self.replicates),
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "replicate_count": len(self.replicates),
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping({"h1_precondition": self.h1_precondition, **self.summary}))
        return row


def run_finite_h2_branch_warning_audit(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    interaction_separation_threshold: float = 0.05,
    terminal_window: int = 5,
) -> tuple[FiniteH2BranchWarningCell, ...]:
    """Evaluate finite H2 first-passage ordering within finite H1 branch pairs.

    The existing finite H1 audit is run first.  Only replicate pairs marked
    `finite_h1_mechanism_supported=True` contribute branch-specific H2 warning
    comparisons.  This precondition requires both finite interaction branch
    separation and a potential high-trait switch.  Other replicates remain in
    the artifact, but are excluded from H2 denominators rather than recoded as
    warning failures.
    """
    h1_cells = run_finite_h1_branch_separation_audit(
        spec,
        scenarios=scenarios,
        interaction_separation_threshold=interaction_separation_threshold,
        terminal_window=terminal_window,
    )
    return tuple(_cell_from_h1(cell) for cell in h1_cells)


def write_finite_h2_branch_warning_artifacts(
    cells: Iterable[FiniteH2BranchWarningCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat branch-conditioned H2 summaries and complete records."""
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    fieldnames = sorted({key for row in rows for key in row})
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)


def _cell_from_h1(h1: FiniteH1BranchSeparationCell) -> FiniteH2BranchWarningCell:
    replicates = tuple(_replicate_from_h1(pair) for pair in h1.replicates)
    precondition = {
        "canonical_branch_dependent_high_trait_mode": h1.canonical_h1.branch_dependent_high_trait_mode,
        "finite_pair_available_probability": h1.summary["finite_pair_available_probability"],
        "finite_h1_mechanism_supported_probability": h1.summary["finite_h1_mechanism_supported_probability"],
        "interaction_separation_threshold": h1.interaction_separation_threshold,
        "terminal_window": h1.terminal_window,
    }
    return FiniteH2BranchWarningCell(
        experiment_id=h1.experiment_id,
        profile=h1.profile,
        scenario_id=h1.scenario_id,
        parameters=h1.parameters,
        h1_precondition=precondition,
        replicates=replicates,
        summary=_summarise_cell(replicates),
    )


def _replicate_from_h1(pair: FiniteH1BranchPair) -> FiniteH2BranchWarningReplicate:
    if pair.finite_h1_mechanism_supported is not True:
        return FiniteH2BranchWarningReplicate(
            replicate_index=pair.replicate_index,
            seed=pair.seed,
            finite_h1_mechanism_supported=pair.finite_h1_mechanism_supported,
            low_start=None,
            high_start=None,
        )
    if pair.low_start is None or pair.high_start is None:
        raise RuntimeError("finite H1 mechanism support requires both branch summaries")
    return FiniteH2BranchWarningReplicate(
        replicate_index=pair.replicate_index,
        seed=pair.seed,
        finite_h1_mechanism_supported=True,
        low_start=_comparisons_for_branch("low_start", pair.low_start),
        high_start=_comparisons_for_branch("high_start", pair.high_start),
    )


def _comparisons_for_branch(branch_id: str, summary: ReplicateSummary) -> tuple[BranchWarningComparison, ...]:
    trait_event = _event_by_name(summary.events, "tau_trait_realised")
    trait_time = summary.tau_trait_realised
    return tuple(
        _comparison(
            branch_id,
            warning_id,
            getattr(summary, attribute),
            _event_by_name(summary.events, attribute),
            trait_time,
            trait_event,
        )
        for warning_id, attribute in _WARNING_ATTRIBUTES
    )


def _comparison(
    branch_id: str,
    warning_id: str,
    warning_time: int | None,
    warning_event: FirstPassageEvent | None,
    trait_time: int | None,
    trait_event: FirstPassageEvent | None,
) -> BranchWarningComparison:
    valid = warning_time is not None and trait_time is not None
    leads = None if not valid else warning_time < trait_time
    lead_time = None if not valid else trait_time - warning_time
    return BranchWarningComparison(
        branch_id=branch_id,
        warning_id=warning_id,
        warning_time=warning_time,
        trait_loss_time=trait_time,
        warning_threshold=None if warning_event is None else warning_event.threshold,
        warning_aggregation_rule=None if warning_event is None else warning_event.aggregation_rule,
        trait_threshold=None if trait_event is None else trait_event.threshold,
        trait_aggregation_rule=None if trait_event is None else trait_event.aggregation_rule,
        valid_pair=valid,
        censored=not valid,
        warning_leads=leads,
        lead_time_trait_minus_warning=lead_time,
    )


def _event_by_name(events: Sequence[FirstPassageEvent], name: str) -> FirstPassageEvent | None:
    return next((event for event in events if event.name == name), None)


def _summarise_cell(replicates: Sequence[FiniteH2BranchWarningReplicate]) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    preconditioned = tuple(rep for rep in replicates if rep.finite_h1_mechanism_supported is True)
    return {
        "h1_conditioning": {
            "total_replicates": len(replicates),
            "finite_h1_mechanism_precondition_count": len(preconditioned),
            "finite_h1_mechanism_precondition_probability": len(preconditioned) / len(replicates),
            "unavailable_or_precondition_failed_count": len(replicates) - len(preconditioned),
        },
        "low_start": _summarise_branch(preconditioned, "low_start"),
        "high_start": _summarise_branch(preconditioned, "high_start"),
    }


def _summarise_branch(
    preconditioned: Sequence[FiniteH2BranchWarningReplicate],
    branch_id: str,
) -> dict[str, object]:
    comparisons = tuple(
        comparison
        for replicate in preconditioned
        for comparison in _required_branch(replicate, branch_id)
    )
    by_warning = {
        warning_id: tuple(comparison for comparison in comparisons if comparison.warning_id == warning_id)
        for warning_id, _ in _WARNING_ATTRIBUTES
    }
    return {
        "h1_preconditioned_replicates": len(preconditioned),
        "warning_vs_realised_trait": {
            warning_id: _summarise_warning(values) for warning_id, values in by_warning.items()
        },
    }


def _summarise_warning(values: Sequence[BranchWarningComparison]) -> dict[str, object]:
    if not values:
        return {
            "valid_pair_count": 0,
            "censored_pair_count": 0,
            "lead_count": 0,
            "lead_probability_conditional_on_valid_pair": None,
            "lead_probability_across_h1_preconditioned_replicates": None,
            "lead_time_trait_minus_warning": _empty_summary(),
            "warning_threshold": None,
            "warning_aggregation_rule": None,
            "trait_threshold": None,
            "trait_aggregation_rule": None,
        }
    valid = tuple(value for value in values if value.valid_pair)
    lead_count = sum(value.warning_leads is True for value in valid)
    return {
        "valid_pair_count": len(valid),
        "censored_pair_count": sum(value.censored for value in values),
        "lead_count": lead_count,
        "lead_probability_conditional_on_valid_pair": None if not valid else lead_count / len(valid),
        "lead_probability_across_h1_preconditioned_replicates": lead_count / len(values),
        "lead_time_trait_minus_warning": _summary(
            value.lead_time_trait_minus_warning for value in valid if value.lead_time_trait_minus_warning is not None
        ),
        "warning_threshold": values[0].warning_threshold,
        "warning_aggregation_rule": values[0].warning_aggregation_rule,
        "trait_threshold": values[0].trait_threshold,
        "trait_aggregation_rule": values[0].trait_aggregation_rule,
    }


def _required_branch(
    replicate: FiniteH2BranchWarningReplicate,
    branch_id: str,
) -> tuple[BranchWarningComparison, ...]:
    value = getattr(replicate, branch_id)
    if value is None:
        raise RuntimeError("H1-preconditioned replicate lacks branch-specific H2 comparisons")
    return value


def _summary(values: Iterable[int]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values)
    if not observed:
        return _empty_summary()
    return {
        "mean": sum(observed) / len(observed),
        "median": median(observed),
        "minimum": min(observed),
        "maximum": max(observed),
    }


def _empty_summary() -> dict[str, None]:
    return {"mean": None, "median": None, "minimum": None, "maximum": None}


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
