"""Censoring-aware event rows for finite-bin phase-diagram artifacts.

This module leaves the legacy ``CellResult.to_csv_row`` contract untouched and
adds a parallel row format in which every warning-versus-trait comparison
retains its denominator.  It is intended for standard/full experiment outputs.
"""
from __future__ import annotations

from typing import Iterable, Sequence

from causal_model.first_passage_reporting import comparison_as_dict, compare_first_passage_times
from causal_model.multipatch_criticality_experiments import CellResult, ReplicateSummary


WARNING_EVENTS: tuple[str, ...] = (
    "tau_H_alpha",
    "tau_H_gamma",
    "tau_FST",
    "tau_allele_loss",
)
REFERENCE_EVENT = "tau_trait_realised"


def event_comparison_for_replicates(
    replicates: Sequence[ReplicateSummary],
    warning_event: str,
    *,
    reference_event: str = REFERENCE_EVENT,
) -> dict[str, object]:
    """Summarise one warning event against realised trait loss across replicates."""
    if warning_event not in WARNING_EVENTS:
        raise ValueError(f"unknown warning event: {warning_event}")
    warnings = tuple(getattr(replicate, warning_event) for replicate in replicates)
    references = tuple(getattr(replicate, reference_event) for replicate in replicates)
    return comparison_as_dict(compare_first_passage_times(warnings, references))


def censoring_aware_cell_row(cell: CellResult) -> dict[str, object]:
    """Return a flat artifact row with all H2 comparison denominators visible."""
    row: dict[str, object] = {
        "experiment_id": cell.experiment_id,
        "profile": cell.profile,
        "scenario_id": cell.scenario_id,
        "replicate_count": len(cell.replicates),
    }
    row.update(cell.parameters.as_dict())
    for warning_event in WARNING_EVENTS:
        comparison = event_comparison_for_replicates(cell.replicates, warning_event)
        prefix = f"{warning_event}_vs_{REFERENCE_EVENT}"
        for key, value in comparison.items():
            if key == "time_differences":
                row[f"{prefix}.{key}"] = ";".join(str(item) for item in value)
            else:
                row[f"{prefix}.{key}"] = value
    return row


def censoring_aware_phase_diagram_rows(results: Iterable[CellResult]) -> tuple[dict[str, object], ...]:
    """Generate one censoring-aware event row per parameter/scenario cell."""
    return tuple(censoring_aware_cell_row(cell) for cell in results)
