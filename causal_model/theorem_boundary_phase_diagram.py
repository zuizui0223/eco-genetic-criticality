"""Phase-diagram artifacts that retain H1 theorem scope beside stochastic outcomes.

This wrapper deliberately reuses the existing multipatch experiment grid, seed
schedule, and first-passage summaries.  Each replicate additionally records the
H1 canonical-update residual and named assumptions that have been relaxed.
Thus an outcome map can be read as a theorem-limit result, a controlled
departure, or an unrestricted stochastic result without conflating the three.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.h1_theorem_boundary_audit import H1TheoremBoundaryAudit, audit_h1_theorem_boundary
from causal_model.multipatch_criticality_dynamics import simulate
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    ReplicateSummary,
    default_scenarios,
    derived_seed,
    parameter_grid,
    parameters_for_cell,
    summarise_replicate,
)

DEPARTURE_LABELS = (
    "density_not_one",
    "trait_feedback_enabled",
    "allele_feedback_enabled",
    "support_not_equal_interaction",
    "migration_enabled",
    "multiple_patches",
)


@dataclass(frozen=True)
class TheoremBoundaryReplicate:
    """One stochastic outcome paired with its H1 theorem-boundary audit."""

    replicate_index: int
    seed: int
    outcome: ReplicateSummary
    audit: H1TheoremBoundaryAudit

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "outcome": self.outcome.as_dict(),
            "h1_theorem_boundary": self.audit.as_dict(),
        }


@dataclass(frozen=True)
class TheoremBoundaryCell:
    """One landscape × parameter cell with outcome and scope summaries."""

    experiment_id: str
    profile: str
    scenario_id: str
    parameters: ParameterCell
    replicates: tuple[TheoremBoundaryReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": asdict(self.parameters),
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
        row.update(_flatten_mapping(self.summary))
        return row


def run_theorem_boundary_phase_diagram(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    tolerance: float = 1e-12,
) -> tuple[TheoremBoundaryCell, ...]:
    """Run the declared stochastic grid with a theorem-boundary audit per replicate."""
    selected_scenarios = tuple(scenarios) if scenarios is not None else default_scenarios(spec)
    if not selected_scenarios:
        raise ValueError("scenarios must be nonempty")
    cells: list[TheoremBoundaryCell] = []
    for scenario in selected_scenarios:
        for cell in parameter_grid(spec):
            replicates = tuple(
                _run_boundary_replicate(spec, scenario, cell, replicate_index, tolerance=tolerance)
                for replicate_index in range(spec.replicates)
            )
            cells.append(
                TheoremBoundaryCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    scenario_id=scenario.scenario_id,
                    parameters=cell,
                    replicates=replicates,
                    summary=_summarise_boundary_replicates(replicates),
                )
            )
    return tuple(cells)


def write_theorem_boundary_phase_artifacts(
    cells: Iterable[TheoremBoundaryCell], *, csv_path: str | Path, json_path: str | Path) -> None:
    """Write flat CSV and full JSON artifacts for reproducible robustness figures."""
    values = tuple(cells)
    if not values:
        raise ValueError("cells must be nonempty")
    csv_destination = Path(csv_path)
    json_destination = Path(json_path)
    csv_destination.parent.mkdir(parents=True, exist_ok=True)
    json_destination.parent.mkdir(parents=True, exist_ok=True)
    rows = [cell.to_csv_row() for cell in values]
    fieldnames = sorted({key for row in rows for key in row})
    with csv_destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with json_destination.open("w", encoding="utf-8") as handle:
        json.dump([cell.as_dict() for cell in values], handle, indent=2, sort_keys=True)


def _run_boundary_replicate(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    replicate_index: int,
    *,
    tolerance: float,
) -> TheoremBoundaryReplicate:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    result = simulate(parameters_for_cell(spec, scenario, cell, seed=seed))
    outcome = summarise_replicate(
        result,
        replicate_index=replicate_index,
        seed=seed,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
    )
    return TheoremBoundaryReplicate(
        replicate_index=replicate_index,
        seed=seed,
        outcome=outcome,
        audit=audit_h1_theorem_boundary(result, tolerance=tolerance),
    )


def _summarise_boundary_replicates(replicates: Sequence[TheoremBoundaryReplicate]) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    maximum_residuals = tuple(rep.audit.maximum_canonical_update_residual for rep in replicates)
    mean_residuals = tuple(rep.audit.mean_canonical_update_residual for rep in replicates)
    outcomes = tuple(rep.outcome for rep in replicates)
    return {
        "scope": {
            "patchwise_canonical_update_probability": _probability(
                rep.audit.patchwise_canonical_update_certified for rep in replicates
            ),
            "single_patch_canonical_theorem_limit_probability": _probability(
                rep.audit.single_patch_canonical_theorem_limit_certified for rep in replicates
            ),
            "maximum_canonical_update_residual": _summary(maximum_residuals),
            "mean_canonical_update_residual": _summary(mean_residuals),
            "maximum_density_deviation_from_one": _summary(
                rep.audit.maximum_density_deviation_from_one for rep in replicates
            ),
            "maximum_support_deviation_from_interaction": _summary(
                rep.audit.maximum_support_deviation_from_interaction for rep in replicates
            ),
            "departure_probabilities": {
                label: _probability(label in rep.audit.departure_labels for rep in replicates)
                for label in DEPARTURE_LABELS
            },
        },
        "outcomes": {
            "realised_high_trait_persistence_final": _probability(
                rep.realised_high_trait_persists for rep in outcomes
            ),
            "potential_high_trait_viability_final": _probability(
                rep.potential_high_trait_viable for rep in outcomes
            ),
            "genetic_lead_H_alpha_conditional": _conditional_lead_probability(outcomes, "tau_H_alpha"),
            "genetic_lead_H_gamma_conditional": _conditional_lead_probability(outcomes, "tau_H_gamma"),
            "allele_loss_lead_conditional": _conditional_lead_probability(outcomes, "tau_allele_loss"),
            "valid_H_alpha_trait_pairs": _valid_pair_count(outcomes, "tau_H_alpha"),
            "censored_H_alpha_or_trait_pairs": len(outcomes) - _valid_pair_count(outcomes, "tau_H_alpha"),
        },
    }


def _conditional_lead_probability(replicates: Sequence[ReplicateSummary], warning_attribute: str) -> float | None:
    valid = tuple(
        getattr(rep, warning_attribute) < rep.tau_trait_realised
        for rep in replicates
        if getattr(rep, warning_attribute) is not None and rep.tau_trait_realised is not None
    )
    return None if not valid else _probability(valid)


def _valid_pair_count(replicates: Sequence[ReplicateSummary], warning_attribute: str) -> int:
    return sum(
        getattr(rep, warning_attribute) is not None and rep.tau_trait_realised is not None
        for rep in replicates
    )


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _summary(values: Iterable[float]) -> dict[str, float]:
    observed = tuple(float(value) for value in values)
    if not observed:
        raise ValueError("values must be nonempty")
    return {
        "mean": sum(observed) / len(observed),
        "median": median(observed),
        "minimum": min(observed),
        "maximum": max(observed),
    }


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
