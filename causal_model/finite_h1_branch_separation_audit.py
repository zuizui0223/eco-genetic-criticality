"""Audit canonical H1 branch separation in the finite coupled simulator.

The canonical H1 theorem proves bistability only for its one-state reduction.
The finite-bin multipatch model adds density variation, trait/allele feedback,
finite recruitment, and optional migration. This module therefore tests a
predeclared finite-model analogue rather than claiming a second theorem.

For each canonical-H1 cell with a low and high stable branch, otherwise matched
finite simulations start from the canonical low versus high interaction values.
They share landscape, trait/allele initial state, parameter cell, replicate
index, and random seed. The audit asks whether terminal interaction remains
separated and whether potential high-trait viability switches with the initial
branch.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import CanonicalH1Certificate, canonical_h1_certificate
from causal_model.h1_theorem_boundary_audit import H1TheoremBoundaryAudit, audit_h1_theorem_boundary
from causal_model.multipatch_criticality_dynamics import SimulationResult, simulate
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


@dataclass(frozen=True)
class FiniteH1BranchPair:
    """One same-seed finite low-start/high-start comparison."""

    replicate_index: int
    seed: int
    low_initial_interaction: float | None
    high_initial_interaction: float | None
    low_start: ReplicateSummary | None
    high_start: ReplicateSummary | None
    low_scope: H1TheoremBoundaryAudit | None
    high_scope: H1TheoremBoundaryAudit | None
    low_terminal_interaction_mean: float | None
    high_terminal_interaction_mean: float | None
    terminal_interaction_difference_high_minus_low: float | None
    terminal_high_trait_mass_difference_high_minus_low: float | None
    terminal_local_effective_size_difference_high_minus_low: float | None
    potential_high_trait_switch: bool | None
    interaction_branch_separation_supported: bool | None
    finite_h1_mechanism_supported: bool | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "low_initial_interaction": self.low_initial_interaction,
            "high_initial_interaction": self.high_initial_interaction,
            "low_start": None if self.low_start is None else self.low_start.as_dict(),
            "high_start": None if self.high_start is None else self.high_start.as_dict(),
            "low_scope": None if self.low_scope is None else self.low_scope.as_dict(),
            "high_scope": None if self.high_scope is None else self.high_scope.as_dict(),
            "finite_branch_audit": {
                "low_terminal_interaction_mean": self.low_terminal_interaction_mean,
                "high_terminal_interaction_mean": self.high_terminal_interaction_mean,
                "terminal_interaction_difference_high_minus_low": self.terminal_interaction_difference_high_minus_low,
                "terminal_high_trait_mass_difference_high_minus_low": self.terminal_high_trait_mass_difference_high_minus_low,
                "terminal_local_effective_size_difference_high_minus_low": self.terminal_local_effective_size_difference_high_minus_low,
                "potential_high_trait_switch": self.potential_high_trait_switch,
                "interaction_branch_separation_supported": self.interaction_branch_separation_supported,
                "finite_h1_mechanism_supported": self.finite_h1_mechanism_supported,
            },
        }


@dataclass(frozen=True)
class FiniteH1BranchSeparationCell:
    """Finite branch-pair evidence for one scenario and parameter cell."""

    experiment_id: str
    profile: str
    scenario_id: str
    parameters: ParameterCell
    canonical_h1: CanonicalH1Certificate
    interaction_separation_threshold: float
    terminal_window: int
    replicates: tuple[FiniteH1BranchPair, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": asdict(self.parameters),
            "canonical_h1": asdict(self.canonical_h1),
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "terminal_window": self.terminal_window,
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
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "terminal_window": self.terminal_window,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_finite_h1_branch_separation_audit(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    interaction_separation_threshold: float = 0.05,
    terminal_window: int = 5,
) -> tuple[FiniteH1BranchSeparationCell, ...]:
    """Run same-seed low-start/high-start finite tests across declared cells.

    `finite_h1_mechanism_supported` requires all of the following in a replicate:
    a canonical branch-dependent high-trait certificate, terminal mean interaction
    separation above the declared threshold, and a potential high-trait switch
    (high-start viable; low-start non-viable). Realised trait and effective-size
    contrasts are retained as finite-model outcomes rather than theorem premises.
    """
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    if terminal_window < 1:
        raise ValueError("terminal_window must be positive")
    selected = tuple(default_scenarios(spec) if scenarios is None else scenarios)
    if not selected:
        raise ValueError("scenarios must be nonempty")
    if len({scenario.scenario_id for scenario in selected}) != len(selected):
        raise ValueError("scenario identifiers must be unique")

    cells: list[FiniteH1BranchSeparationCell] = []
    for scenario in selected:
        for cell in parameter_grid(spec):
            base = parameters_for_cell(spec, scenario, cell, seed=spec.master_seed)
            canonical = canonical_h1_certificate(
                feedback_strength=cell.interaction_feedback,
                area=scenario.patch_areas[0],
                area_reference=cell.area_reference,
                barrier=cell.interaction_barrier,
                trait_parameters=base,
            )
            pairs = tuple(
                _run_branch_pair(
                    spec,
                    scenario,
                    cell,
                    replicate_index=index,
                    canonical=canonical,
                    interaction_separation_threshold=interaction_separation_threshold,
                    terminal_window=terminal_window,
                )
                for index in range(spec.replicates)
            )
            cells.append(
                FiniteH1BranchSeparationCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    scenario_id=scenario.scenario_id,
                    parameters=cell,
                    canonical_h1=canonical,
                    interaction_separation_threshold=interaction_separation_threshold,
                    terminal_window=terminal_window,
                    replicates=pairs,
                    summary=_summarise_pairs(pairs, canonical),
                )
            )
    return tuple(cells)


def write_finite_h1_branch_separation_artifacts(
    cells: Iterable[FiniteH1BranchSeparationCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat summaries and complete same-seed low/high branch records."""
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


def _run_branch_pair(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    *,
    replicate_index: int,
    canonical: CanonicalH1Certificate,
    interaction_separation_threshold: float,
    terminal_window: int,
) -> FiniteH1BranchPair:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    if not canonical.branch_dependent_high_trait_mode:
        return FiniteH1BranchPair(
            replicate_index=replicate_index,
            seed=seed,
            low_initial_interaction=None,
            high_initial_interaction=None,
            low_start=None,
            high_start=None,
            low_scope=None,
            high_scope=None,
            low_terminal_interaction_mean=None,
            high_terminal_interaction_mean=None,
            terminal_interaction_difference_high_minus_low=None,
            terminal_high_trait_mass_difference_high_minus_low=None,
            terminal_local_effective_size_difference_high_minus_low=None,
            potential_high_trait_switch=None,
            interaction_branch_separation_supported=None,
            finite_h1_mechanism_supported=None,
        )

    assert canonical.low_stable_branch is not None
    assert canonical.high_stable_branch is not None
    low_initial = canonical.low_stable_branch.interaction
    high_initial = canonical.high_stable_branch.interaction
    base = parameters_for_cell(spec, scenario, cell, seed=seed)
    low_result = simulate(_with_uniform_initial_interaction(base, low_initial))
    high_result = simulate(_with_uniform_initial_interaction(base, high_initial))
    low_summary = _summary_for_result(spec, low_result, replicate_index, seed)
    high_summary = _summary_for_result(spec, high_result, replicate_index, seed)
    low_tail = _terminal_interaction_mean(low_result, terminal_window)
    high_tail = _terminal_interaction_mean(high_result, terminal_window)
    interaction_difference = high_tail - low_tail
    trait_difference = high_summary.realised_high_trait_mass_mean - low_summary.realised_high_trait_mass_mean
    effective_size_difference = _mean(high_summary.final_effective_size_by_patch) - _mean(
        low_summary.final_effective_size_by_patch
    )
    potential_switch = high_summary.potential_high_trait_viable and not low_summary.potential_high_trait_viable
    separated = interaction_difference > interaction_separation_threshold
    return FiniteH1BranchPair(
        replicate_index=replicate_index,
        seed=seed,
        low_initial_interaction=low_initial,
        high_initial_interaction=high_initial,
        low_start=low_summary,
        high_start=high_summary,
        low_scope=audit_h1_theorem_boundary(low_result),
        high_scope=audit_h1_theorem_boundary(high_result),
        low_terminal_interaction_mean=low_tail,
        high_terminal_interaction_mean=high_tail,
        terminal_interaction_difference_high_minus_low=interaction_difference,
        terminal_high_trait_mass_difference_high_minus_low=trait_difference,
        terminal_local_effective_size_difference_high_minus_low=effective_size_difference,
        potential_high_trait_switch=potential_switch,
        interaction_branch_separation_supported=separated,
        finite_h1_mechanism_supported=separated and potential_switch,
    )


def _with_uniform_initial_interaction(parameters, interaction: float):
    return replace(parameters, initial_interaction=tuple(interaction for _ in parameters.patch_areas))


def _summary_for_result(
    spec: ExperimentSpec,
    result: SimulationResult,
    replicate_index: int,
    seed: int,
) -> ReplicateSummary:
    return summarise_replicate(
        result,
        replicate_index=replicate_index,
        seed=seed,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
    )


def _terminal_interaction_mean(result: SimulationResult, terminal_window: int) -> float:
    snapshots = result.snapshots[-min(terminal_window, len(result.snapshots)) :]
    return _mean(tuple(_mean(snapshot.interaction) for snapshot in snapshots))


def _summarise_pairs(
    pairs: Sequence[FiniteH1BranchPair],
    canonical: CanonicalH1Certificate,
) -> dict[str, object]:
    available = tuple(pair for pair in pairs if pair.finite_h1_mechanism_supported is not None)
    summary: dict[str, object] = {
        "canonical_context": _canonical_context(canonical),
        "finite_pair_available_probability": len(available) / len(pairs),
    }
    if not available:
        summary.update(
            {
                "interaction_branch_separation_probability": None,
                "potential_high_trait_switch_probability": None,
                "finite_h1_mechanism_supported_probability": None,
                "terminal_interaction_difference_high_minus_low": _empty_summary(),
                "terminal_high_trait_mass_difference_high_minus_low": _empty_summary(),
                "terminal_local_effective_size_difference_high_minus_low": _empty_summary(),
                "h1_theorem_scope": {"low_start": None, "high_start": None},
            }
        )
        return summary

    summary.update(
        {
            "interaction_branch_separation_probability": _probability(
                bool(pair.interaction_branch_separation_supported) for pair in available
            ),
            "potential_high_trait_switch_probability": _probability(
                bool(pair.potential_high_trait_switch) for pair in available
            ),
            "finite_h1_mechanism_supported_probability": _probability(
                bool(pair.finite_h1_mechanism_supported) for pair in available
            ),
            "terminal_interaction_difference_high_minus_low": _summary(
                _required(pair.terminal_interaction_difference_high_minus_low) for pair in available
            ),
            "terminal_high_trait_mass_difference_high_minus_low": _summary(
                _required(pair.terminal_high_trait_mass_difference_high_minus_low) for pair in available
            ),
            "terminal_local_effective_size_difference_high_minus_low": _summary(
                _required(pair.terminal_local_effective_size_difference_high_minus_low) for pair in available
            ),
            "h1_theorem_scope": {
                "low_start": _summarise_scope(_required(pair.low_scope) for pair in available),
                "high_start": _summarise_scope(_required(pair.high_scope) for pair in available),
            },
        }
    )
    return summary


def _canonical_context(certificate: CanonicalH1Certificate) -> dict[str, object]:
    return {
        "gain": certificate.bifurcation.gain,
        "strict_bistability_certified": certificate.bifurcation.strict_bistability_certified,
        "branch_dependent_high_trait_mode": certificate.branch_dependent_high_trait_mode,
        "low_branch_interaction": None if certificate.low_stable_branch is None else certificate.low_stable_branch.interaction,
        "high_branch_interaction": None if certificate.high_stable_branch is None else certificate.high_stable_branch.interaction,
    }


def _summarise_scope(audits: Iterable[H1TheoremBoundaryAudit]) -> dict[str, object]:
    values = tuple(audits)
    return {
        "patchwise_canonical_update_probability": _probability(
            audit.patchwise_canonical_update_certified for audit in values
        ),
        "single_patch_canonical_theorem_limit_probability": _probability(
            audit.single_patch_canonical_theorem_limit_certified for audit in values
        ),
        "maximum_canonical_update_residual": _summary(
            audit.maximum_canonical_update_residual for audit in values
        ),
    }


def _required(value):
    if value is None:
        raise RuntimeError("unexpected missing value for available finite branch pair")
    return value


def _mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("values must be nonempty")
    return sum(float(value) for value in values) / len(values)


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _summary(values: Iterable[float]) -> dict[str, float]:
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
