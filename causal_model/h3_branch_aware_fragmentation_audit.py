"""Branch-aware finite H3 audit for fragmentation and migration modulation.

H3 compares matched total-area landscapes.  H1 and H2 add an important
qualification: spatial subdivision may matter differently depending on whether
the underlying finite system can retain distinct low/high interaction histories
and whether genetic warning events occur on those histories.

This Type S audit derives low/high initial interactions from the *one-large*
canonical H1 branches, then applies those same branch starts to one-large,
equal-isolated, and equal-migrating landscapes using matched parameter cells,
replicate indices, and seeds.  A replicate enters H3 branch-aware contrasts only
when the one-large pair itself retains the finite H1 mechanism predicate.

The multipatch model's migration is an allele-frequency mixing closure.  This
module therefore reports migration modulation of outcomes; it never labels an
outcome demographic rescue or recolonisation.  Those terms remain specific to
the separate extinction--recolonisation lifecycle model.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import CanonicalH1Certificate, canonical_h1_certificate
from causal_model.h1_theorem_boundary_audit import H1TheoremBoundaryAudit, audit_h1_theorem_boundary
from causal_model.multipatch_criticality_dynamics import FirstPassageEvent, SimulationResult, simulate
from causal_model.multipatch_criticality_experiments import (
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
    SCENARIO_ONE_LARGE,
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

_SCENARIOS = (SCENARIO_ONE_LARGE, SCENARIO_EQUAL_ISOLATED, SCENARIO_EQUAL_MIGRATING)
_BRANCHES = ("low_start", "high_start")
_WARNINGS = (("H_alpha", "tau_H_alpha"), ("H_gamma", "tau_H_gamma"), ("allele_loss", "tau_allele_loss"))


@dataclass(frozen=True)
class WarningOrder:
    """One genetic-warning first-passage comparison against realised trait loss."""

    warning_id: str
    warning_time: int | None
    trait_loss_time: int | None
    valid_pair: bool
    censored: bool
    warning_leads: bool | None
    lead_time_trait_minus_warning: int | None
    warning_threshold: float | int | None
    warning_aggregation_rule: str | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class BranchLandscapeOutcome:
    """A branch-start result for one matched H3 landscape."""

    scenario_id: str
    branch_id: str
    summary: ReplicateSummary
    terminal_interaction_mean: float
    terminal_local_effective_size_mean: float
    terminal_potential_high_trait_viable: bool
    warning_orders: tuple[WarningOrder, ...]
    h1_scope: H1TheoremBoundaryAudit

    def as_dict(self) -> dict[str, object]:
        return {
            "scenario_id": self.scenario_id,
            "branch_id": self.branch_id,
            "summary": self.summary.as_dict(),
            "terminal_interaction_mean": self.terminal_interaction_mean,
            "terminal_local_effective_size_mean": self.terminal_local_effective_size_mean,
            "terminal_potential_high_trait_viable": self.terminal_potential_high_trait_viable,
            "warning_orders": [warning.as_dict() for warning in self.warning_orders],
            "h1_scope": self.h1_scope.as_dict(),
        }


@dataclass(frozen=True)
class BranchLandscapeContrast:
    """Outcome differences for a named landscape comparison within one branch."""

    branch_id: str
    comparison_id: str
    interaction_difference: float
    local_effective_size_difference: float
    realised_high_trait_mass_difference: float
    h_alpha_difference: float
    h_gamma_difference: float
    fst_difference: float | None
    potential_high_trait_viability_changed: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


@dataclass(frozen=True)
class H3BranchAwareReplicate:
    """One matched low/high branch and landscape set for one finite replicate."""

    replicate_index: int
    seed: int
    low_initial_interaction: float | None
    high_initial_interaction: float | None
    one_large_finite_h1_mechanism_supported: bool | None
    outcomes: Mapping[str, Mapping[str, BranchLandscapeOutcome]] | None
    branch_retention: Mapping[str, bool] | None
    contrasts: tuple[BranchLandscapeContrast, ...] | None

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "low_initial_interaction": self.low_initial_interaction,
            "high_initial_interaction": self.high_initial_interaction,
            "one_large_finite_h1_mechanism_supported": self.one_large_finite_h1_mechanism_supported,
            "outcomes": None
            if self.outcomes is None
            else {
                scenario_id: {branch_id: outcome.as_dict() for branch_id, outcome in branch_map.items()}
                for scenario_id, branch_map in self.outcomes.items()
            },
            "branch_retention": None if self.branch_retention is None else dict(self.branch_retention),
            "contrasts": None if self.contrasts is None else [contrast.as_dict() for contrast in self.contrasts],
        }


@dataclass(frozen=True)
class H3BranchAwareCell:
    """Branch-aware H3 outcomes for one shared parameter cell."""

    experiment_id: str
    profile: str
    parameters: ParameterCell
    canonical_one_large_h1: CanonicalH1Certificate
    interaction_separation_threshold: float
    replicates: tuple[H3BranchAwareReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "parameters": asdict(self.parameters),
            "canonical_one_large_h1": asdict(self.canonical_one_large_h1),
            "interaction_separation_threshold": self.interaction_separation_threshold,
            "replicate_count": len(self.replicates),
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "replicate_count": len(self.replicates),
            "interaction_separation_threshold": self.interaction_separation_threshold,
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_h3_branch_aware_fragmentation_audit(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    interaction_separation_threshold: float = 0.05,
) -> tuple[H3BranchAwareCell, ...]:
    """Run matched H3 landscapes from the same one-large low/high branch starts.

    The one-large landscape supplies the finite H1 precondition.  It must retain
    low/high interaction separation above `interaction_separation_threshold` and
    show high-start potential high-trait viability while low-start lacks it.
    Only preconditioned replicates contribute H3 contrast and warning-order
    denominators.  Other replicates remain recorded as unavailable rather than
    being interpreted as fragmentation or migration failures.
    """
    if interaction_separation_threshold < 0.0:
        raise ValueError("interaction_separation_threshold must be non-negative")
    scenario_map = _resolve_scenarios(spec, scenarios)
    one_large = scenario_map[SCENARIO_ONE_LARGE]
    cells: list[H3BranchAwareCell] = []
    for cell in parameter_grid(spec):
        base = parameters_for_cell(spec, one_large, cell, seed=spec.master_seed)
        canonical = canonical_h1_certificate(
            feedback_strength=cell.interaction_feedback,
            area=one_large.patch_areas[0],
            area_reference=cell.area_reference,
            barrier=cell.interaction_barrier,
            trait_parameters=base,
        )
        replicates = tuple(
            _run_replicate(
                spec,
                scenario_map,
                cell,
                canonical=canonical,
                replicate_index=index,
                interaction_separation_threshold=interaction_separation_threshold,
            )
            for index in range(spec.replicates)
        )
        cells.append(
            H3BranchAwareCell(
                experiment_id=spec.experiment_id,
                profile=spec.profile,
                parameters=cell,
                canonical_one_large_h1=canonical,
                interaction_separation_threshold=interaction_separation_threshold,
                replicates=replicates,
                summary=_summarise_cell(replicates, canonical),
            )
        )
    return tuple(cells)


def write_h3_branch_aware_fragmentation_artifacts(
    cells: Iterable[H3BranchAwareCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat summaries and full matched branch-aware H3 records."""
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


def _resolve_scenarios(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None,
) -> Mapping[str, LandscapeScenario]:
    selected = tuple(default_scenarios(spec) if scenarios is None else scenarios)
    by_id = {scenario.scenario_id: scenario for scenario in selected}
    if len(by_id) != len(selected):
        raise ValueError("scenario identifiers must be unique")
    missing = sorted(set(_SCENARIOS).difference(by_id))
    extra = sorted(set(by_id).difference(_SCENARIOS))
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing required scenarios: {', '.join(missing)}")
        if extra:
            details.append(f"unexpected scenarios: {', '.join(extra)}")
        raise ValueError("; ".join(details))
    if {scenario.total_area for scenario in by_id.values()} != {spec.total_area}:
        raise ValueError("all H3 scenarios must retain spec.total_area")
    return by_id


def _run_replicate(
    spec: ExperimentSpec,
    scenarios: Mapping[str, LandscapeScenario],
    cell: ParameterCell,
    *,
    canonical: CanonicalH1Certificate,
    replicate_index: int,
    interaction_separation_threshold: float,
) -> H3BranchAwareReplicate:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    if not canonical.branch_dependent_high_trait_mode:
        return H3BranchAwareReplicate(
            replicate_index=replicate_index,
            seed=seed,
            low_initial_interaction=None,
            high_initial_interaction=None,
            one_large_finite_h1_mechanism_supported=None,
            outcomes=None,
            branch_retention=None,
            contrasts=None,
        )
    assert canonical.low_stable_branch is not None
    assert canonical.high_stable_branch is not None
    low_initial = canonical.low_stable_branch.interaction
    high_initial = canonical.high_stable_branch.interaction
    outcomes: dict[str, dict[str, BranchLandscapeOutcome]] = {}
    for scenario_id in _SCENARIOS:
        outcomes[scenario_id] = {
            "low_start": _run_outcome(spec, scenarios[scenario_id], cell, seed, "low_start", low_initial, replicate_index),
            "high_start": _run_outcome(spec, scenarios[scenario_id], cell, seed, "high_start", high_initial, replicate_index),
        }
    large_low = outcomes[SCENARIO_ONE_LARGE]["low_start"]
    large_high = outcomes[SCENARIO_ONE_LARGE]["high_start"]
    one_large_precondition = (
        large_high.terminal_interaction_mean - large_low.terminal_interaction_mean > interaction_separation_threshold
        and large_high.terminal_potential_high_trait_viable
        and not large_low.terminal_potential_high_trait_viable
    )
    if not one_large_precondition:
        return H3BranchAwareReplicate(
            replicate_index=replicate_index,
            seed=seed,
            low_initial_interaction=low_initial,
            high_initial_interaction=high_initial,
            one_large_finite_h1_mechanism_supported=False,
            outcomes=outcomes,
            branch_retention=None,
            contrasts=None,
        )
    branch_retention = {
        scenario_id: _branch_retained(outcomes[scenario_id], interaction_separation_threshold)
        for scenario_id in _SCENARIOS
    }
    contrasts = tuple(
        contrast
        for branch_id in _BRANCHES
        for contrast in (
            _contrast(
                branch_id,
                "equal_isolated_minus_one_large",
                outcomes[SCENARIO_EQUAL_ISOLATED][branch_id],
                outcomes[SCENARIO_ONE_LARGE][branch_id],
            ),
            _contrast(
                branch_id,
                "equal_migrating_minus_equal_isolated",
                outcomes[SCENARIO_EQUAL_MIGRATING][branch_id],
                outcomes[SCENARIO_EQUAL_ISOLATED][branch_id],
            ),
        )
    )
    return H3BranchAwareReplicate(
        replicate_index=replicate_index,
        seed=seed,
        low_initial_interaction=low_initial,
        high_initial_interaction=high_initial,
        one_large_finite_h1_mechanism_supported=True,
        outcomes=outcomes,
        branch_retention=branch_retention,
        contrasts=contrasts,
    )


def _run_outcome(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    seed: int,
    branch_id: str,
    initial_interaction: float,
    replicate_index: int,
) -> BranchLandscapeOutcome:
    parameters = parameters_for_cell(spec, scenario, cell, seed=seed)
    parameters = parameters.__class__(**{**asdict(parameters), "initial_interaction": tuple(initial_interaction for _ in parameters.patch_areas)})
    result = simulate(parameters)
    summary = summarise_replicate(
        result,
        replicate_index=replicate_index,
        seed=seed,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
    )
    final = result.snapshots[-1]
    return BranchLandscapeOutcome(
        scenario_id=scenario.scenario_id,
        branch_id=branch_id,
        summary=summary,
        terminal_interaction_mean=_mean(final.interaction),
        terminal_local_effective_size_mean=_mean(final.effective_size),
        terminal_potential_high_trait_viable=any(
            item.high_trait_component_present for item in final.trait_space
        ),
        warning_orders=_warning_orders(summary),
        h1_scope=audit_h1_theorem_boundary(result),
    )


def _warning_orders(summary: ReplicateSummary) -> tuple[WarningOrder, ...]:
    trait_event = _event_by_name(summary.events, "tau_trait_realised")
    values = []
    for warning_id, attribute in _WARNINGS:
        warning_event = _event_by_name(summary.events, attribute)
        warning_time = getattr(summary, attribute)
        trait_time = summary.tau_trait_realised
        valid = warning_time is not None and trait_time is not None
        values.append(
            WarningOrder(
                warning_id=warning_id,
                warning_time=warning_time,
                trait_loss_time=trait_time,
                valid_pair=valid,
                censored=not valid,
                warning_leads=None if not valid else warning_time < trait_time,
                lead_time_trait_minus_warning=None if not valid else trait_time - warning_time,
                warning_threshold=None if warning_event is None else warning_event.threshold,
                warning_aggregation_rule=None if warning_event is None else warning_event.aggregation_rule,
            )
        )
    return tuple(values)


def _event_by_name(events: Sequence[FirstPassageEvent], name: str) -> FirstPassageEvent | None:
    return next((event for event in events if event.name == name), None)


def _branch_retained(
    outcome_map: Mapping[str, BranchLandscapeOutcome],
    threshold: float,
) -> bool:
    high = outcome_map["high_start"]
    low = outcome_map["low_start"]
    return (
        high.terminal_interaction_mean - low.terminal_interaction_mean > threshold
        and high.terminal_potential_high_trait_viable
        and not low.terminal_potential_high_trait_viable
    )


def _contrast(
    branch_id: str,
    comparison_id: str,
    left: BranchLandscapeOutcome,
    right: BranchLandscapeOutcome,
) -> BranchLandscapeContrast:
    left_summary = left.summary
    right_summary = right.summary
    return BranchLandscapeContrast(
        branch_id=branch_id,
        comparison_id=comparison_id,
        interaction_difference=left.terminal_interaction_mean - right.terminal_interaction_mean,
        local_effective_size_difference=left.terminal_local_effective_size_mean - right.terminal_local_effective_size_mean,
        realised_high_trait_mass_difference=left_summary.realised_high_trait_mass_mean - right_summary.realised_high_trait_mass_mean,
        h_alpha_difference=left_summary.h_alpha - right_summary.h_alpha,
        h_gamma_difference=left_summary.h_gamma - right_summary.h_gamma,
        fst_difference=None
        if left_summary.fst is None or right_summary.fst is None
        else left_summary.fst - right_summary.fst,
        potential_high_trait_viability_changed=(
            left.terminal_potential_high_trait_viable != right.terminal_potential_high_trait_viable
        ),
    )


def _summarise_cell(
    replicates: Sequence[H3BranchAwareReplicate],
    canonical: CanonicalH1Certificate,
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    available = tuple(rep for rep in replicates if rep.one_large_finite_h1_mechanism_supported is True)
    summary: dict[str, object] = {
        "canonical_one_large_context": {
            "gain": canonical.bifurcation.gain,
            "strict_bistability_certified": canonical.bifurcation.strict_bistability_certified,
            "branch_dependent_high_trait_mode": canonical.branch_dependent_high_trait_mode,
        },
        "one_large_h1_conditioning": {
            "total_replicates": len(replicates),
            "finite_h1_precondition_count": len(available),
            "finite_h1_precondition_probability": len(available) / len(replicates),
            "unavailable_or_precondition_failed_count": len(replicates) - len(available),
        },
    }
    summary["branch_retention"] = {
        scenario_id: None if not available else _probability(_required(rep.branch_retention)[scenario_id] for rep in available)
        for scenario_id in _SCENARIOS
    }
    summary["fragmentation_and_migration_contrasts"] = _summarise_contrasts(available)
    summary["branch_conditioned_warning_order"] = _summarise_warnings(available)
    summary["h1_theorem_scope"] = _summarise_scope(available)
    return summary


def _summarise_contrasts(available: Sequence[H3BranchAwareReplicate]) -> dict[str, object]:
    output: dict[str, object] = {}
    for branch_id in _BRANCHES:
        for comparison_id in ("equal_isolated_minus_one_large", "equal_migrating_minus_equal_isolated"):
            contrasts = tuple(
                contrast
                for rep in available
                for contrast in _required(rep.contrasts)
                if contrast.branch_id == branch_id and contrast.comparison_id == comparison_id
            )
            key = f"{branch_id}.{comparison_id}"
            output[key] = _contrast_summary(contrasts)
    return output


def _contrast_summary(values: Sequence[BranchLandscapeContrast]) -> dict[str, object]:
    if not values:
        return {
            "replicate_count": 0,
            "interaction_difference": _empty_summary(),
            "local_effective_size_difference": _empty_summary(),
            "realised_high_trait_mass_difference": _empty_summary(),
            "h_alpha_difference": _empty_summary(),
            "h_gamma_difference": _empty_summary(),
            "fst_difference": _empty_summary(),
            "potential_high_trait_viability_change_probability": None,
        }
    return {
        "replicate_count": len(values),
        "interaction_difference": _summary(value.interaction_difference for value in values),
        "local_effective_size_difference": _summary(value.local_effective_size_difference for value in values),
        "realised_high_trait_mass_difference": _summary(value.realised_high_trait_mass_difference for value in values),
        "h_alpha_difference": _summary(value.h_alpha_difference for value in values),
        "h_gamma_difference": _summary(value.h_gamma_difference for value in values),
        "fst_difference": _summary(value.fst_difference for value in values if value.fst_difference is not None),
        "potential_high_trait_viability_change_probability": _probability(
            value.potential_high_trait_viability_changed for value in values
        ),
    }


def _summarise_warnings(available: Sequence[H3BranchAwareReplicate]) -> dict[str, object]:
    output: dict[str, object] = {}
    for scenario_id in _SCENARIOS:
        for branch_id in _BRANCHES:
            warnings = tuple(
                warning
                for rep in available
                for warning in _required(rep.outcomes)[scenario_id][branch_id].warning_orders
            )
            for warning_id, _ in _WARNINGS:
                values = tuple(warning for warning in warnings if warning.warning_id == warning_id)
                output[f"{scenario_id}.{branch_id}.{warning_id}"] = _warning_summary(values)
    return output


def _warning_summary(values: Sequence[WarningOrder]) -> dict[str, object]:
    if not values:
        return {
            "h1_preconditioned_replicates": 0,
            "valid_pair_count": 0,
            "censored_pair_count": 0,
            "lead_count": 0,
            "lead_probability_conditional_on_valid_pair": None,
            "lead_probability_across_h1_preconditioned_replicates": None,
            "lead_time_trait_minus_warning": _empty_summary(),
        }
    valid = tuple(value for value in values if value.valid_pair)
    lead_count = sum(value.warning_leads is True for value in valid)
    return {
        "h1_preconditioned_replicates": len(values),
        "valid_pair_count": len(valid),
        "censored_pair_count": sum(value.censored for value in values),
        "lead_count": lead_count,
        "lead_probability_conditional_on_valid_pair": None if not valid else lead_count / len(valid),
        "lead_probability_across_h1_preconditioned_replicates": lead_count / len(values),
        "lead_time_trait_minus_warning": _summary(
            value.lead_time_trait_minus_warning for value in valid if value.lead_time_trait_minus_warning is not None
        ),
    }


def _summarise_scope(available: Sequence[H3BranchAwareReplicate]) -> dict[str, object]:
    output: dict[str, object] = {}
    for scenario_id in _SCENARIOS:
        for branch_id in _BRANCHES:
            audits = tuple(
                _required(rep.outcomes)[scenario_id][branch_id].h1_scope for rep in available
            )
            key = f"{scenario_id}.{branch_id}"
            output[key] = {
                "replicate_count": len(audits),
                "patchwise_canonical_update_probability": None
                if not audits
                else _probability(audit.patchwise_canonical_update_certified for audit in audits),
                "single_patch_canonical_theorem_limit_probability": None
                if not audits
                else _probability(audit.single_patch_canonical_theorem_limit_certified for audit in audits),
                "maximum_canonical_update_residual": _summary(
                    audit.maximum_canonical_update_residual for audit in audits
                ),
            }
    return output


def _mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("values must be nonempty")
    return sum(float(value) for value in values) / len(values)


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _summary(values: Iterable[float]) -> dict[str, float | None]:
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


def _required(value):
    if value is None:
        raise RuntimeError("unexpected unavailable branch-aware H3 value")
    return value


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
