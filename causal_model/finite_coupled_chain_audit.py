"""Joint finite-model audit for the declared H1--H2--H3 chain.

The canonical coupled certificate composes separate theorems only after two
explicit ecological closures are supplied.  This module does not assume those
closures.  Instead, it runs the existing finite-bin multipatch simulator under
matched one-large, equal-isolated, and equal-migrating landscapes and records
whether the proposed chain is jointly observed within the declared stochastic
model:

    equal isolation lowers mean interaction
    -> lowers mean local effective size
    -> lowers realised high-trait mass
    -> H-alpha warning precedes realised trait loss in the isolated replicate.

All landscape variants in a replicate use the same seed, parameter cell, and
initialisation rule.  The audit reports the full denominator, valid first-
passage-pair denominator, and counterexamples.  It is Type S evidence for this
finite model, not a theorem and not evidence that every fragmented ecosystem
follows this chain.
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
from causal_model.multipatch_criticality_dynamics import simulate
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

_REQUIRED_SCENARIOS = (
    SCENARIO_ONE_LARGE,
    SCENARIO_EQUAL_ISOLATED,
    SCENARIO_EQUAL_MIGRATING,
)


@dataclass(frozen=True)
class FiniteCoupledChainReplicate:
    """Matched finite outcomes and predeclared joint-chain predicates."""

    replicate_index: int
    seed: int
    one_large: ReplicateSummary
    equal_isolated: ReplicateSummary
    equal_migrating: ReplicateSummary
    h1_scope: Mapping[str, H1TheoremBoundaryAudit]
    fragmentation_mean_interaction_difference: float
    fragmentation_mean_local_effective_size_difference: float
    fragmentation_high_trait_mass_difference: float
    fragmentation_h_alpha_difference: float
    isolated_h_alpha_leads_trait: bool | None
    finite_chain_supported: bool
    migration_mean_interaction_difference: float
    migration_mean_local_effective_size_difference: float
    migration_high_trait_mass_difference: float
    migration_h_alpha_difference: float
    migration_increases_interaction_relative_to_isolation: bool
    migration_increases_local_effective_size_relative_to_isolation: bool
    migration_increases_high_trait_mass_relative_to_isolation: bool
    migration_increases_h_alpha_relative_to_isolation: bool

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "outcomes": {
                SCENARIO_ONE_LARGE: self.one_large.as_dict(),
                SCENARIO_EQUAL_ISOLATED: self.equal_isolated.as_dict(),
                SCENARIO_EQUAL_MIGRATING: self.equal_migrating.as_dict(),
            },
            "h1_scope": {scenario_id: audit.as_dict() for scenario_id, audit in self.h1_scope.items()},
            "fragmentation": {
                "mean_interaction_difference_isolated_minus_large": self.fragmentation_mean_interaction_difference,
                "mean_local_effective_size_difference_isolated_minus_large": self.fragmentation_mean_local_effective_size_difference,
                "high_trait_mass_difference_isolated_minus_large": self.fragmentation_high_trait_mass_difference,
                "h_alpha_difference_isolated_minus_large": self.fragmentation_h_alpha_difference,
                "isolated_h_alpha_leads_trait": self.isolated_h_alpha_leads_trait,
                "finite_chain_supported": self.finite_chain_supported,
            },
            "migration_relative_to_isolation": {
                "mean_interaction_difference_migrating_minus_isolated": self.migration_mean_interaction_difference,
                "mean_local_effective_size_difference_migrating_minus_isolated": self.migration_mean_local_effective_size_difference,
                "high_trait_mass_difference_migrating_minus_isolated": self.migration_high_trait_mass_difference,
                "h_alpha_difference_migrating_minus_isolated": self.migration_h_alpha_difference,
                "increases_interaction": self.migration_increases_interaction_relative_to_isolation,
                "increases_local_effective_size": self.migration_increases_local_effective_size_relative_to_isolation,
                "increases_high_trait_mass": self.migration_increases_high_trait_mass_relative_to_isolation,
                "increases_h_alpha": self.migration_increases_h_alpha_relative_to_isolation,
            },
        }


@dataclass(frozen=True)
class FiniteCoupledChainCell:
    """One parameter cell with canonical context and matched finite-chain outcomes."""

    experiment_id: str
    profile: str
    parameters: ParameterCell
    canonical_h1_context: Mapping[str, CanonicalH1Certificate]
    replicates: tuple[FiniteCoupledChainReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "parameters": asdict(self.parameters),
            "canonical_h1_context": {
                scenario_id: asdict(certificate)
                for scenario_id, certificate in self.canonical_h1_context.items()
            },
            "replicate_count": len(self.replicates),
            "replicates": [replicate.as_dict() for replicate in self.replicates],
            "summary": dict(self.summary),
        }

    def to_csv_row(self) -> dict[str, object]:
        row: dict[str, object] = {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "replicate_count": len(self.replicates),
        }
        row.update(asdict(self.parameters))
        row.update(_flatten_mapping(self.summary))
        return row


def run_finite_coupled_chain_audit(
    spec: ExperimentSpec,
    *,
    scenarios: Sequence[LandscapeScenario] | None = None,
    tolerance: float = 1e-12,
) -> tuple[FiniteCoupledChainCell, ...]:
    """Test the proposed fragmentation-to-warning chain in matched finite runs.

    The returned joint predicate is deliberately stringent and predeclared.  It
    requires, in the same replicate, that equal isolation relative to one large
    lowers mean interaction, mean *local* effective size, and realised high-trait
    mass, and that an H-alpha warning precedes realised trait loss in the
    isolated landscape.  A censored first-passage pair cannot support the joint
    predicate and is reported separately in the summary.
    """
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    selected = _resolve_scenarios(spec, scenarios)
    cells: list[FiniteCoupledChainCell] = []
    for cell in parameter_grid(spec):
        canonical_context = _canonical_h1_context(spec, selected, cell)
        replicates = tuple(
            _run_matched_replicate(spec, selected, cell, index, tolerance=tolerance)
            for index in range(spec.replicates)
        )
        cells.append(
            FiniteCoupledChainCell(
                experiment_id=spec.experiment_id,
                profile=spec.profile,
                parameters=cell,
                canonical_h1_context=canonical_context,
                replicates=replicates,
                summary=_summarise_cell(replicates, canonical_context),
            )
        )
    return tuple(cells)


def write_finite_coupled_chain_artifacts(
    cells: Iterable[FiniteCoupledChainCell],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat cell summaries and full matched finite-chain records."""
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
    values = tuple(default_scenarios(spec) if scenarios is None else scenarios)
    by_id = {scenario.scenario_id: scenario for scenario in values}
    if len(by_id) != len(values):
        raise ValueError("scenario identifiers must be unique")
    missing = sorted(set(_REQUIRED_SCENARIOS).difference(by_id))
    extra = sorted(set(by_id).difference(_REQUIRED_SCENARIOS))
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing required scenarios: {', '.join(missing)}")
        if extra:
            details.append(f"unexpected scenarios: {', '.join(extra)}")
        raise ValueError("; ".join(details))
    total_areas = {scenario.total_area for scenario in by_id.values()}
    if len(total_areas) != 1 or total_areas != {spec.total_area}:
        raise ValueError("all scenarios must retain spec.total_area")
    return by_id


def _canonical_h1_context(
    spec: ExperimentSpec,
    scenarios: Mapping[str, LandscapeScenario],
    cell: ParameterCell,
) -> dict[str, CanonicalH1Certificate]:
    context: dict[str, CanonicalH1Certificate] = {}
    for scenario_id in _REQUIRED_SCENARIOS:
        scenario = scenarios[scenario_id]
        local_area = scenario.patch_areas[0]
        parameters = parameters_for_cell(spec, scenario, cell, seed=spec.master_seed)
        context[scenario_id] = canonical_h1_certificate(
            feedback_strength=cell.interaction_feedback,
            area=local_area,
            area_reference=cell.area_reference,
            barrier=cell.interaction_barrier,
            trait_parameters=parameters,
        )
    return context


def _run_matched_replicate(
    spec: ExperimentSpec,
    scenarios: Mapping[str, LandscapeScenario],
    cell: ParameterCell,
    replicate_index: int,
    *,
    tolerance: float,
) -> FiniteCoupledChainReplicate:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    outcomes: dict[str, ReplicateSummary] = {}
    scopes: dict[str, H1TheoremBoundaryAudit] = {}
    for scenario_id in _REQUIRED_SCENARIOS:
        result = simulate(parameters_for_cell(spec, scenarios[scenario_id], cell, seed=seed))
        outcomes[scenario_id] = summarise_replicate(
            result,
            replicate_index=replicate_index,
            seed=seed,
            h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
            h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
            fst_warning_threshold=spec.fst_warning_threshold,
            allele_loss_threshold=spec.allele_loss_threshold,
        )
        scopes[scenario_id] = audit_h1_theorem_boundary(result)

    large = outcomes[SCENARIO_ONE_LARGE]
    isolated = outcomes[SCENARIO_EQUAL_ISOLATED]
    migrating = outcomes[SCENARIO_EQUAL_MIGRATING]
    fragmentation_interaction = _mean(isolated.final_q_by_patch) - _mean(large.final_q_by_patch)
    fragmentation_effective_size = _mean(isolated.final_effective_size_by_patch) - _mean(large.final_effective_size_by_patch)
    fragmentation_trait = isolated.realised_high_trait_mass_mean - large.realised_high_trait_mass_mean
    fragmentation_h_alpha = isolated.h_alpha - large.h_alpha
    isolated_lead = _strict_h_alpha_lead(isolated)
    chain = (
        fragmentation_interaction < -tolerance
        and fragmentation_effective_size < -tolerance
        and fragmentation_trait < -tolerance
        and isolated_lead is True
    )
    migration_interaction = _mean(migrating.final_q_by_patch) - _mean(isolated.final_q_by_patch)
    migration_effective_size = _mean(migrating.final_effective_size_by_patch) - _mean(isolated.final_effective_size_by_patch)
    migration_trait = migrating.realised_high_trait_mass_mean - isolated.realised_high_trait_mass_mean
    migration_h_alpha = migrating.h_alpha - isolated.h_alpha
    return FiniteCoupledChainReplicate(
        replicate_index=replicate_index,
        seed=seed,
        one_large=large,
        equal_isolated=isolated,
        equal_migrating=migrating,
        h1_scope=scopes,
        fragmentation_mean_interaction_difference=fragmentation_interaction,
        fragmentation_mean_local_effective_size_difference=fragmentation_effective_size,
        fragmentation_high_trait_mass_difference=fragmentation_trait,
        fragmentation_h_alpha_difference=fragmentation_h_alpha,
        isolated_h_alpha_leads_trait=isolated_lead,
        finite_chain_supported=chain,
        migration_mean_interaction_difference=migration_interaction,
        migration_mean_local_effective_size_difference=migration_effective_size,
        migration_high_trait_mass_difference=migration_trait,
        migration_h_alpha_difference=migration_h_alpha,
        migration_increases_interaction_relative_to_isolation=migration_interaction > tolerance,
        migration_increases_local_effective_size_relative_to_isolation=migration_effective_size > tolerance,
        migration_increases_high_trait_mass_relative_to_isolation=migration_trait > tolerance,
        migration_increases_h_alpha_relative_to_isolation=migration_h_alpha > tolerance,
    )


def _summarise_cell(
    replicates: Sequence[FiniteCoupledChainReplicate],
    canonical_context: Mapping[str, CanonicalH1Certificate],
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    isolated_leads = tuple(
        replicate.isolated_h_alpha_leads_trait
        for replicate in replicates
        if replicate.isolated_h_alpha_leads_trait is not None
    )
    return {
        "canonical_h1_context": {
            scenario_id: {
                "gain": certificate.bifurcation.gain,
                "strict_bistability_certified": certificate.bifurcation.strict_bistability_certified,
                "branch_dependent_high_trait_mode": certificate.branch_dependent_high_trait_mode,
            }
            for scenario_id, certificate in canonical_context.items()
        },
        "fragmentation_relative_to_one_large": {
            "mean_interaction_difference_isolated_minus_large": _summary(
                replicate.fragmentation_mean_interaction_difference for replicate in replicates
            ),
            "mean_local_effective_size_difference_isolated_minus_large": _summary(
                replicate.fragmentation_mean_local_effective_size_difference for replicate in replicates
            ),
            "high_trait_mass_difference_isolated_minus_large": _summary(
                replicate.fragmentation_high_trait_mass_difference for replicate in replicates
            ),
            "h_alpha_difference_isolated_minus_large": _summary(
                replicate.fragmentation_h_alpha_difference for replicate in replicates
            ),
            "interaction_lower_probability": _probability(
                replicate.fragmentation_mean_interaction_difference < 0.0 for replicate in replicates
            ),
            "local_effective_size_lower_probability": _probability(
                replicate.fragmentation_mean_local_effective_size_difference < 0.0 for replicate in replicates
            ),
            "high_trait_mass_lower_probability": _probability(
                replicate.fragmentation_high_trait_mass_difference < 0.0 for replicate in replicates
            ),
            "h_alpha_lower_probability": _probability(
                replicate.fragmentation_h_alpha_difference < 0.0 for replicate in replicates
            ),
            "valid_isolated_H_alpha_trait_pairs": len(isolated_leads),
            "censored_isolated_H_alpha_or_trait_pairs": len(replicates) - len(isolated_leads),
            "isolated_H_alpha_lead_conditional": None if not isolated_leads else _probability(isolated_leads),
            "finite_chain_supported_probability": _probability(
                replicate.finite_chain_supported for replicate in replicates
            ),
        },
        "migration_relative_to_isolation": {
            "mean_interaction_difference_migrating_minus_isolated": _summary(
                replicate.migration_mean_interaction_difference for replicate in replicates
            ),
            "mean_local_effective_size_difference_migrating_minus_isolated": _summary(
                replicate.migration_mean_local_effective_size_difference for replicate in replicates
            ),
            "high_trait_mass_difference_migrating_minus_isolated": _summary(
                replicate.migration_high_trait_mass_difference for replicate in replicates
            ),
            "h_alpha_difference_migrating_minus_isolated": _summary(
                replicate.migration_h_alpha_difference for replicate in replicates
            ),
            "interaction_increase_probability": _probability(
                replicate.migration_increases_interaction_relative_to_isolation for replicate in replicates
            ),
            "local_effective_size_increase_probability": _probability(
                replicate.migration_increases_local_effective_size_relative_to_isolation for replicate in replicates
            ),
            "high_trait_mass_increase_probability": _probability(
                replicate.migration_increases_high_trait_mass_relative_to_isolation for replicate in replicates
            ),
            "h_alpha_increase_probability": _probability(
                replicate.migration_increases_h_alpha_relative_to_isolation for replicate in replicates
            ),
        },
        "h1_theorem_scope": {
            scenario_id: _summarise_h1_scope(
                tuple(replicate.h1_scope[scenario_id] for replicate in replicates)
            )
            for scenario_id in _REQUIRED_SCENARIOS
        },
    }


def _strict_h_alpha_lead(summary: ReplicateSummary) -> bool | None:
    if summary.tau_H_alpha is None or summary.tau_trait_realised is None:
        return None
    return summary.tau_H_alpha < summary.tau_trait_realised


def _summarise_h1_scope(audits: Sequence[H1TheoremBoundaryAudit]) -> dict[str, object]:
    return {
        "patchwise_canonical_update_probability": _probability(
            audit.patchwise_canonical_update_certified for audit in audits
        ),
        "single_patch_canonical_theorem_limit_probability": _probability(
            audit.single_patch_canonical_theorem_limit_certified for audit in audits
        ),
        "maximum_canonical_update_residual": _summary(
            audit.maximum_canonical_update_residual for audit in audits
        ),
    }


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
