"""Paired within-simulator ablations for eco-genetic baseline comparisons.

These are *mechanistic ablations*, not claims to reproduce every external
trait-only or population-genetic model.  All variants retain the declared finite
life cycle, landscapes, initial state, and random seed.  They differ only in the
channels allowed to support interaction q and in whether allele state enters
trait recruitment.

For a full support signal
    alpha * q + beta * x_high + gamma * p,
we hold total support weight fixed by replacing a removed channel with current
interaction q.  Thus trait-only uses
    (alpha + gamma) * q + beta * x_high,
and genetic-only uses
    (alpha + beta) * q + gamma * p.
This prevents a simple loss of total input amplitude from masquerading as a
mechanistic ablation effect.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.multipatch_criticality_dynamics import DynamicsParameters, simulate
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    LandscapeScenario,
    ParameterCell,
    ReplicateSummary,
    default_scenarios,
    derived_seed,
    parameter_grid,
    parameters_for_cell,
    quick_profile,
    standard_profile,
    summarise_replicate,
)

BASELINE_TRAIT_ONLY = "trait_only"
BASELINE_GENETIC_ONLY = "genetic_only"
BASELINE_FULL_ECO_GENETIC = "full_eco_genetic"
BASELINE_IDS = (
    BASELINE_TRAIT_ONLY,
    BASELINE_GENETIC_ONLY,
    BASELINE_FULL_ECO_GENETIC,
)


@dataclass(frozen=True)
class BaselineDefinition:
    """Declared causal channels retained by one within-simulator ablation."""

    baseline_id: str
    label: str
    interaction_support: str
    genotype_trait_recruitment: str
    interpretation: str

    def as_dict(self) -> dict[str, str]:
        return asdict(self)


@dataclass(frozen=True)
class PairedBaselineReplicate:
    """Matched outcomes from all ablations using one common random seed."""

    replicate_index: int
    seed: int
    outcomes: Mapping[str, ReplicateSummary]

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "outcomes": {baseline_id: outcome.as_dict() for baseline_id, outcome in self.outcomes.items()},
        }


@dataclass(frozen=True)
class PairedBaselineCell:
    """One landscape × parameter cell with paired baseline outcomes and contrasts."""

    experiment_id: str
    profile: str
    scenario_id: str
    parameters: ParameterCell
    baseline_definitions: tuple[BaselineDefinition, ...]
    replicates: tuple[PairedBaselineReplicate, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": asdict(self.parameters),
            "baseline_definitions": [definition.as_dict() for definition in self.baseline_definitions],
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


def comparison_quick_profile() -> ExperimentSpec:
    """Small coupled profile suitable for baseline-comparison smoke tests.

    The general quick profile intentionally uses a minimal legacy closure with no
    explicit trait/allele q-feedback weights.  This comparison profile instead
    keeps its tiny grid while borrowing the finite coupled closure from the
    standard profile.
    """
    quick = quick_profile()
    coupled = standard_profile().base_parameters
    return replace(
        quick,
        experiment_id="paired_baseline_comparison",
        base_parameters=replace(coupled, trait_grid_size=21),
    )


def resolved_feedback_weights(parameters: DynamicsParameters) -> tuple[float, float, float]:
    """Return exactly the q/trait/allele weights used by the simulator."""
    if parameters.q_feedback_alpha is None and parameters.q_feedback_gamma_allele is None:
        return (
            parameters.interaction_memory_weight,
            parameters.q_feedback_beta_trait,
            1.0 - parameters.interaction_memory_weight,
        )
    alpha = parameters.interaction_memory_weight if parameters.q_feedback_alpha is None else parameters.q_feedback_alpha
    gamma = 0.0 if parameters.q_feedback_gamma_allele is None else parameters.q_feedback_gamma_allele
    return alpha, parameters.q_feedback_beta_trait, gamma


def baseline_definition(parameters: DynamicsParameters, baseline_id: str) -> BaselineDefinition:
    """Describe one causal ablation using resolved support weights."""
    alpha, beta, gamma = resolved_feedback_weights(parameters)
    if baseline_id == BASELINE_TRAIT_ONLY:
        return BaselineDefinition(
            baseline_id=baseline_id,
            label="Trait-only feedback ablation",
            interaction_support=f"({alpha + gamma:g})*q + ({beta:g})*x_high",
            genotype_trait_recruitment="resident_trait_only",
            interpretation="Allele state is retained as an observed drifting/selected state but cannot enter q or recruit trait bins.",
        )
    if baseline_id == BASELINE_GENETIC_ONLY:
        return BaselineDefinition(
            baseline_id=baseline_id,
            label="Genetic-only feedback ablation",
            interaction_support=f"({alpha + beta:g})*q + ({gamma:g})*p",
            genotype_trait_recruitment="resident_trait_only",
            interpretation="Trait occupancy remains an ecological response to q but cannot support q or receive allele-linked recruitment.",
        )
    if baseline_id == BASELINE_FULL_ECO_GENETIC:
        return BaselineDefinition(
            baseline_id=baseline_id,
            label="Full eco-genetic model",
            interaction_support=f"({alpha:g})*q + ({beta:g})*x_high + ({gamma:g})*p",
            genotype_trait_recruitment=parameters.genotype_trait_recruitment,
            interpretation="Declared coupled feedback and genotype-linked trait recruitment are retained exactly.",
        )
    raise ValueError(f"unknown baseline_id: {baseline_id}")


def baseline_parameters(parameters: DynamicsParameters, baseline_id: str) -> DynamicsParameters:
    """Return one calibrated ablation while preserving total q-support weight."""
    alpha, beta, gamma = resolved_feedback_weights(parameters)
    if baseline_id == BASELINE_TRAIT_ONLY:
        return replace(
            parameters,
            q_feedback_alpha=alpha + gamma,
            q_feedback_beta_trait=beta,
            q_feedback_gamma_allele=0.0,
            genotype_trait_recruitment="resident_trait_only",
        )
    if baseline_id == BASELINE_GENETIC_ONLY:
        return replace(
            parameters,
            q_feedback_alpha=alpha + beta,
            q_feedback_beta_trait=0.0,
            q_feedback_gamma_allele=gamma,
            genotype_trait_recruitment="resident_trait_only",
        )
    if baseline_id == BASELINE_FULL_ECO_GENETIC:
        return parameters
    raise ValueError(f"unknown baseline_id: {baseline_id}")


def run_paired_baseline_comparisons(
    spec: ExperimentSpec,
    scenarios: Sequence[LandscapeScenario] | None = None,
    *,
    baseline_ids: Sequence[str] = BASELINE_IDS,
) -> tuple[PairedBaselineCell, ...]:
    """Run matched trait-only, genetic-only, and full-model ablations.

    Each baseline in a replicate receives the same derived seed and differs only
    through ``baseline_parameters``.  Shared seeds make pairwise contrasts more
    stable; they do not make the stochastic paths identical after mechanisms
    cause the trajectories to diverge.
    """
    selected_scenarios = tuple(scenarios) if scenarios is not None else default_scenarios(spec)
    selected_baselines = _validated_baseline_ids(baseline_ids)
    if not selected_scenarios:
        raise ValueError("scenarios must be nonempty")
    cells: list[PairedBaselineCell] = []
    for scenario in selected_scenarios:
        for cell in parameter_grid(spec):
            common_parameters = parameters_for_cell(spec, scenario, cell, seed=spec.master_seed)
            definitions = tuple(baseline_definition(common_parameters, baseline_id) for baseline_id in selected_baselines)
            replicates = tuple(
                _run_paired_replicate(spec, scenario, cell, index, selected_baselines)
                for index in range(spec.replicates)
            )
            cells.append(
                PairedBaselineCell(
                    experiment_id=spec.experiment_id,
                    profile=spec.profile,
                    scenario_id=scenario.scenario_id,
                    parameters=cell,
                    baseline_definitions=definitions,
                    replicates=replicates,
                    summary=_summarise_paired_replicates(replicates, selected_baselines),
                )
            )
    return tuple(cells)


def write_paired_baseline_artifacts(
    cells: Iterable[PairedBaselineCell], *, csv_path: str | Path, json_path: str | Path
) -> None:
    """Write flat comparison summaries and full paired replicate records."""
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


def _run_paired_replicate(
    spec: ExperimentSpec,
    scenario: LandscapeScenario,
    cell: ParameterCell,
    replicate_index: int,
    baseline_ids: Sequence[str],
) -> PairedBaselineReplicate:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    common = parameters_for_cell(spec, scenario, cell, seed=seed)
    outcomes: dict[str, ReplicateSummary] = {}
    for baseline_id in baseline_ids:
        result = simulate(baseline_parameters(common, baseline_id))
        outcomes[baseline_id] = summarise_replicate(
            result,
            replicate_index=replicate_index,
            seed=seed,
            h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
            h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
            fst_warning_threshold=spec.fst_warning_threshold,
            allele_loss_threshold=spec.allele_loss_threshold,
        )
    return PairedBaselineReplicate(replicate_index=replicate_index, seed=seed, outcomes=outcomes)


def _summarise_paired_replicates(
    replicates: Sequence[PairedBaselineReplicate], baseline_ids: Sequence[str]
) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    by_baseline = {
        baseline_id: tuple(replicate.outcomes[baseline_id] for replicate in replicates)
        for baseline_id in baseline_ids
    }
    summary: dict[str, object] = {
        "models": {baseline_id: _summarise_model(outcomes) for baseline_id, outcomes in by_baseline.items()}
    }
    if BASELINE_FULL_ECO_GENETIC in by_baseline:
        full = by_baseline[BASELINE_FULL_ECO_GENETIC]
        contrasts: dict[str, object] = {}
        for baseline_id in (BASELINE_TRAIT_ONLY, BASELINE_GENETIC_ONLY):
            if baseline_id in by_baseline:
                contrasts[f"full_minus_{baseline_id}"] = _paired_contrast(full, by_baseline[baseline_id])
        summary["paired_contrasts"] = contrasts
    return summary


def _summarise_model(outcomes: Sequence[ReplicateSummary]) -> dict[str, float | int | None]:
    return {
        "realised_high_trait_persistence_final": _probability(outcome.realised_high_trait_persists for outcome in outcomes),
        "potential_high_trait_viability_final": _probability(outcome.potential_high_trait_viable for outcome in outcomes),
        "realised_high_trait_mass_mean_final": _mean(outcome.realised_high_trait_mass_mean for outcome in outcomes),
        "final_h_alpha_mean": _mean(outcome.h_alpha for outcome in outcomes),
        "final_h_gamma_mean": _mean(outcome.h_gamma for outcome in outcomes),
        "final_high_allele_frequency_mean": _mean(
            _mean(outcome.final_p_by_patch) for outcome in outcomes
        ),
        "genetic_lead_H_alpha_conditional": _conditional_lead_probability(outcomes, "tau_H_alpha"),
        "valid_H_alpha_trait_pairs": _valid_pair_count(outcomes, "tau_H_alpha"),
        "censored_H_alpha_or_trait_pairs": len(outcomes) - _valid_pair_count(outcomes, "tau_H_alpha"),
    }


def _paired_contrast(full: Sequence[ReplicateSummary], baseline: Sequence[ReplicateSummary]) -> dict[str, float | int]:
    if len(full) != len(baseline):
        raise ValueError("paired outcomes must have equal length")
    persistence_difference = tuple(
        int(a.realised_high_trait_persists) - int(b.realised_high_trait_persists)
        for a, b in zip(full, baseline)
    )
    mass_difference = tuple(
        a.realised_high_trait_mass_mean - b.realised_high_trait_mass_mean
        for a, b in zip(full, baseline)
    )
    alpha_difference = tuple(a.h_alpha - b.h_alpha for a, b in zip(full, baseline))
    allele_difference = tuple(
        _mean(a.final_p_by_patch) - _mean(b.final_p_by_patch)
        for a, b in zip(full, baseline)
    )
    return {
        "replicate_count": len(full),
        "realised_high_trait_persistence_difference_mean": _mean(persistence_difference),
        "realised_high_trait_mass_difference_mean": _mean(mass_difference),
        "final_h_alpha_difference_mean": _mean(alpha_difference),
        "final_high_allele_frequency_difference_mean": _mean(allele_difference),
        "full_trait_mass_greater_probability": _probability(value > 0.0 for value in mass_difference),
        "full_trait_mass_tie_probability": _probability(value == 0.0 for value in mass_difference),
        "full_trait_mass_less_probability": _probability(value < 0.0 for value in mass_difference),
    }


def _validated_baseline_ids(baseline_ids: Sequence[str]) -> tuple[str, ...]:
    values = tuple(dict.fromkeys(baseline_ids))
    if not values:
        raise ValueError("baseline_ids must be nonempty")
    unknown = sorted(set(values).difference(BASELINE_IDS))
    if unknown:
        raise ValueError(f"unknown baseline_ids: {', '.join(unknown)}")
    return values


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


def _mean(values: Iterable[float]) -> float:
    observed = tuple(float(value) for value in values)
    if not observed:
        raise ValueError("values must be nonempty")
    return sum(observed) / len(observed)


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
