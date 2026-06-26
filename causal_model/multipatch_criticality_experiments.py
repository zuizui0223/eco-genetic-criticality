"""Reproducible phase-diagram experiments for realised occupancy hypotheses.

This is a simulation/reporting layer. It does not alter the theorem layer, and it
keeps potential viability, realised occupancy, allele persistence, and diversity
statistics distinct.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from itertools import product
from statistics import median
from typing import Iterable, Mapping, Sequence

from causal_model.multipatch_criticality_dynamics import (
    DynamicsParameters,
    FirstPassageEvent,
    SimulationResult,
    first_passage_events,
    simulate,
    tau_FST,
    tau_H_alpha,
    tau_H_gamma,
    tau_allele_loss,
    tau_trait_potential,
    tau_trait_realised,
)

SCENARIO_ONE_LARGE = "one_large"
SCENARIO_EQUAL_ISOLATED = "equal_isolated"
SCENARIO_EQUAL_MIGRATING = "equal_migrating"
PROFILE_QUICK = "quick"
PROFILE_STANDARD = "standard"
PROFILE_FULL = "full"


@dataclass(frozen=True)
class ExperimentSpec:
    """Immutable specification for declared stochastic experiments."""

    experiment_id: str = "realised_occupancy_phase_diagram"
    profile: str = PROFILE_QUICK
    total_area: float = 4.0
    patch_count: int = 4
    generations: int = 20
    replicates: int = 4
    master_seed: int = 1
    area_reference_values: tuple[float, ...] = (1.0,)
    interaction_feedback_values: tuple[float, ...] = (4.0,)
    interaction_barrier_values: tuple[float, ...] = (0.45, 0.65)
    migration_rate: float = 0.1
    h_alpha_warning_threshold: float = 0.2
    h_gamma_warning_threshold: float = 0.2
    fst_warning_threshold: float = 0.2
    allele_loss_threshold: float = 0.0
    base_parameters: DynamicsParameters = DynamicsParameters(patch_areas=(1.0,))

    def __post_init__(self) -> None:
        if self.profile not in {PROFILE_QUICK, PROFILE_STANDARD, PROFILE_FULL}:
            raise ValueError("unknown experiment profile")
        if self.total_area <= 0.0 or self.patch_count < 1:
            raise ValueError("total_area must be positive and patch_count at least one")
        if self.generations < 1 or self.replicates < 1:
            raise ValueError("generations and replicates must be positive")
        for value, label in (
            (self.migration_rate, "migration_rate"),
            (self.h_alpha_warning_threshold, "h_alpha_warning_threshold"),
            (self.h_gamma_warning_threshold, "h_gamma_warning_threshold"),
            (self.fst_warning_threshold, "fst_warning_threshold"),
            (self.allele_loss_threshold, "allele_loss_threshold"),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{label} must lie in [0, 1]")
        if not all((self.area_reference_values, self.interaction_feedback_values, self.interaction_barrier_values)):
            raise ValueError("parameter-grid values must be nonempty")


@dataclass(frozen=True)
class ParameterCell:
    cell_index: int
    area_reference: float
    interaction_feedback: float
    interaction_barrier: float

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class LandscapeScenario:
    scenario_id: str
    patch_areas: tuple[float, ...]
    migration_rate: float

    @property
    def total_area(self) -> float:
        return sum(self.patch_areas)


@dataclass(frozen=True)
class ReplicateSummary:
    """One replicate with raw first-passage times and explicit event metadata."""

    replicate_index: int
    seed: int
    final_q_by_patch: tuple[float, ...]
    final_population_by_patch: tuple[int, ...]
    final_effective_size_by_patch: tuple[float, ...]
    final_p_by_patch: tuple[float, ...]
    final_high_trait_abundance_by_patch: tuple[int, ...]
    h_alpha: float
    h_gamma: float
    fst: float | None
    realised_high_trait_patch_fraction: float
    realised_high_trait_mass_mean: float
    realised_high_trait_abundance_mean: float
    potential_high_trait_patch_fraction: float
    potential_high_trait_viable: bool
    realised_high_trait_persists: bool
    tau_trait_potential: int | None
    tau_trait_realised: int | None
    tau_allele_loss: int | None
    tau_H_alpha: int | None
    tau_H_gamma: int | None
    tau_FST: int | None
    events: tuple[FirstPassageEvent, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "replicate_index": self.replicate_index,
            "seed": self.seed,
            "final_q_by_patch": list(self.final_q_by_patch),
            "final_population_by_patch": list(self.final_population_by_patch),
            "final_effective_size_by_patch": list(self.final_effective_size_by_patch),
            "final_p_by_patch": list(self.final_p_by_patch),
            "final_high_trait_abundance_by_patch": list(self.final_high_trait_abundance_by_patch),
            "h_alpha": self.h_alpha,
            "h_gamma": self.h_gamma,
            "fst": self.fst,
            "realised_high_trait_patch_fraction": self.realised_high_trait_patch_fraction,
            "realised_high_trait_mass_mean": self.realised_high_trait_mass_mean,
            "realised_high_trait_abundance_mean": self.realised_high_trait_abundance_mean,
            "potential_high_trait_patch_fraction": self.potential_high_trait_patch_fraction,
            "potential_high_trait_viable": self.potential_high_trait_viable,
            "realised_high_trait_persists": self.realised_high_trait_persists,
            "tau_trait_potential": self.tau_trait_potential,
            "tau_trait_realised": self.tau_trait_realised,
            "tau_allele_loss": self.tau_allele_loss,
            "tau_H_alpha": self.tau_H_alpha,
            "tau_H_gamma": self.tau_H_gamma,
            "tau_FST": self.tau_FST,
            "events": [asdict(event) for event in self.events],
        }


@dataclass(frozen=True)
class CellResult:
    experiment_id: str
    profile: str
    scenario_id: str
    parameters: ParameterCell
    replicates: tuple[ReplicateSummary, ...]
    summary: Mapping[str, object]

    def as_dict(self) -> dict[str, object]:
        return {
            "experiment_id": self.experiment_id,
            "profile": self.profile,
            "scenario_id": self.scenario_id,
            "parameters": self.parameters.as_dict(),
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
        row.update(self.parameters.as_dict())
        row.update(_flatten_mapping(self.summary))
        return row


def _finite_base_parameters(grid_size: int = 31) -> DynamicsParameters:
    """Declared finite-bin, two-kernel, coupled-feedback experiment closure."""
    return replace(
        DynamicsParameters(patch_areas=(1.0,)),
        trait_grid_size=grid_size,
        trait_occupancy_mode="finite_trait_bin_recruitment",
        genotype_trait_recruitment="two_kernel_recruitment",
        inheritance_weight=0.5,
        q_feedback_alpha=0.6,
        q_feedback_beta_trait=0.3,
        q_feedback_gamma_allele=0.1,
    )


def quick_profile() -> ExperimentSpec:
    """Tiny legacy/deterministic profile for tests and examples only."""
    return ExperimentSpec(
        profile=PROFILE_QUICK,
        generations=8,
        replicates=3,
        area_reference_values=(1.0,),
        interaction_feedback_values=(3.5,),
        interaction_barrier_values=(0.4, 0.7),
        base_parameters=replace(DynamicsParameters(patch_areas=(1.0,)), trait_grid_size=21),
    )


def standard_profile() -> ExperimentSpec:
    """Moderate finite-bin profile for interpretable local phase diagrams."""
    return ExperimentSpec(
        profile=PROFILE_STANDARD,
        generations=30,
        replicates=20,
        area_reference_values=(0.8, 1.0, 1.2),
        interaction_feedback_values=(3.0, 4.5, 6.0),
        interaction_barrier_values=(0.35, 0.5, 0.65),
        base_parameters=_finite_base_parameters(),
    )


def full_profile() -> ExperimentSpec:
    """Opt-in finite-bin profile; never execute during normal tests."""
    return ExperimentSpec(
        profile=PROFILE_FULL,
        generations=80,
        replicates=100,
        area_reference_values=(0.6, 0.8, 1.0, 1.2, 1.4),
        interaction_feedback_values=(2.5, 3.5, 4.5, 5.5, 6.5),
        interaction_barrier_values=(0.3, 0.4, 0.5, 0.6, 0.7),
        base_parameters=_finite_base_parameters(grid_size=41),
    )


def scenario_one_large(spec: ExperimentSpec) -> LandscapeScenario:
    return LandscapeScenario(SCENARIO_ONE_LARGE, (spec.total_area,), 0.0)


def scenario_equal_isolated(spec: ExperimentSpec) -> LandscapeScenario:
    area = spec.total_area / spec.patch_count
    return LandscapeScenario(SCENARIO_EQUAL_ISOLATED, tuple(area for _ in range(spec.patch_count)), 0.0)


def scenario_equal_migrating(spec: ExperimentSpec) -> LandscapeScenario:
    return LandscapeScenario(SCENARIO_EQUAL_MIGRATING, scenario_equal_isolated(spec).patch_areas, spec.migration_rate)


def default_scenarios(spec: ExperimentSpec) -> tuple[LandscapeScenario, ...]:
    return (scenario_one_large(spec), scenario_equal_isolated(spec), scenario_equal_migrating(spec))


def parameter_grid(spec: ExperimentSpec) -> tuple[ParameterCell, ...]:
    return tuple(
        ParameterCell(index, *values)
        for index, values in enumerate(product(spec.area_reference_values, spec.interaction_feedback_values, spec.interaction_barrier_values))
    )


def derived_seed(master_seed: int, cell_index: int, replicate_index: int) -> int:
    if cell_index < 0 or replicate_index < 0:
        raise ValueError("cell_index and replicate_index must be nonnegative")
    return (master_seed * 1_000_003 + cell_index * 10_007 + replicate_index * 101 + 17) % (2**31 - 1)


def parameters_for_cell(spec: ExperimentSpec, scenario: LandscapeScenario, cell: ParameterCell, *, seed: int) -> DynamicsParameters:
    return replace(
        spec.base_parameters,
        patch_areas=scenario.patch_areas,
        generations=spec.generations,
        area_reference=cell.area_reference,
        interaction_feedback=cell.interaction_feedback,
        interaction_barrier=cell.interaction_barrier,
        migration_rate=scenario.migration_rate,
        random_seed=seed,
    )


def run_parameter_grid(spec: ExperimentSpec, scenarios: Sequence[LandscapeScenario] | None = None) -> tuple[CellResult, ...]:
    selected = tuple(scenarios) if scenarios is not None else default_scenarios(spec)
    output: list[CellResult] = []
    for scenario in selected:
        for cell in parameter_grid(spec):
            reps = tuple(_run_replicate(spec, scenario, cell, index) for index in range(spec.replicates))
            output.append(CellResult(spec.experiment_id, spec.profile, scenario.scenario_id, cell, reps, _aggregate_replicates(reps)))
    return tuple(output)


def results_to_csv_rows(results: Iterable[CellResult]) -> tuple[dict[str, object], ...]:
    return tuple(result.to_csv_row() for result in results)


def _run_replicate(spec: ExperimentSpec, scenario: LandscapeScenario, cell: ParameterCell, replicate_index: int) -> ReplicateSummary:
    seed = derived_seed(spec.master_seed, cell.cell_index, replicate_index)
    result = simulate(parameters_for_cell(spec, scenario, cell, seed=seed))
    return summarise_replicate(
        result,
        replicate_index=replicate_index,
        seed=seed,
        h_alpha_warning_threshold=spec.h_alpha_warning_threshold,
        h_gamma_warning_threshold=spec.h_gamma_warning_threshold,
        fst_warning_threshold=spec.fst_warning_threshold,
        allele_loss_threshold=spec.allele_loss_threshold,
    )


def summarise_replicate(
    result: SimulationResult,
    *,
    replicate_index: int,
    seed: int,
    h_alpha_warning_threshold: float,
    h_gamma_warning_threshold: float,
    fst_warning_threshold: float,
    allele_loss_threshold: float = 0.0,
) -> ReplicateSummary:
    final = result.snapshots[-1]
    patch_count = len(final.population)
    realised = final.trait_occupancy
    realised_flags = tuple(summary.realised_high_trait_occupied for summary in realised)
    realised_mass = tuple(summary.high_trait_mass for summary in realised)
    realised_abundance = tuple(summary.high_trait_abundance for summary in realised)
    potential_flags = tuple(summary.high_trait_component_present for summary in final.trait_space)
    events = first_passage_events(
        result,
        h_alpha_threshold=h_alpha_warning_threshold,
        h_gamma_threshold=h_gamma_warning_threshold,
        fst_threshold=fst_warning_threshold,
        allele_threshold=allele_loss_threshold,
    )
    times = {event.name: event.time for event in events}
    return ReplicateSummary(
        replicate_index=replicate_index,
        seed=seed,
        final_q_by_patch=final.interaction,
        final_population_by_patch=final.population,
        final_effective_size_by_patch=final.effective_size,
        final_p_by_patch=final.high_allele_frequency,
        final_high_trait_abundance_by_patch=realised_abundance,
        h_alpha=final.h_alpha,
        h_gamma=final.h_gamma,
        fst=final.fst,
        realised_high_trait_patch_fraction=sum(realised_flags) / patch_count,
        realised_high_trait_mass_mean=sum(realised_mass) / patch_count,
        realised_high_trait_abundance_mean=sum(realised_abundance) / patch_count,
        potential_high_trait_patch_fraction=sum(potential_flags) / patch_count,
        potential_high_trait_viable=any(potential_flags),
        realised_high_trait_persists=any(realised_flags),
        tau_trait_potential=times["tau_trait_potential"],
        tau_trait_realised=times["tau_trait_realised"],
        tau_allele_loss=times["tau_allele_loss"],
        tau_H_alpha=times["tau_H_alpha"],
        tau_H_gamma=times["tau_H_gamma"],
        tau_FST=times["tau_FST"],
        events=events,
    )


def _aggregate_replicates(replicates: Sequence[ReplicateSummary]) -> dict[str, object]:
    if not replicates:
        raise ValueError("replicates must be nonempty")
    warnings = ("tau_H_alpha", "tau_H_gamma", "tau_FST", "tau_allele_loss")
    return {
        "uncertainty": "empirical 2.5%, 50%, and 97.5% quantiles across replicates",
        "valid_event_pair_counts": {
            f"{name}_vs_trait_realised": _valid_pair_count(replicates, name, "tau_trait_realised") for name in warnings
        },
        "censored_event_counts": {name: _censored_count(replicates, name) for name in (
            "tau_trait_potential", "tau_trait_realised", "tau_allele_loss", "tau_H_alpha", "tau_H_gamma", "tau_FST"
        )},
        "metrics": {
            "H_alpha": _metric_summary(rep.h_alpha for rep in replicates),
            "H_gamma": _metric_summary(rep.h_gamma for rep in replicates),
            "F_ST": _metric_summary(rep.fst for rep in replicates if rep.fst is not None),
            "realised_high_trait_occupancy": _metric_summary(rep.realised_high_trait_patch_fraction for rep in replicates),
            "realised_high_trait_abundance": _metric_summary(rep.realised_high_trait_abundance_mean for rep in replicates),
            "potential_high_trait_occupancy": _metric_summary(rep.potential_high_trait_patch_fraction for rep in replicates),
        },
        "probabilities": {
            "realised_high_trait_persistence_final": _probability(rep.realised_high_trait_persists for rep in replicates),
            "potential_high_trait_viability_final": _probability(rep.potential_high_trait_viable for rep in replicates),
            "genetic_lead_H_alpha": _lead_probability(replicates, "tau_H_alpha"),
            "genetic_lead_H_gamma": _lead_probability(replicates, "tau_H_gamma"),
            "genetic_lead_FST": _lead_probability(replicates, "tau_FST"),
            "allele_loss_lead": _lead_probability(replicates, "tau_allele_loss"),
        },
        "first_passage": {
            **{name: _event_summary(replicates, name) for name in (
                "tau_trait_potential", "tau_trait_realised", "tau_allele_loss", "tau_H_alpha", "tau_H_gamma", "tau_FST"
            )},
            **{f"{name}_minus_tau_trait_realised": _difference_summary(replicates, name) for name in warnings},
        },
    }


def _metric_summary(values: Iterable[float]) -> dict[str, float | None]:
    observed = tuple(float(value) for value in values)
    if not observed:
        return {"mean": None, "q025": None, "median": None, "q975": None}
    return {"mean": sum(observed) / len(observed), "q025": _quantile(observed, 0.025), "median": _quantile(observed, 0.5), "q975": _quantile(observed, 0.975)}


def _event_summary(replicates: Sequence[ReplicateSummary], attribute: str) -> dict[str, object]:
    values = tuple(getattr(rep, attribute) for rep in replicates)
    observed = tuple(float(value) for value in values if value is not None)
    metadata = tuple(event for rep in replicates for event in rep.events if event.name == attribute)
    return {
        "values": list(values),
        "median": None if not observed else median(observed),
        "q025": None if not observed else _quantile(observed, 0.025),
        "q975": None if not observed else _quantile(observed, 0.975),
        "censored_count": sum(value is None for value in values),
        "threshold": None if not metadata else metadata[0].threshold,
        "aggregation_rule": None if not metadata else metadata[0].aggregation_rule,
    }


def _difference_summary(replicates: Sequence[ReplicateSummary], warning_attribute: str) -> dict[str, object]:
    values = tuple(
        getattr(rep, warning_attribute) - rep.tau_trait_realised
        for rep in replicates
        if getattr(rep, warning_attribute) is not None and rep.tau_trait_realised is not None
    )
    return {
        "values": list(values),
        "median": None if not values else median(values),
        "q025": None if not values else _quantile(values, 0.025),
        "q975": None if not values else _quantile(values, 0.975),
        "censored_count": len(replicates) - len(values),
    }


def _lead_probability(replicates: Sequence[ReplicateSummary], warning_attribute: str) -> float:
    valid = tuple(
        getattr(rep, warning_attribute) < rep.tau_trait_realised
        for rep in replicates
        if getattr(rep, warning_attribute) is not None and rep.tau_trait_realised is not None
    )
    return 0.0 if not valid else sum(valid) / len(valid)


def _valid_pair_count(replicates: Sequence[ReplicateSummary], left: str, right: str) -> int:
    return sum(getattr(rep, left) is not None and getattr(rep, right) is not None for rep in replicates)


def _censored_count(replicates: Sequence[ReplicateSummary], attribute: str) -> int:
    return sum(getattr(rep, attribute) is None for rep in replicates)


def _probability(values: Iterable[bool]) -> float:
    observed = tuple(values)
    return sum(observed) / len(observed)


def _quantile(values: Sequence[float], probability: float) -> float:
    if not 0.0 <= probability <= 1.0:
        raise ValueError("probability must lie in [0, 1]")
    ordered = sorted(float(value) for value in values)
    if len(ordered) == 1:
        return ordered[0]
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] * (1.0 - position + lower) + ordered[upper] * (position - lower)


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
