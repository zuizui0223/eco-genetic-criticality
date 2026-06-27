"""Reproducible fixed-total-capacity phase diagrams for the finite H3 lifecycle."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence
import csv
import json

from causal_model.h3_landscape_presets import (
    LandscapePreset,
    equal_complete_network,
    equal_isolated,
    one_large,
)
from causal_model.network_h3_experiments import H3EnsembleSummary, simulate_h3_ensemble
from causal_model.network_h3_lifecycle import NetworkLifecycleParameters, PatchState


@dataclass(frozen=True)
class H3PhaseDiagramSpec:
    """Predeclared grid and initial composition for matched-capacity H3 runs."""

    total_capacity: int
    patch_count: int
    generations: int = 40
    replicates: int = 100
    adult_survival_grid: tuple[float, ...] = (0.6, 0.8)
    emigration_grid: tuple[float, ...] = (0.0, 0.1, 0.3)
    recruitment_per_adult: float = 0.7
    high_trait_recruitment_multiplier: float = 1.0
    persistence_threshold: int = 2
    colonisation_threshold: int = 2
    initial_occupancy_fraction: float = 0.5
    initial_high_trait_fraction: float = 0.5
    initial_high_allele_frequency: float = 0.5
    random_seed: int = 1

    def __post_init__(self) -> None:
        if self.total_capacity < 1 or self.patch_count < 1:
            raise ValueError("total_capacity and patch_count must be positive")
        if self.total_capacity % self.patch_count:
            raise ValueError("total_capacity must divide exactly across patch_count")
        if self.generations < 1 or self.replicates < 1:
            raise ValueError("generations and replicates must be positive")
        for name, values in (
            ("adult_survival_grid", self.adult_survival_grid),
            ("emigration_grid", self.emigration_grid),
        ):
            if not values or any(not 0.0 <= float(value) <= 1.0 for value in values):
                raise ValueError(f"{name} must be nonempty values in [0, 1]")
        for name, value in (
            ("initial_occupancy_fraction", self.initial_occupancy_fraction),
            ("initial_high_trait_fraction", self.initial_high_trait_fraction),
            ("initial_high_allele_frequency", self.initial_high_allele_frequency),
        ):
            if not 0.0 <= float(value) <= 1.0:
                raise ValueError(f"{name} must lie in [0, 1]")
        if self.recruitment_per_adult < 0.0 or self.high_trait_recruitment_multiplier < 0.0:
            raise ValueError("recruitment parameters must be non-negative")
        if self.persistence_threshold < 1 or self.colonisation_threshold < 1:
            raise ValueError("thresholds must be positive")


@dataclass(frozen=True)
class H3PhaseDiagramRow:
    """One landscape × life-cycle cell, including a full declared denominator."""

    scenario_id: str
    total_capacity: int
    patch_count: int
    patch_capacity: int
    adult_survival_probability: float
    emigration_probability: float
    generations: int
    replicates: int
    persistence_threshold: int
    colonisation_threshold: int
    metapopulation_extinction_probability: float
    realised_high_trait_loss_probability: float
    recolonisation_probability: float
    rescue_probability: float
    final_occupied_patch_count_median: float
    final_high_trait_patch_count_median: float
    final_h_alpha_median: float
    final_h_gamma_median: float
    final_fst_median: float | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def standard_h3_landscapes(spec: H3PhaseDiagramSpec) -> tuple[LandscapePreset, ...]:
    """Return the three predeclared fixed-total-capacity H3 comparisons."""
    return (
        one_large(spec.total_capacity),
        equal_isolated(spec.total_capacity, spec.patch_count),
        equal_complete_network(spec.total_capacity, spec.patch_count),
    )


def initial_states_for_preset(preset: LandscapePreset, spec: H3PhaseDiagramSpec) -> tuple[PatchState, ...]:
    """Create matched-composition starting states proportional to patch capacity."""
    states: list[PatchState] = []
    for capacity in preset.capacities:
        population = round(capacity * spec.initial_occupancy_fraction)
        high_trait = round(population * spec.initial_high_trait_fraction)
        allele_copies = round(2 * population * spec.initial_high_allele_frequency)
        states.append(PatchState(population, high_trait, allele_copies))
    return tuple(states)


def _parameters(
    preset: LandscapePreset,
    spec: H3PhaseDiagramSpec,
    survival: float,
    emigration: float,
    seed: int,
) -> NetworkLifecycleParameters:
    return NetworkLifecycleParameters(
        capacities=preset.capacities,
        source_to_destination_kernel=preset.kernel,
        generations=spec.generations,
        adult_survival_probability=survival,
        emigration_probability=emigration,
        recruitment_per_adult=spec.recruitment_per_adult,
        high_trait_recruitment_multiplier=spec.high_trait_recruitment_multiplier,
        persistence_threshold=spec.persistence_threshold,
        colonisation_threshold=spec.colonisation_threshold,
        random_seed=seed,
    )


def _row(
    preset: LandscapePreset,
    spec: H3PhaseDiagramSpec,
    survival: float,
    emigration: float,
    cell_index: int,
) -> H3PhaseDiagramRow:
    parameters = _parameters(
        preset,
        spec,
        survival,
        emigration,
        spec.random_seed + cell_index * 10_000,
    )
    _, summary = simulate_h3_ensemble(
        parameters,
        initial_states_for_preset(preset, spec),
        replicates=spec.replicates,
    )
    return H3PhaseDiagramRow(
        scenario_id=preset.name,
        total_capacity=preset.total_capacity,
        patch_count=len(preset.capacities),
        patch_capacity=preset.capacities[0],
        adult_survival_probability=survival,
        emigration_probability=emigration,
        generations=spec.generations,
        replicates=spec.replicates,
        persistence_threshold=spec.persistence_threshold,
        colonisation_threshold=spec.colonisation_threshold,
        metapopulation_extinction_probability=summary.metapopulation_extinction_probability,
        realised_high_trait_loss_probability=summary.realised_high_trait_loss_probability,
        recolonisation_probability=summary.recolonisation_probability,
        rescue_probability=summary.rescue_probability,
        final_occupied_patch_count_median=summary.final_occupied_patch_count_median,
        final_high_trait_patch_count_median=summary.final_high_trait_patch_count_median,
        final_h_alpha_median=summary.final_h_alpha_median,
        final_h_gamma_median=summary.final_h_gamma_median,
        final_fst_median=summary.final_fst_median,
    )


def run_h3_phase_diagram(
    spec: H3PhaseDiagramSpec,
    *,
    landscapes: Sequence[LandscapePreset] | None = None,
) -> tuple[H3PhaseDiagramRow, ...]:
    """Run every predeclared H3 landscape and parameter cell reproducibly."""
    declared_landscapes = tuple(standard_h3_landscapes(spec) if landscapes is None else landscapes)
    if not declared_landscapes:
        raise ValueError("at least one landscape is required")
    if any(preset.total_capacity != spec.total_capacity for preset in declared_landscapes):
        raise ValueError("all landscapes must match spec.total_capacity")
    rows = []
    for cell_index, (preset, survival, emigration) in enumerate(
        product(declared_landscapes, spec.adult_survival_grid, spec.emigration_grid)
    ):
        rows.append(_row(preset, spec, float(survival), float(emigration), cell_index))
    return tuple(rows)


def write_h3_phase_diagram_artifacts(
    rows: Iterable[H3PhaseDiagramRow],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    """Write flat CSV and self-describing JSON artifacts from declared rows."""
    values = tuple(rows)
    if not values:
        raise ValueError("rows must be nonempty")
    dictionaries = tuple(row.as_dict() for row in values)
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dictionaries[0]))
        writer.writeheader()
        writer.writerows(dictionaries)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump(dictionaries, handle, ensure_ascii=False, indent=2)
