"""Canonical-H1-targeted design for finite H1--H3 validation campaigns.

The ordinary experiment profiles use a Cartesian grid of raw interaction
barriers.  Because canonical H1 depends on the local multiplier A/A_ref, the
same raw barrier can be inside the one-large bistable interval for one
(A_ref, kappa) pair and outside it for another.  This module builds an explicit
*design layer* around the existing unified finite validation campaign.

For each (area_reference, interaction_feedback) pair in a profile, it derives
the strict one-large canonical bistable interval and selects declared fractional
positions inside that interval.  Each pair is run as a separate subcampaign
with a deterministic pair-specific seed.  All target positions, raw barriers,
subcampaign specs, raw artifacts, and the aggregate ledger are retained.

The design targets the canonical interaction mechanism; it does not filter on
finite outcomes or on branch-dependent trait viability.  Those results remain
part of the campaign evidence.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from causal_model.canonical_h1_bifurcation import canonical_bistable_barrier_interval
from causal_model.finite_validation_campaign import (
    FiniteValidationCampaign,
    FiniteValidationCampaignPaths,
    run_finite_validation_campaign,
    write_finite_validation_campaign_artifacts,
)
from causal_model.multipatch_criticality_experiments import ExperimentSpec

DEFAULT_INSIDE_POSITIONS = (0.25, 0.5, 0.75)


@dataclass(frozen=True)
class H1TargetedDesignCell:
    """One one-large canonical-bistability target expressed in raw parameters."""

    design_cell_index: int
    pair_index: int
    position_index: int
    area_reference: float
    interaction_feedback: float
    normalized_barrier_position: float
    bistable_barrier_lower: float
    bistable_barrier_upper: float
    interaction_barrier: float
    subcampaign_seed: int

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


@dataclass(frozen=True)
class H1TargetedSubcampaign:
    """One fixed (A_ref, kappa) family of inside-interval barrier positions."""

    pair_index: int
    area_reference: float
    interaction_feedback: float
    interval: tuple[float, float]
    seed: int
    design_cells: tuple[H1TargetedDesignCell, ...]
    spec: ExperimentSpec
    campaign: FiniteValidationCampaign


@dataclass(frozen=True)
class H1TargetedValidationCampaign:
    """Aggregate of all one-large-H1-targeted finite validation subcampaigns."""

    base_spec: ExperimentSpec
    normalized_positions: tuple[float, ...]
    subcampaigns: tuple[H1TargetedSubcampaign, ...]
    ledger: tuple[Mapping[str, object], ...]

    @property
    def design_cell_count(self) -> int:
        return sum(len(subcampaign.design_cells) for subcampaign in self.subcampaigns)


@dataclass(frozen=True)
class H1TargetedCampaignPaths:
    """Top-level outputs plus the directory containing subcampaign artifacts."""

    design_json: Path
    ledger_csv: Path
    ledger_json: Path
    manifest_json: Path
    subcampaign_root: Path

    def as_dict(self) -> dict[str, Path]:
        return {
            "design_json": self.design_json,
            "ledger_csv": self.ledger_csv,
            "ledger_json": self.ledger_json,
            "manifest_json": self.manifest_json,
            "subcampaign_root": self.subcampaign_root,
        }


def build_h1_targeted_design(
    spec: ExperimentSpec,
    *,
    normalized_positions: Sequence[float] = DEFAULT_INSIDE_POSITIONS,
) -> tuple[H1TargetedDesignCell, ...]:
    """Map declared inside-interval positions to raw barriers for one-large H1.

    A target at position r uses theta = lower + r * (upper - lower), with
    0 < r < 1.  Thus every emitted raw barrier lies strictly inside the
    one-large canonical bistable interval for its own (A_ref, kappa) pair.
    """
    positions = _validate_positions(normalized_positions)
    design: list[H1TargetedDesignCell] = []
    design_index = 0
    pair_index = 0
    for area_reference in spec.area_reference_values:
        for interaction_feedback in spec.interaction_feedback_values:
            interval = canonical_bistable_barrier_interval(
                interaction_feedback,
                spec.total_area,
                area_reference,
            )
            if interval is None:
                raise ValueError(
                    "one-large canonical H1 has no strict bistable interval for "
                    f"area_reference={area_reference}, interaction_feedback={interaction_feedback}"
                )
            lower, upper = interval
            seed = _subcampaign_seed(spec.master_seed, pair_index)
            for position_index, position in enumerate(positions):
                design.append(
                    H1TargetedDesignCell(
                        design_cell_index=design_index,
                        pair_index=pair_index,
                        position_index=position_index,
                        area_reference=area_reference,
                        interaction_feedback=interaction_feedback,
                        normalized_barrier_position=position,
                        bistable_barrier_lower=lower,
                        bistable_barrier_upper=upper,
                        interaction_barrier=lower + position * (upper - lower),
                        subcampaign_seed=seed,
                    )
                )
                design_index += 1
            pair_index += 1
    return tuple(design)


def run_h1_targeted_validation_campaign(
    spec: ExperimentSpec,
    *,
    normalized_positions: Sequence[float] = DEFAULT_INSIDE_POSITIONS,
    interaction_separation_threshold: float = 0.05,
    terminal_window: int = 5,
    hysteresis_barrier_points: int = 7,
    hysteresis_barrier_padding: float = 0.1,
    hysteresis_stage_generations: int = 10,
    hysteresis_low_state_threshold: float = 0.25,
    hysteresis_high_state_threshold: float = 0.75,
    coupled_tolerance: float = 1e-12,
) -> H1TargetedValidationCampaign:
    """Run all finite audits over one-large canonical-H1 inside positions.

    No finite outcome is used to choose or remove a design cell. The only design
    criterion is strict canonical bistability for the one-large interaction map.
    """
    design = build_h1_targeted_design(spec, normalized_positions=normalized_positions)
    grouped = _group_design_by_pair(design)
    subcampaigns: list[H1TargetedSubcampaign] = []
    aggregate_ledger: list[Mapping[str, object]] = []
    for pair_index, cells in grouped:
        first = cells[0]
        sub_spec = replace(
            spec,
            area_reference_values=(first.area_reference,),
            interaction_feedback_values=(first.interaction_feedback,),
            interaction_barrier_values=tuple(cell.interaction_barrier for cell in cells),
            master_seed=first.subcampaign_seed,
        )
        campaign = run_finite_validation_campaign(
            sub_spec,
            interaction_separation_threshold=interaction_separation_threshold,
            terminal_window=terminal_window,
            hysteresis_barrier_points=hysteresis_barrier_points,
            hysteresis_barrier_padding=hysteresis_barrier_padding,
            hysteresis_stage_generations=hysteresis_stage_generations,
            hysteresis_low_state_threshold=hysteresis_low_state_threshold,
            hysteresis_high_state_threshold=hysteresis_high_state_threshold,
            coupled_tolerance=coupled_tolerance,
        )
        interval = (first.bistable_barrier_lower, first.bistable_barrier_upper)
        subcampaigns.append(
            H1TargetedSubcampaign(
                pair_index=pair_index,
                area_reference=first.area_reference,
                interaction_feedback=first.interaction_feedback,
                interval=interval,
                seed=first.subcampaign_seed,
                design_cells=cells,
                spec=sub_spec,
                campaign=campaign,
            )
        )
        for design_cell, ledger_row in zip(cells, campaign.ledger, strict=True):
            aggregate_ledger.append(
                {
                    "design": design_cell.as_dict(),
                    "finite_validation": ledger_row,
                }
            )
    return H1TargetedValidationCampaign(
        base_spec=spec,
        normalized_positions=tuple(normalized_positions),
        subcampaigns=tuple(subcampaigns),
        ledger=tuple(aggregate_ledger),
    )


def targeted_output_paths(output_dir: str | Path, prefix: str) -> H1TargetedCampaignPaths:
    """Return deterministic top-level destinations for an H1-targeted campaign."""
    directory = Path(output_dir)
    return H1TargetedCampaignPaths(
        design_json=directory / f"{prefix}.design.json",
        ledger_csv=directory / f"{prefix}.ledger.csv",
        ledger_json=directory / f"{prefix}.ledger.json",
        manifest_json=directory / f"{prefix}.manifest.json",
        subcampaign_root=directory / f"{prefix}.subcampaigns",
    )


def write_h1_targeted_validation_campaign(
    campaign: H1TargetedValidationCampaign,
    *,
    output_dir: str | Path,
    prefix: str,
    audit_arguments: Mapping[str, object],
    code_revision: str | None = None,
) -> H1TargetedCampaignPaths:
    """Write all subcampaign artifacts, design map, aggregate ledger, and manifest."""
    paths = targeted_output_paths(output_dir, prefix)
    paths.subcampaign_root.mkdir(parents=True, exist_ok=True)
    subcampaign_outputs: dict[str, Mapping[str, str]] = {}
    for subcampaign in campaign.subcampaigns:
        name = f"pair_{subcampaign.pair_index:03d}"
        output = write_finite_validation_campaign_artifacts(
            subcampaign.campaign,
            output_dir=paths.subcampaign_root / name,
            prefix=name,
        )
        subcampaign_outputs[name] = {key: str(value) for key, value in output.as_dict().items()}
    _write_json(
        paths.design_json,
        {
            "base_spec": asdict(campaign.base_spec),
            "normalized_positions": list(campaign.normalized_positions),
            "design_cells": [record["design"] for record in campaign.ledger],
        },
    )
    _write_ledger(campaign.ledger, csv_path=paths.ledger_csv, json_path=paths.ledger_json)
    _write_json(
        paths.manifest_json,
        {
            "runner": "causal_model.h1_targeted_validation_campaign",
            "campaign": "one_large_canonical_h1_targeted_v1",
            "code_revision": code_revision,
            "base_spec": asdict(campaign.base_spec),
            "normalized_positions": list(campaign.normalized_positions),
            "subcampaign_count": len(campaign.subcampaigns),
            "design_cell_count": campaign.design_cell_count,
            "audit_arguments": dict(audit_arguments),
            "selection_policy": (
                "Design cells are selected only by strict one-large canonical H1 bistability at declared normalized barrier "
                "positions. No finite branch outcome, trait outcome, warning order, censoring status, or spatial result is used "
                "to exclude a cell."
            ),
            "subcampaigns": [
                {
                    "pair_index": sub.pair_index,
                    "area_reference": sub.area_reference,
                    "interaction_feedback": sub.interaction_feedback,
                    "bistable_interval": list(sub.interval),
                    "seed": sub.seed,
                    "design_cells": [cell.as_dict() for cell in sub.design_cells],
                    "spec": asdict(sub.spec),
                    "outputs": subcampaign_outputs[f"pair_{sub.pair_index:03d}"],
                }
                for sub in campaign.subcampaigns
            ],
            "outputs": {key: str(value) for key, value in paths.as_dict().items()},
        },
    )
    return paths


def _validate_positions(values: Sequence[float]) -> tuple[float, ...]:
    positions = tuple(float(value) for value in values)
    if not positions:
        raise ValueError("normalized_positions must be nonempty")
    if any(not 0.0 < value < 1.0 for value in positions):
        raise ValueError("normalized barrier positions must lie strictly inside (0, 1)")
    if len(set(positions)) != len(positions):
        raise ValueError("normalized barrier positions must be unique")
    return positions


def _group_design_by_pair(
    design: Sequence[H1TargetedDesignCell],
) -> tuple[tuple[int, tuple[H1TargetedDesignCell, ...]], ...]:
    grouped: dict[int, list[H1TargetedDesignCell]] = {}
    for cell in design:
        grouped.setdefault(cell.pair_index, []).append(cell)
    return tuple((pair_index, tuple(cells)) for pair_index, cells in sorted(grouped.items()))


def _subcampaign_seed(master_seed: int, pair_index: int) -> int:
    if pair_index < 0:
        raise ValueError("pair_index must be nonnegative")
    return (master_seed * 1_000_003 + pair_index * 10_007 + 911) % (2**31 - 1)


def _write_ledger(
    ledger: Iterable[Mapping[str, object]],
    *,
    csv_path: Path,
    json_path: Path,
) -> None:
    values = tuple(ledger)
    if not values:
        raise ValueError("ledger must be nonempty")
    rows = [_flatten_mapping(value) for value in values]
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    _write_json(json_path, list(values))


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
