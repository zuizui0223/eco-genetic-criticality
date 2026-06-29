"""Unified reproducible validation campaign for finite H1--H3 audits.

This module intentionally orchestrates existing audits rather than adding a new
mechanistic closure.  Every audit receives the same immutable ExperimentSpec,
the same declared landscapes, and therefore the same parameter grid and seed
rule.  It writes all cells, including canonical-H1-inapplicable cells, finite
precondition failures, censored first-passage pairs, and counterexamples.

The resulting ledger is a navigation layer across five complementary Type S
checks:

* finite H1 branch separation;
* finite H1 stateful hysteresis continuation;
* H2 warning ordering conditional on finite H1 branches;
* the H1--H2--H3 coupled-chain audit; and
* branch-aware H3 fragmentation and migration modulation.

It is not a meta-theorem.  The campaign keeps the theorem scope and finite
simulation scope explicit in each underlying artifact.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

from causal_model.finite_coupled_chain_audit import (
    FiniteCoupledChainCell,
    run_finite_coupled_chain_audit,
    write_finite_coupled_chain_artifacts,
)
from causal_model.finite_h1_branch_separation_audit import (
    FiniteH1BranchSeparationCell,
    run_finite_h1_branch_separation_audit,
    write_finite_h1_branch_separation_artifacts,
)
from causal_model.finite_h1_hysteresis_audit import (
    FiniteH1HysteresisCell,
    run_finite_h1_hysteresis_audit,
    write_finite_h1_hysteresis_artifacts,
)
from causal_model.finite_h2_branch_warning_audit import (
    FiniteH2BranchWarningCell,
    run_finite_h2_branch_warning_audit,
    write_finite_h2_branch_warning_artifacts,
)
from causal_model.h3_branch_aware_fragmentation_audit import (
    H3BranchAwareCell,
    run_h3_branch_aware_fragmentation_audit,
    write_h3_branch_aware_fragmentation_artifacts,
)
from causal_model.multipatch_criticality_experiments import (
    ExperimentSpec,
    ParameterCell,
    default_scenarios,
    parameter_grid,
)


@dataclass(frozen=True)
class FiniteValidationCampaign:
    """All finite H1--H3 audit cells generated from one declared experiment spec."""

    spec: ExperimentSpec
    scenario_ids: tuple[str, ...]
    h1_branch_cells: tuple[FiniteH1BranchSeparationCell, ...]
    h1_hysteresis_cells: tuple[FiniteH1HysteresisCell, ...]
    h2_warning_cells: tuple[FiniteH2BranchWarningCell, ...]
    coupled_chain_cells: tuple[FiniteCoupledChainCell, ...]
    h3_branch_aware_cells: tuple[H3BranchAwareCell, ...]
    ledger: tuple[Mapping[str, object], ...]

    @property
    def cell_count(self) -> int:
        return len(self.ledger)


@dataclass(frozen=True)
class FiniteValidationCampaignPaths:
    """All deterministic output paths produced by a campaign writer."""

    h1_branch_csv: Path
    h1_branch_json: Path
    h1_hysteresis_csv: Path
    h1_hysteresis_json: Path
    h2_warning_csv: Path
    h2_warning_json: Path
    coupled_chain_csv: Path
    coupled_chain_json: Path
    h3_branch_aware_csv: Path
    h3_branch_aware_json: Path
    ledger_csv: Path
    ledger_json: Path

    def as_dict(self) -> dict[str, Path]:
        return {
            "h1_branch_csv": self.h1_branch_csv,
            "h1_branch_json": self.h1_branch_json,
            "h1_hysteresis_csv": self.h1_hysteresis_csv,
            "h1_hysteresis_json": self.h1_hysteresis_json,
            "h2_warning_csv": self.h2_warning_csv,
            "h2_warning_json": self.h2_warning_json,
            "coupled_chain_csv": self.coupled_chain_csv,
            "coupled_chain_json": self.coupled_chain_json,
            "h3_branch_aware_csv": self.h3_branch_aware_csv,
            "h3_branch_aware_json": self.h3_branch_aware_json,
            "ledger_csv": self.ledger_csv,
            "ledger_json": self.ledger_json,
        }


def run_finite_validation_campaign(
    spec: ExperimentSpec,
    *,
    interaction_separation_threshold: float = 0.05,
    terminal_window: int = 5,
    hysteresis_barrier_points: int = 7,
    hysteresis_barrier_padding: float = 0.1,
    hysteresis_stage_generations: int = 10,
    hysteresis_low_state_threshold: float = 0.25,
    hysteresis_high_state_threshold: float = 0.75,
    coupled_tolerance: float = 1e-12,
) -> FiniteValidationCampaign:
    """Run every finite audit on one shared parameter grid and scenario set.

    The shared spec is deliberate.  The campaign is designed to distinguish
    genuine scope differences across audits from accidental differences in
    parameter grids, replicate counts, or seed rules.
    """
    scenarios = default_scenarios(spec)
    h1_branch_cells = run_finite_h1_branch_separation_audit(
        spec,
        scenarios=scenarios,
        interaction_separation_threshold=interaction_separation_threshold,
        terminal_window=terminal_window,
    )
    h1_hysteresis_cells = run_finite_h1_hysteresis_audit(
        spec,
        scenarios=scenarios,
        barrier_points=hysteresis_barrier_points,
        barrier_padding=hysteresis_barrier_padding,
        stage_generations=hysteresis_stage_generations,
        interaction_separation_threshold=interaction_separation_threshold,
        low_state_threshold=hysteresis_low_state_threshold,
        high_state_threshold=hysteresis_high_state_threshold,
    )
    h2_warning_cells = run_finite_h2_branch_warning_audit(
        spec,
        scenarios=scenarios,
        interaction_separation_threshold=interaction_separation_threshold,
        terminal_window=terminal_window,
    )
    coupled_chain_cells = run_finite_coupled_chain_audit(
        spec,
        scenarios=scenarios,
        tolerance=coupled_tolerance,
    )
    h3_branch_aware_cells = run_h3_branch_aware_fragmentation_audit(
        spec,
        scenarios=scenarios,
        interaction_separation_threshold=interaction_separation_threshold,
    )
    ledger = _build_ledger(
        spec,
        h1_branch_cells=h1_branch_cells,
        h1_hysteresis_cells=h1_hysteresis_cells,
        h2_warning_cells=h2_warning_cells,
        coupled_chain_cells=coupled_chain_cells,
        h3_branch_aware_cells=h3_branch_aware_cells,
    )
    return FiniteValidationCampaign(
        spec=spec,
        scenario_ids=tuple(scenario.scenario_id for scenario in scenarios),
        h1_branch_cells=h1_branch_cells,
        h1_hysteresis_cells=h1_hysteresis_cells,
        h2_warning_cells=h2_warning_cells,
        coupled_chain_cells=coupled_chain_cells,
        h3_branch_aware_cells=h3_branch_aware_cells,
        ledger=ledger,
    )


def campaign_output_paths(output_dir: str | Path, prefix: str) -> FiniteValidationCampaignPaths:
    """Return all campaign destinations without creating or replacing files."""
    directory = Path(output_dir)
    return FiniteValidationCampaignPaths(
        h1_branch_csv=directory / f"{prefix}.h1_branch.csv",
        h1_branch_json=directory / f"{prefix}.h1_branch.json",
        h1_hysteresis_csv=directory / f"{prefix}.h1_hysteresis.csv",
        h1_hysteresis_json=directory / f"{prefix}.h1_hysteresis.json",
        h2_warning_csv=directory / f"{prefix}.h2_warning.csv",
        h2_warning_json=directory / f"{prefix}.h2_warning.json",
        coupled_chain_csv=directory / f"{prefix}.coupled_chain.csv",
        coupled_chain_json=directory / f"{prefix}.coupled_chain.json",
        h3_branch_aware_csv=directory / f"{prefix}.h3_branch_aware.csv",
        h3_branch_aware_json=directory / f"{prefix}.h3_branch_aware.json",
        ledger_csv=directory / f"{prefix}.ledger.csv",
        ledger_json=directory / f"{prefix}.ledger.json",
    )


def write_finite_validation_campaign_artifacts(
    campaign: FiniteValidationCampaign,
    *,
    output_dir: str | Path,
    prefix: str,
) -> FiniteValidationCampaignPaths:
    """Write every audit artifact plus a cell-aligned cross-audit ledger."""
    paths = campaign_output_paths(output_dir, prefix)
    write_finite_h1_branch_separation_artifacts(
        campaign.h1_branch_cells,
        csv_path=paths.h1_branch_csv,
        json_path=paths.h1_branch_json,
    )
    write_finite_h1_hysteresis_artifacts(
        campaign.h1_hysteresis_cells,
        csv_path=paths.h1_hysteresis_csv,
        json_path=paths.h1_hysteresis_json,
    )
    write_finite_h2_branch_warning_artifacts(
        campaign.h2_warning_cells,
        csv_path=paths.h2_warning_csv,
        json_path=paths.h2_warning_json,
    )
    write_finite_coupled_chain_artifacts(
        campaign.coupled_chain_cells,
        csv_path=paths.coupled_chain_csv,
        json_path=paths.coupled_chain_json,
    )
    write_h3_branch_aware_fragmentation_artifacts(
        campaign.h3_branch_aware_cells,
        csv_path=paths.h3_branch_aware_csv,
        json_path=paths.h3_branch_aware_json,
    )
    _write_ledger(campaign.ledger, csv_path=paths.ledger_csv, json_path=paths.ledger_json)
    return paths


def _build_ledger(
    spec: ExperimentSpec,
    *,
    h1_branch_cells: Iterable[FiniteH1BranchSeparationCell],
    h1_hysteresis_cells: Iterable[FiniteH1HysteresisCell],
    h2_warning_cells: Iterable[FiniteH2BranchWarningCell],
    coupled_chain_cells: Iterable[FiniteCoupledChainCell],
    h3_branch_aware_cells: Iterable[H3BranchAwareCell],
) -> tuple[Mapping[str, object], ...]:
    branch_by_key = _scenario_summary_map(h1_branch_cells)
    hysteresis_by_key = _scenario_summary_map(h1_hysteresis_cells)
    h2_by_key = _scenario_summary_map(h2_warning_cells)
    chain_by_key = _summary_map(coupled_chain_cells)
    h3_by_key = _summary_map(h3_branch_aware_cells)
    rows: list[Mapping[str, object]] = []
    for parameters in parameter_grid(spec):
        key = _parameter_key(parameters)
        h1 = branch_by_key.get(key, {})
        hysteresis = hysteresis_by_key.get(key, {})
        h2 = h2_by_key.get(key, {})
        if key not in chain_by_key or key not in h3_by_key:
            raise RuntimeError("campaign audit outputs do not cover the declared parameter grid")
        rows.append(
            {
                "parameters": parameters.as_dict(),
                "audit_coverage": {
                    "h1_branch_scenario_count": len(h1),
                    "h1_hysteresis_scenario_count": len(hysteresis),
                    "h2_warning_scenario_count": len(h2),
                    "coupled_chain_available": True,
                    "h3_branch_aware_available": True,
                },
                "h1_branch_by_scenario": h1,
                "h1_hysteresis_by_scenario": hysteresis,
                "h2_warning_by_scenario": h2,
                "coupled_chain": chain_by_key[key],
                "h3_branch_aware": h3_by_key[key],
            }
        )
    return tuple(rows)


def _scenario_summary_map(cells: Iterable[object]) -> dict[tuple[int, float, float, float], dict[str, Mapping[str, object]]]:
    result: dict[tuple[int, float, float, float], dict[str, Mapping[str, object]]] = {}
    for cell in cells:
        key = _parameter_key(cell.parameters)
        scenario_map = result.setdefault(key, {})
        if cell.scenario_id in scenario_map:
            raise RuntimeError("duplicate scenario summary for one campaign parameter cell")
        scenario_map[cell.scenario_id] = dict(cell.summary)
    return result


def _summary_map(cells: Iterable[object]) -> dict[tuple[int, float, float, float], Mapping[str, object]]:
    result: dict[tuple[int, float, float, float], Mapping[str, object]] = {}
    for cell in cells:
        key = _parameter_key(cell.parameters)
        if key in result:
            raise RuntimeError("duplicate campaign summary for one parameter cell")
        result[key] = dict(cell.summary)
    return result


def _parameter_key(parameters: ParameterCell) -> tuple[int, float, float, float]:
    return (
        parameters.cell_index,
        parameters.area_reference,
        parameters.interaction_feedback,
        parameters.interaction_barrier,
    )


def _write_ledger(
    ledger: Iterable[Mapping[str, object]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> None:
    values = tuple(ledger)
    if not values:
        raise ValueError("ledger must be nonempty")
    csv_target = Path(csv_path)
    json_target = Path(json_path)
    csv_target.parent.mkdir(parents=True, exist_ok=True)
    json_target.parent.mkdir(parents=True, exist_ok=True)
    rows = [_flatten_mapping(row) for row in values]
    fieldnames = sorted({key for row in rows for key in row})
    with csv_target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with json_target.open("w", encoding="utf-8") as handle:
        json.dump(list(values), handle, indent=2, sort_keys=True)


def _flatten_mapping(mapping: Mapping[str, object], prefix: str = "") -> dict[str, object]:
    flat: dict[str, object] = {}
    for key, value in mapping.items():
        name = f"{prefix}{key}" if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(_flatten_mapping(value, name))
        else:
            flat[name] = value
    return flat
