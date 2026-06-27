"""Turn theorem-boundary CSV artifacts into figure-ready robustness reports.

The report layer never reruns the simulation.  It only reads the flat artifact
written by ``theorem_boundary_phase_diagram`` and produces one heat map per
scenario, area-reference value, and metric.  This separation keeps simulation,
statistical aggregation, and visualisation independently reproducible.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PARAMETER_COLUMNS = (
    "scenario_id",
    "area_reference",
    "interaction_feedback",
    "interaction_barrier",
)
DEFAULT_METRICS = (
    "outcomes.realised_high_trait_persistence_final",
    "outcomes.genetic_lead_H_alpha_conditional",
    "scope.maximum_canonical_update_residual.mean",
    "scope.single_patch_canonical_theorem_limit_probability",
)
METRIC_TITLES = {
    "outcomes.realised_high_trait_persistence_final": "Realised high-trait persistence probability",
    "outcomes.genetic_lead_H_alpha_conditional": "Conditional Hα genetic-lead probability",
    "scope.maximum_canonical_update_residual.mean": "Mean maximum canonical-update residual",
    "scope.single_patch_canonical_theorem_limit_probability": "Canonical H1 theorem-limit probability",
}


@dataclass(frozen=True)
class ReportFigure:
    """Metadata for a single saved heat map."""

    metric: str
    scenario_id: str
    area_reference: float
    path: Path
    observed_cells: int
    missing_cells: int


def load_phase_artifact(csv_path: str | Path) -> pd.DataFrame:
    """Load and validate a flat theorem-boundary phase-diagram CSV artifact."""
    frame = pd.read_csv(csv_path)
    required = set(PARAMETER_COLUMNS)
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"artifact is missing required columns: {', '.join(missing)}")
    if frame.empty:
        raise ValueError("artifact has no phase-diagram rows")
    for name in ("area_reference", "interaction_feedback", "interaction_barrier"):
        frame[name] = pd.to_numeric(frame[name], errors="raise")
    return frame


def available_metrics(frame: pd.DataFrame) -> tuple[str, ...]:
    """Return default report metrics available in this artifact."""
    return tuple(metric for metric in DEFAULT_METRICS if metric in frame.columns)


def metric_matrix(
    frame: pd.DataFrame,
    *,
    scenario_id: str,
    area_reference: float,
    metric: str,
) -> pd.DataFrame:
    """Build interaction-feedback × barrier matrix, retaining missing values as NA."""
    if metric not in frame.columns:
        raise ValueError(f"metric not found in artifact: {metric}")
    subset = frame.loc[
        (frame["scenario_id"] == scenario_id) & (frame["area_reference"] == area_reference),
        ["interaction_feedback", "interaction_barrier", metric],
    ].copy()
    if subset.empty:
        raise ValueError("no rows match scenario_id and area_reference")
    subset[metric] = pd.to_numeric(subset[metric], errors="coerce")
    if subset.duplicated(["interaction_feedback", "interaction_barrier"]).any():
        raise ValueError("artifact has duplicate parameter cells for this scenario and area reference")
    return subset.pivot(index="interaction_feedback", columns="interaction_barrier", values=metric).sort_index().sort_index(axis=1)


def write_theorem_boundary_report(
    csv_path: str | Path,
    *,
    output_dir: str | Path,
    metrics: Sequence[str] | None = None,
    dpi: int = 200,
) -> tuple[ReportFigure, ...]:
    """Write one annotated heat map per scenario, area reference, and metric."""
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    frame = load_phase_artifact(csv_path)
    selected_metrics = tuple(metrics) if metrics is not None else available_metrics(frame)
    if not selected_metrics:
        raise ValueError("none of the requested report metrics occur in the artifact")
    unknown = sorted(set(selected_metrics).difference(frame.columns))
    if unknown:
        raise ValueError(f"requested metrics not found in artifact: {', '.join(unknown)}")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    figures: list[ReportFigure] = []
    for scenario_id in sorted(frame["scenario_id"].dropna().unique()):
        scenario = frame.loc[frame["scenario_id"] == scenario_id]
        for area_reference in sorted(scenario["area_reference"].unique()):
            for metric in selected_metrics:
                matrix = metric_matrix(
                    frame,
                    scenario_id=str(scenario_id),
                    area_reference=float(area_reference),
                    metric=metric,
                )
                filename = _safe_filename(metric, scenario_id=str(scenario_id), area_reference=float(area_reference))
                path = destination / filename
                _write_heatmap(
                    matrix,
                    path=path,
                    title=_title(metric, scenario_id=str(scenario_id), area_reference=float(area_reference)),
                    dpi=dpi,
                )
                figures.append(
                    ReportFigure(
                        metric=metric,
                        scenario_id=str(scenario_id),
                        area_reference=float(area_reference),
                        path=path,
                        observed_cells=int(matrix.notna().sum().sum()),
                        missing_cells=int(matrix.isna().sum().sum()),
                    )
                )
    _write_index(destination, figures=figures, source=Path(csv_path))
    return tuple(figures)


def _write_heatmap(matrix: pd.DataFrame, *, path: Path, title: str, dpi: int) -> None:
    values = matrix.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)
    figure, axis = plt.subplots(figsize=(7, 5))
    image = axis.imshow(masked, aspect="auto")
    axis.set_title(title)
    axis.set_xlabel("Interaction barrier")
    axis.set_ylabel("Interaction feedback")
    axis.set_xticks(np.arange(len(matrix.columns)), labels=[f"{value:g}" for value in matrix.columns])
    axis.set_yticks(np.arange(len(matrix.index)), labels=[f"{value:g}" for value in matrix.index])
    for row, feedback in enumerate(matrix.index):
        for column, barrier in enumerate(matrix.columns):
            value = matrix.loc[feedback, barrier]
            label = "NA" if pd.isna(value) else f"{value:.3g}"
            axis.text(column, row, label, ha="center", va="center")
    figure.colorbar(image, ax=axis, label="Metric value")
    figure.tight_layout()
    figure.savefig(path, dpi=dpi)
    plt.close(figure)


def _title(metric: str, *, scenario_id: str, area_reference: float) -> str:
    return f"{METRIC_TITLES.get(metric, metric)}\nscenario={scenario_id}; area reference={area_reference:g}"


def _safe_filename(metric: str, *, scenario_id: str, area_reference: float) -> str:
    safe_metric = metric.replace(".", "_")
    safe_scenario = scenario_id.replace("/", "_")
    return f"{safe_metric}__scenario-{safe_scenario}__area-reference-{area_reference:g}.png"


def _write_index(destination: Path, *, figures: Iterable[ReportFigure], source: Path) -> None:
    rows = ["# Theorem-boundary phase-diagram report", "", f"Source CSV: `{source}`", ""]
    rows.extend([
        "Each figure displays one outcome or theorem-scope metric over the interaction-feedback × barrier grid.",
        "Cells marked `NA` lack a valid value; for conditional genetic-lead metrics this usually means no uncensored event pair, not a lead probability of zero.",
        "",
        "| Metric | Scenario | Area reference | Observed cells | Missing cells | Figure |",
        "|---|---:|---:|---:|---:|---|",
    ])
    for figure in figures:
        rows.append(
            f"| `{figure.metric}` | {figure.scenario_id} | {figure.area_reference:g} | {figure.observed_cells} | {figure.missing_cells} | [{figure.path.name}]({figure.path.name}) |"
        )
    (destination / "REPORT.md").write_text("\n".join(rows) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Generate figure-ready reports from theorem-boundary phase-diagram CSV artifacts.")
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--metric", action="append", dest="metrics", help="Metric column to plot. Repeat for several metrics.")
    parser.add_argument("--dpi", type=int, default=200)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    figures = write_theorem_boundary_report(
        args.csv_path,
        output_dir=args.output_dir,
        metrics=args.metrics,
        dpi=args.dpi,
    )
    print(f"Wrote {len(figures)} figures and REPORT.md to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
