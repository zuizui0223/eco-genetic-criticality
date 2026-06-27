import pandas as pd
import pytest

from causal_model.theorem_boundary_report import (
    available_metrics,
    load_phase_artifact,
    metric_matrix,
    write_theorem_boundary_report,
)


def _artifact(path):
    pd.DataFrame(
        {
            "scenario_id": ["one_large", "one_large", "one_large", "one_large", "equal_migrating", "equal_migrating"],
            "area_reference": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "interaction_feedback": [3.0, 3.0, 5.0, 5.0, 3.0, 5.0],
            "interaction_barrier": [0.4, 0.6, 0.4, 0.6, 0.4, 0.4],
            "outcomes.realised_high_trait_persistence_final": [1.0, 0.5, 0.5, 0.0, 0.5, 0.25],
            "outcomes.genetic_lead_H_alpha_conditional": [1.0, None, 0.5, 0.0, None, 1.0],
            "scope.maximum_canonical_update_residual.mean": [0.0, 0.0, 0.1, 0.2, 0.3, 0.4],
            "scope.single_patch_canonical_theorem_limit_probability": [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        }
    ).to_csv(path, index=False)


def test_metric_matrix_preserves_conditional_lead_missingness(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    frame = load_phase_artifact(csv_path)
    matrix = metric_matrix(
        frame,
        scenario_id="one_large",
        area_reference=1.0,
        metric="outcomes.genetic_lead_H_alpha_conditional",
    )

    assert matrix.loc[3.0, 0.4] == 1.0
    assert pd.isna(matrix.loc[3.0, 0.6])
    assert matrix.loc[5.0, 0.6] == 0.0


def test_report_writes_one_figure_per_metric_scenario_and_area_reference(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    figures = write_theorem_boundary_report(
        csv_path,
        output_dir=tmp_path / "report",
        metrics=(
            "outcomes.realised_high_trait_persistence_final",
            "outcomes.genetic_lead_H_alpha_conditional",
        ),
        dpi=100,
    )

    assert len(figures) == 4
    assert all(figure.path.exists() for figure in figures)
    report = (tmp_path / "report" / "REPORT.md").read_text(encoding="utf-8")
    assert "conditional genetic-lead metrics" in report
    assert "Missing cells" in report


def test_report_rejects_unknown_metrics_and_missing_required_columns(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    with pytest.raises(ValueError, match="requested metrics"):
        write_theorem_boundary_report(csv_path, output_dir=tmp_path / "report", metrics=("missing",))

    incomplete = tmp_path / "incomplete.csv"
    pd.DataFrame({"scenario_id": ["one_large"]}).to_csv(incomplete, index=False)
    with pytest.raises(ValueError, match="required columns"):
        load_phase_artifact(incomplete)


def test_default_metric_selection_only_returns_columns_present(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    frame = load_phase_artifact(csv_path)

    assert "outcomes.realised_high_trait_persistence_final" in available_metrics(frame)
    assert len(available_metrics(frame)) == 4
