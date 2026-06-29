import pandas as pd
import pytest

from causal_model.paired_baseline_report import (
    available_metrics,
    load_paired_baseline_artifact,
    metric_matrix,
    write_paired_baseline_report,
)


def _artifact(path):
    pd.DataFrame(
        {
            "scenario_id": [
                "one_large",
                "one_large",
                "one_large",
                "one_large",
                "equal_migrating",
                "equal_migrating",
            ],
            "area_reference": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            "interaction_feedback": [3.0, 3.0, 5.0, 5.0, 3.0, 5.0],
            "interaction_barrier": [0.4, 0.6, 0.4, 0.6, 0.4, 0.4],
            "paired_contrasts.full_minus_trait_only.realised_high_trait_mass_difference_mean": [
                0.4,
                None,
                0.1,
                -0.2,
                0.3,
                0.0,
            ],
            "paired_contrasts.full_minus_genetic_only.realised_high_trait_mass_difference_mean": [
                0.2,
                0.1,
                0.0,
                -0.1,
                0.4,
                0.3,
            ],
            "paired_contrasts.full_minus_trait_only.final_h_alpha_difference_mean": [
                0.05,
                0.02,
                0.01,
                -0.03,
                0.1,
                0.07,
            ],
            "scope.full_eco_genetic.single_patch_canonical_theorem_limit_probability": [
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        }
    ).to_csv(path, index=False)


def test_metric_matrix_preserves_missing_paired_contrast_values(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    frame = load_paired_baseline_artifact(csv_path)
    matrix = metric_matrix(
        frame,
        scenario_id="one_large",
        area_reference=1.0,
        metric="paired_contrasts.full_minus_trait_only.realised_high_trait_mass_difference_mean",
    )

    assert matrix.loc[3.0, 0.4] == pytest.approx(0.4)
    assert pd.isna(matrix.loc[3.0, 0.6])
    assert matrix.loc[5.0, 0.6] == pytest.approx(-0.2)


def test_report_writes_one_figure_per_metric_scenario_and_area_reference(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    figures = write_paired_baseline_report(
        csv_path,
        output_dir=tmp_path / "report",
        metrics=(
            "paired_contrasts.full_minus_trait_only.realised_high_trait_mass_difference_mean",
            "scope.full_eco_genetic.single_patch_canonical_theorem_limit_probability",
        ),
        dpi=100,
    )

    assert len(figures) == 4
    assert all(figure.path.exists() for figure in figures)
    report = (tmp_path / "report" / "REPORT.md").read_text(encoding="utf-8")
    assert "positive values mean the full eco-genetic model is larger" in report
    assert "theorem limit" in report
    assert "Missing cells" in report


def test_report_rejects_unknown_metrics_and_missing_required_columns(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    with pytest.raises(ValueError, match="requested metrics"):
        write_paired_baseline_report(csv_path, output_dir=tmp_path / "report", metrics=("missing",))

    incomplete = tmp_path / "incomplete.csv"
    pd.DataFrame({"scenario_id": ["one_large"]}).to_csv(incomplete, index=False)
    with pytest.raises(ValueError, match="required columns"):
        load_paired_baseline_artifact(incomplete)


def test_default_metric_selection_only_returns_columns_present(tmp_path):
    csv_path = tmp_path / "artifact.csv"
    _artifact(csv_path)
    frame = load_paired_baseline_artifact(csv_path)

    metrics = available_metrics(frame)
    assert "paired_contrasts.full_minus_trait_only.realised_high_trait_mass_difference_mean" in metrics
    assert "scope.full_eco_genetic.single_patch_canonical_theorem_limit_probability" in metrics
    assert "scope.trait_only.single_patch_canonical_theorem_limit_probability" not in metrics
