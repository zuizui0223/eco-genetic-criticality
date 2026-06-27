from types import SimpleNamespace

import pytest

from causal_model.censoring_aware_phase_diagram import (
    censoring_aware_cell_row,
    event_comparison_for_replicates,
)


def _replicate(**times):
    defaults = {
        "tau_H_alpha": None,
        "tau_H_gamma": None,
        "tau_FST": None,
        "tau_allele_loss": None,
        "tau_trait_realised": None,
    }
    defaults.update(times)
    return SimpleNamespace(**defaults)


def test_event_comparison_reports_pair_observability_and_both_probabilities():
    replicates = (
        _replicate(tau_H_alpha=1, tau_trait_realised=3),
        _replicate(tau_H_alpha=5, tau_trait_realised=4),
        _replicate(tau_H_alpha=None, tau_trait_realised=2),
    )
    summary = event_comparison_for_replicates(replicates, "tau_H_alpha")

    assert summary["valid_pair_count"] == 2
    assert summary["censored_pair_count"] == 1
    assert summary["valid_pair_probability"] == pytest.approx(2 / 3)
    assert summary["conditional_lead_probability"] == pytest.approx(0.5)
    assert summary["unconditional_observed_lead_fraction"] == pytest.approx(1 / 3)


def test_cell_row_flattens_all_warning_comparisons_without_terminal_imputation():
    replicates = (
        _replicate(tau_H_alpha=1, tau_H_gamma=None, tau_FST=4, tau_allele_loss=2, tau_trait_realised=3),
        _replicate(tau_H_alpha=None, tau_H_gamma=2, tau_FST=None, tau_allele_loss=5, tau_trait_realised=4),
    )
    parameters = SimpleNamespace(as_dict=lambda: {"area_reference": 1.0})
    cell = SimpleNamespace(
        experiment_id="test",
        profile="test",
        scenario_id="test",
        replicates=replicates,
        parameters=parameters,
    )

    row = censoring_aware_cell_row(cell)
    alpha_prefix = "tau_H_alpha_vs_tau_trait_realised"
    assert row["replicate_count"] == 2
    assert row["area_reference"] == 1.0
    assert row[f"{alpha_prefix}.valid_pair_count"] == 1
    assert row[f"{alpha_prefix}.conditional_lead_probability"] == pytest.approx(1.0)
    assert row[f"{alpha_prefix}.unconditional_observed_lead_fraction"] == pytest.approx(0.5)
    assert row[f"{alpha_prefix}.time_differences"] == "-2"
