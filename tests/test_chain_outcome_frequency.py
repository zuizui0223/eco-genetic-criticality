from dataclasses import replace

from causal_model.branch_locked_h1_h2_h3_chain import _run_outcome
from causal_model.multipatch_criticality_experiments import (
    default_scenarios,
    parameter_grid,
    quick_profile,
)


def test_outcome_summary_reads_terminal_high_frequency():
    spec = replace(quick_profile(), generations=1)
    scenario = default_scenarios(spec)[0]
    cell = parameter_grid(spec)[0]
    outcome = _run_outcome(
        spec,
        scenario,
        cell,
        seed=101,
        replicate_index=0,
        branch_id="high_start",
        initial_interaction=0.5,
    )
    assert 0.0 <= outcome.terminal_allele_frequency_mean <= 1.0
