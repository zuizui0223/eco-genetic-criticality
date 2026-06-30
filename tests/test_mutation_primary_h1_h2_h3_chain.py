from causal_model.mutation_h1_primary_domain import primary_analysis_cells
from causal_model.mutation_primary_h1_h2_h3_chain import _h2


def test_h2_requires_baseline_eligibility_and_both_observed_events():
    assert _h2("H_alpha", 4, 8, True).warning_leads is True
    censored = _h2("H_alpha", 4, None, True)
    assert censored.valid_pair is False
    assert censored.censored is True
    assert censored.warning_leads is None
    baseline_ineligible = _h2("H_gamma", 1, 8, False)
    assert baseline_ineligible.valid_pair is False
    assert baseline_ineligible.warning_leads is None


def test_primary_chain_uses_frozen_domain_not_all_validation_cells():
    assert len(primary_analysis_cells()) == 12
