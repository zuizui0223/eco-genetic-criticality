from causal_model.multipatch_criticality_experiments import quick_profile
from causal_model.mutation_primary_h1_h2_h3_runtime import _certificate_adapter


def test_primary_chain_certificate_adapter_matches_keyword_only_api():
    spec = quick_profile()
    certificate = _certificate_adapter(6.0, 4.0, 1.0, 1.0, spec.base_parameters)
    assert certificate.bifurcation.feedback_strength == 6.0
    assert certificate.bifurcation.area == 4.0
