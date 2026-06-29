from causal_model.branch_locked_h1_h2_h3_chain import _validate_master_seeds


def test_seed_ensemble_validation_accepts_two_distinct_values():
    assert _validate_master_seeds((101, 202)) == (101, 202)
