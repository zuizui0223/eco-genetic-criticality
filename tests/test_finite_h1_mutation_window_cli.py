from causal_model.finite_h1_mutation_window_audit_cli import build_parser


def test_campaign_role_defaults_to_screen_and_accepts_independent_validation_label():
    assert build_parser().parse_args([]).campaign_role == "screen_v1"
    assert (
        build_parser().parse_args(["--campaign-role", "independent_validation_v1"]).campaign_role
        == "independent_validation_v1"
    )
