import pytest

from causal_model.colonization_metapopulation_abm import ColonizationParameters, ColonizationRegime
from causal_model.colonization_recruitment_factorization import (
    ColonizationRecruitmentContext,
    one_step_recruitment_factors,
    require_theorem_interior,
)


def _params(*, extinction_rate: float = 0.0) -> ColonizationParameters:
    return ColonizationParameters(
        dispersal_cost=0.2,
        fecundity=0.45,
        extinction_rate=extinction_rate,
        density_threshold=0.8,
        mutation_rate=0.01,
        mutation_std=0.05,
        distance_decay=1.0,
        resource_replenishment=0.1,
        base_survival=0.9,
        max_age=8,
        benefit_saturation=0.0,
    )


def _context(*, trait: float = 0.6, age: int = 2, corridor_available: bool = True):
    return ColonizationRecruitmentContext(
        trait=trait,
        age=age,
        local_density=0.4,
        resource=0.5,
        mate_success=0.6,
        local_room=0.6,
        expected_target_room=0.8,
        corridor_available=corridor_available,
    )


def test_one_step_expected_recruitment_is_exact_product_of_declared_factors():
    factors = one_step_recruitment_factors(_context(), _params(), ColonizationRegime())

    # survival = 0.9 * [1 - 0.6 * (2/8)^2]
    expected_survival = 0.9 * (1.0 - 0.6 * (2.0 / 8.0) ** 2)
    # conception = 0.45 + .20*.6 + .30*.5 - .2*.6
    expected_conception = 0.45 + 0.20 * 0.6 + 0.30 * 0.5 - 0.2 * 0.6
    # linear benefit_shape gives d=z=.6; extinction rate is zero in this case.
    expected_settlement = 0.6 * 1.0 * 0.8 + 0.4 * 0.6

    assert factors.survival_probability == pytest.approx(expected_survival)
    assert factors.conception_probability == pytest.approx(expected_conception)
    assert factors.patch_persistence_probability == 1.0
    assert factors.settlement_factor == pytest.approx(expected_settlement)
    assert factors.local_reproductive_factor == pytest.approx(
        expected_survival * expected_conception
    )
    assert factors.expected_juvenile_recruitment == pytest.approx(
        expected_survival * expected_conception * expected_settlement
    )
    assert factors.factorisation_residual == pytest.approx(0.0)
    assert require_theorem_interior(factors) is factors


def test_patch_extinction_multiplies_end_of_step_settlement_expectation():
    context = _context()
    no_extinction = one_step_recruitment_factors(context, _params(extinction_rate=0.0), ColonizationRegime())
    with_extinction = one_step_recruitment_factors(context, _params(extinction_rate=0.25), ColonizationRegime())

    assert with_extinction.patch_persistence_probability == pytest.approx(0.75)
    assert with_extinction.local_reproductive_factor == pytest.approx(no_extinction.local_reproductive_factor)
    assert with_extinction.settlement_factor == pytest.approx(0.75 * no_extinction.settlement_factor)
    assert with_extinction.expected_juvenile_recruitment == pytest.approx(
        0.75 * no_extinction.expected_juvenile_recruitment
    )


def test_corridor_loss_changes_settlement_factor_without_changing_local_reproductive_factor():
    params = _params()
    before = one_step_recruitment_factors(_context(), params, ColonizationRegime(connectivity_present=1.0))
    after = one_step_recruitment_factors(_context(), params, ColonizationRegime(connectivity_present=0.0))

    assert after.local_reproductive_factor == pytest.approx(before.local_reproductive_factor)
    assert after.settlement_factor < before.settlement_factor
    assert after.expected_juvenile_recruitment < before.expected_juvenile_recruitment

    # N1 instantiated in this life cycle: the E attenuation caused by corridor loss
    # is observationally equivalent in W_recruit to the same attenuation applied
    # to F_local while keeping the pre-loss settlement factor.
    attenuation = after.settlement_factor / before.settlement_factor
    equivalent_local_reproduction = attenuation * before.local_reproductive_factor
    assert after.expected_juvenile_recruitment == pytest.approx(
        equivalent_local_reproduction * before.settlement_factor
    )


def test_unavailable_corridor_matches_zero_connectivity_branch_without_local_fallback():
    params = _params()
    unavailable = one_step_recruitment_factors(
        _context(corridor_available=False), params, ColonizationRegime(connectivity_present=1.0)
    )
    zero_connectivity = one_step_recruitment_factors(
        _context(corridor_available=True), params, ColonizationRegime(connectivity_present=0.0)
    )

    assert unavailable.settlement_factor == pytest.approx(zero_connectivity.settlement_factor)
    # With no corridor, dispersers fail; only non-dispersing offspring use local room.
    assert unavailable.settlement_factor == pytest.approx((1.0 - 0.6) * 0.6)


def test_maximum_age_is_a_valid_biological_boundary_but_not_theorem_interior():
    factors = one_step_recruitment_factors(_context(age=8), _params(), ColonizationRegime())

    assert factors.survival_probability == 0.0
    assert factors.local_reproductive_factor == 0.0
    assert factors.expected_juvenile_recruitment == 0.0
    with pytest.raises(ValueError, match="require positive F and E"):
        require_theorem_interior(factors)


def test_zero_settlement_is_a_valid_boundary_but_not_theorem_interior():
    context = ColonizationRecruitmentContext(
        trait=1.0,
        age=1,
        local_density=0.4,
        resource=0.5,
        mate_success=0.6,
        local_room=0.0,
        expected_target_room=0.0,
        corridor_available=False,
    )
    factors = one_step_recruitment_factors(context, _params(), ColonizationRegime())

    assert factors.settlement_factor == 0.0
    with pytest.raises(ValueError, match="require positive F and E"):
        require_theorem_interior(factors)