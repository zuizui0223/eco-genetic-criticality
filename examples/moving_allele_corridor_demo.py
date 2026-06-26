"""Demonstrate a finite-horizon moving allele-corridor certificate."""
from causal_model.moving_allele_corridor_theory import (
    moving_corridor_certificate,
    selection_shifted_corridor,
)
from causal_model.multipatch_criticality_dynamics import DynamicsParameters


def main() -> None:
    parameters = DynamicsParameters(
        patch_areas=(1.0, 1.0, 1.0),
        high_base=2.0,
        high_interaction_benefit=0.0,
        viability_threshold=1.0,
        selection_strength=0.5,
        effective_fraction=1.0,
        migration_rate=0.8,
    )
    lower, upper = selection_shifted_corridor(
        parameters,
        lower_initial=0.60,
        upper_initial=0.70,
        interaction_lower_bound=0.0,
        interaction_upper_bound=0.0,
        generations=3,
        slack=0.05,
    )
    certificate = moving_corridor_certificate(
        parameters,
        patches=3,
        lower_path=lower,
        upper_path=upper,
        interaction_lower_bound=0.0,
        interaction_upper_bound=0.0,
        population_lower_bound=1000,
    )
    print("Moving allele corridor")
    print(f"  lower path: {tuple(round(value, 4) for value in lower)}")
    print(f"  upper path: {tuple(round(value, 4) for value in upper)}")
    print(f"  horizon retention lower bound: {certificate.horizon_retention_probability_lower_bound:.6f}")
    print(f"  certified: {certificate.certified}")


if __name__ == "__main__":
    main()
