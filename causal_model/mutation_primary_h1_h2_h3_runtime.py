"""Runtime adapter for the mutation primary H1-H2-H3 chain.

The chain module retains a compact positional reconstruction helper while the
canonical certificate API is keyword-only.  This adapter patches only that
module-local name before dispatching the public runner; it does not alter the
canonical theorem code or any legacy audit.
"""
from __future__ import annotations

import causal_model.mutation_primary_h1_h2_h3_chain as chain
from causal_model.canonical_h1_bifurcation import canonical_h1_certificate as _certificate
from causal_model.multipatch_criticality_dynamics import DynamicsParameters


def _certificate_adapter(
    feedback_strength: float,
    area: float,
    area_reference: float,
    barrier: float,
    trait_parameters: DynamicsParameters,
):
    return _certificate(
        feedback_strength=feedback_strength,
        area=area,
        area_reference=area_reference,
        barrier=barrier,
        trait_parameters=trait_parameters,
    )


chain.canonical_h1_certificate = _certificate_adapter

run_mutation_primary_h1_h2_h3_chain = chain.run_mutation_primary_h1_h2_h3_chain
write_mutation_primary_h1_h2_h3_artifacts = chain.write_mutation_primary_h1_h2_h3_artifacts
DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS = chain.DEFAULT_PRIMARY_CHAIN_MASTER_SEEDS
