"""Constitutive-theory registry for the eco-genetic criticality program.

A constitutive theory is a rule set for building related models. This registry is
not a biological simulator and does not prove H1--H3. It makes explicit which
propositions a model instantiates, the model's domain, and the strongest claim
category that model may support.

This prevents a result proved in a canonical map, a conditional submodel, or a
stochastic simulation from being silently promoted to a universal claim.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

ClaimType = Literal["T", "C", "H", "S"]


@dataclass(frozen=True)
class Proposition:
    identifier: str
    statement: str
    role: str


@dataclass(frozen=True)
class ModelInstantiation:
    identifier: str
    domain: str
    propositions: tuple[str, ...]
    strongest_claim: ClaimType
    supports: tuple[str, ...]
    excludes: tuple[str, ...]


PROPOSITIONS: tuple[Proposition, ...] = (
    Proposition(
        "CT1",
        "Patch configuration constrains local interaction support through a declared state equation.",
        "ecological support",
    ),
    Proposition(
        "CT2",
        "Interaction state changes the potential performance or viability of trait states.",
        "trait-space mechanism",
    ),
    Proposition(
        "CT3",
        "Realised trait occupancy is a finite-population state distinct from potential viability.",
        "potential-realised separation",
    ),
    Proposition(
        "CT4",
        "Census size, effective reproductive size, and genetic diversity are distinct state variables.",
        "demographic-genetic separation",
    ),
    Proposition(
        "CT5",
        "Finite recruitment and allele transmission create stochastic first-passage events.",
        "stochasticity",
    ),
    Proposition(
        "CT6",
        "Connectivity can preserve, erode, or rescue lower bounds; its effect must be derived from the stated migration update.",
        "connectivity",
    ),
    Proposition(
        "CT7",
        "Potential viability, realised occupancy, allele persistence, and diversity warnings require separate observables and event definitions.",
        "measurement discipline",
    ),
)

MODEL_INSTANTIATIONS: tuple[ModelInstantiation, ...] = (
    ModelInstantiation(
        "canonical_logistic_criticality",
        "one-dimensional deterministic interaction map with declared performance surface",
        ("CT1", "CT2", "CT7"),
        "T",
        ("H1", "H3"),
        ("finite occupancy", "genetic lead", "migration rescue"),
    ),
    ModelInstantiation(
        "finite_bin_coupled_feedback",
        "finite-bin multipatch stochastic closure with optional two-kernel recruitment",
        ("CT1", "CT2", "CT3", "CT4", "CT5", "CT6", "CT7"),
        "S",
        ("H1", "H2", "H3"),
        ("universal theorem", "empirical ecosystem projection"),
    ),
    ModelInstantiation(
        "conditional_h_alpha_lead",
        "abstract finite-population process with declared diversity-decay and trait-persistence bounds",
        ("CT3", "CT4", "CT5", "CT7"),
        "C",
        ("H2",),
        ("automatic simulator applicability", "universal lead ordering"),
    ),
    ModelInstantiation(
        "stochastic_refuge_no_migration",
        "single-patch finite-bin refuge over a finite horizon with zero migration",
        ("CT1", "CT2", "CT3", "CT4", "CT5", "CT7"),
        "C",
        ("H2",),
        ("arbitrary-migration metapopulation theorem", "infinite-horizon persistence"),
    ),
    ModelInstantiation(
        "migration_safe_lower_bound",
        "convex-combination allele migration update with declared focal and source lower bounds",
        ("CT4", "CT6", "CT7"),
        "T",
        ("H2", "H3"),
        ("source-mean invariance", "whole-network persistence"),
    ),
)


def proposition_ids() -> tuple[str, ...]:
    """Return the unique constitutive-theory proposition identifiers."""
    return tuple(proposition.identifier for proposition in PROPOSITIONS)


def model_by_identifier(identifier: str) -> ModelInstantiation:
    """Retrieve one registered model instantiation."""
    for model in MODEL_INSTANTIATIONS:
        if model.identifier == identifier:
            return model
    raise KeyError(identifier)


def validate_constitutive_registry() -> tuple[str, ...]:
    """Return registry violations; an empty tuple means the blueprint is coherent."""
    identifiers = proposition_ids()
    if len(set(identifiers)) != len(identifiers):
        return ("proposition identifiers must be unique",)
    valid_claims = {"T", "C", "H", "S"}
    violations: list[str] = []
    for model in MODEL_INSTANTIATIONS:
        unknown = sorted(set(model.propositions).difference(identifiers))
        if unknown:
            violations.append(f"{model.identifier} uses unknown propositions: {', '.join(unknown)}")
        if model.strongest_claim not in valid_claims:
            violations.append(f"{model.identifier} has invalid claim type")
        if not model.supports:
            violations.append(f"{model.identifier} must map to at least one hypothesis")
        if not model.excludes:
            violations.append(f"{model.identifier} must state at least one scope exclusion")
    return tuple(violations)
