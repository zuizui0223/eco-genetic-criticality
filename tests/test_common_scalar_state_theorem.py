from __future__ import annotations

import csv
from itertools import product
from pathlib import Path

from causal_model.state_scalarization import (
    common_monotone_scalar_audit,
    scalarization_is_valid,
)

ROOT = Path(__file__).resolve().parents[1]


def brute_scalar_exists(target_values, target_names):
    states = tuple(target_values)
    n = len(states)
    for labels in product(range(n), repeat=n):
        scalar = dict(zip(states, labels))
        if scalarization_is_valid(
            target_values=target_values,
            target_names=target_names,
            scalar_state=scalar,
        ):
            return True
    return False


def test_chain_criterion_matches_exhaustive_scalar_search_on_all_small_binary_tables():
    states = ("a", "b", "c")
    target_names = ("ecology", "genetics")
    # Every 3-state, 2-target table with binary target values.
    for flat in product((0.0, 1.0), repeat=6):
        target_values = {
            state: {
                "ecology": flat[2 * i],
                "genetics": flat[2 * i + 1],
            }
            for i, state in enumerate(states)
        }
        audit = common_monotone_scalar_audit(
            target_values=target_values, target_names=target_names
        )
        oracle = brute_scalar_exists(target_values, target_names)
        assert audit.chain_under_product_order == oracle
        if audit.chain_under_product_order:
            assert audit.scalar_state is not None
            assert scalarization_is_valid(
                target_values=target_values,
                target_names=target_names,
                scalar_state=audit.scalar_state,
            )
        else:
            assert audit.crossing_pair is not None


def test_crossing_pair_is_an_impossibility_certificate():
    targets = {
        "two_patches": {"interaction": 0.8, "trait_mass": 0.2},
        "many_patches": {"interaction": 0.3, "trait_mass": 0.6},
    }
    audit = common_monotone_scalar_audit(
        target_values=targets,
        target_names=("interaction", "trait_mass"),
    )
    assert not audit.chain_under_product_order
    assert audit.scalar_state is None
    assert set(audit.crossing_pair or ()) == set(targets)
    assert not brute_scalar_exists(targets, ("interaction", "trait_mass"))


def test_duplicate_target_vectors_can_share_one_scalar_state():
    targets = {
        "a": {"x": 0.0, "y": 0.0},
        "b": {"x": 0.0, "y": 0.0},
        "c": {"x": 1.0, "y": 2.0},
    }
    audit = common_monotone_scalar_audit(
        target_values=targets, target_names=("x", "y")
    )
    assert audit.chain_under_product_order
    assert audit.scalar_state is not None
    assert audit.scalar_state["a"] == audit.scalar_state["b"]
    assert audit.scalar_state["a"] < audit.scalar_state["c"]


def test_locked_h3_gradient_contains_two_vs_sixteen_crossing_certificate():
    path = ROOT / "artifacts/h3_fragmentation_gradient/h3_fragmentation_gradient_pooled_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {int(row["patch_count"]): row for row in csv.DictReader(handle)}

    def target(row):
        return {
            "interaction": float(row["final_interaction_mean_ratio_to_n1_median"]),
            "effective_size": float(row["final_effective_size_mean_ratio_to_n1_median"]),
            "trait_mass": float(row["realised_high_trait_mass_mean_ratio_to_n1_median"]),
        }

    targets = {2: target(rows[2]), 16: target(rows[16])}
    assert targets[2]["interaction"] > targets[16]["interaction"]
    assert targets[2]["effective_size"] > targets[16]["effective_size"]
    assert targets[2]["trait_mass"] < targets[16]["trait_mass"]

    audit = common_monotone_scalar_audit(
        target_values=targets,
        target_names=("interaction", "effective_size", "trait_mass"),
    )
    assert not audit.chain_under_product_order
    assert set(audit.crossing_pair or ()) == {2, 16}


def test_locked_h3_crossing_already_blocks_scalarization_of_just_two_targets():
    path = ROOT / "artifacts/h3_fragmentation_gradient/h3_fragmentation_gradient_pooled_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = {int(row["patch_count"]): row for row in csv.DictReader(handle)}
    targets = {
        n: {
            "interaction": float(rows[n]["final_interaction_mean_ratio_to_n1_median"]),
            "trait_mass": float(rows[n]["realised_high_trait_mass_mean_ratio_to_n1_median"]),
        }
        for n in (2, 16)
    }
    audit = common_monotone_scalar_audit(
        target_values=targets, target_names=("interaction", "trait_mass")
    )
    assert not audit.chain_under_product_order


def test_restricted_chain_can_lose_scalarizability_when_domain_expands():
    restricted = {
        2: {"interaction": 0.001744, "trait_mass": 0.282918},
        4: {"interaction": 0.001447, "trait_mass": 0.301800},
    }
    # This pair already crosses, so use a truly chain-like restricted toy domain
    # to test the theorem's domain-expansion clause independently of H3 values.
    chain = {
        "low": {"interaction": 0.1, "trait_mass": 0.2},
        "high": {"interaction": 0.2, "trait_mass": 0.3},
    }
    assert common_monotone_scalar_audit(
        target_values=chain, target_names=("interaction", "trait_mass")
    ).chain_under_product_order
    expanded = {
        **chain,
        "crossing": {"interaction": 0.15, "trait_mass": 0.4},
    }
    assert not common_monotone_scalar_audit(
        target_values=expanded, target_names=("interaction", "trait_mass")
    ).chain_under_product_order
