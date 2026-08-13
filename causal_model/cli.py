from __future__ import annotations

import importlib
import sys
from collections.abc import Sequence

COMMANDS: dict[str, str] = {
    "branch-locked-chain": "causal_model.branch_locked_h1_h2_h3_chain_cli:main",
    "finite-validation": "causal_model.finite_validation_campaign_cli:main",
    "h1-boundary": "causal_model.finite_h1_boundary_resolution_audit_cli:main",
    "h1-continuation": "causal_model.finite_h1_continuation_state_audit_cli:main",
    "h1-fragment-projection": "causal_model.finite_h1_fragment_projection_audit_cli:main",
    "h1-hysteresis-duration": "causal_model.finite_h1_hysteresis_duration_audit_cli:main",
    "h1-mutation-window": "causal_model.finite_h1_mutation_window_audit_cli:main",
    "h1-polymorphism": "causal_model.finite_h1_polymorphism_eligibility_audit_cli:main",
    "h1-seed-ensemble": "causal_model.finite_h1_boundary_seed_ensemble_cli:main",
    "h1-sweep-boundary": "causal_model.finite_h1_sweep_boundary_audit_cli:main",
    "h1-targeted-validation": "causal_model.h1_targeted_validation_campaign_cli:main",
    "h2a-secondary-audit": "causal_model.h2_absolute_secondary_audit_cli:main",
    "h2r-calibration": "causal_model.h2r_trait_loss_calibration_cli:main",
    "h2r-independent-validation": "causal_model.h2r_independent_relative_validation_cli:main",
    "h2r-ramp-hold-calibration": "causal_model.h2r_ramp_hold_trait_loss_calibration_cli:main",
    "mutation-primary-chain": "causal_model.mutation_primary_h1_h2_h3_chain_cli:main",
    "paired-baseline": "causal_model.paired_baseline_cli:main",
    "theorem-boundary": "causal_model.theorem_boundary_cli:main",
}


def _print_help() -> None:
    print("usage: egc <command> [args ...]")
    print()
    print("Unified entry point for the frozen eco-genetic-criticality command surface.")
    print("Scientific claims remain locked to the canonical scientific commit.")
    print()
    print("commands:")
    for command in sorted(COMMANDS):
        print(f"  {command}")
    print()
    print("Run 'egc <command> --help' for command-specific arguments.")


def _load(target: str):
    module_name, function_name = target.split(":", 1)
    module = importlib.import_module(module_name)
    return getattr(module, function_name)


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args or args[0] in {"-h", "--help", "help"}:
        _print_help()
        return 0
    if args[0] in {"--list", "list"}:
        for command in sorted(COMMANDS):
            print(command)
        return 0

    command = args.pop(0)
    target = COMMANDS.get(command)
    if target is None:
        print(f"unknown command: {command}", file=sys.stderr)
        print("run 'egc --list' to see available commands", file=sys.stderr)
        return 2

    entry_point = _load(target)
    previous_argv = sys.argv
    sys.argv = [f"egc {command}", *args]
    try:
        result = entry_point()
    finally:
        sys.argv = previous_argv
    return 0 if result is None else int(result)


if __name__ == "__main__":
    raise SystemExit(main())
