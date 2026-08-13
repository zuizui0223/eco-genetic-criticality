from __future__ import annotations

import sys
from types import ModuleType

from causal_model import cli


def test_command_registry_is_explicit_and_unique() -> None:
    assert "theorem-boundary" in cli.COMMANDS
    assert "h2r-independent-validation" in cli.COMMANDS
    assert len(cli.COMMANDS) == len(set(cli.COMMANDS))
    assert all(
        target.startswith("causal_model.") and target.endswith(":main")
        for target in cli.COMMANDS.values()
    )


def test_list_does_not_import_campaign_modules(capsys) -> None:
    assert cli.main(["--list"]) == 0
    output = capsys.readouterr().out
    assert "theorem-boundary" in output
    assert "mutation-primary-chain" in output


def test_unknown_command_returns_usage_error(capsys) -> None:
    assert cli.main(["not-a-command"]) == 2
    assert "unknown command" in capsys.readouterr().err


def test_dispatch_preserves_subcommand_arguments(monkeypatch) -> None:
    module = ModuleType("test_fake_egc_cli")
    calls: list[list[str]] = []

    def fake_main() -> int:
        calls.append(list(sys.argv))
        return 7

    module.main = fake_main  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, module.__name__, module)
    monkeypatch.setitem(cli.COMMANDS, "fake", f"{module.__name__}:main")

    assert cli.main(["fake", "--flag", "value"]) == 7
    assert calls == [["egc fake", "--flag", "value"]]
