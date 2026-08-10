from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tomllib
from pathlib import Path


def _git(repo: Path, *args: str) -> str:
    completed = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _fail(message: str) -> None:
    raise AssertionError(message)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Verify the frozen eco-genetic-criticality release contract."
    )
    parser.add_argument("--repo-root", default=".")
    parser.add_argument(
        "--manifest", default="reproducibility/release_manifest.json"
    )
    args = parser.parse_args()

    root = Path(args.repo_root).resolve()
    manifest_path = root / args.manifest
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    scientific_commit = manifest["scientific_commit"]

    try:
        _git(root, "cat-file", "-e", f"{scientific_commit}^{{commit}}")
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "canonical commit is unavailable; fetch full history before verification"
        ) from exc

    mismatches: list[str] = []
    for relative_path, expected_blob in manifest["canonical_files"].items():
        actual_blob = _git(root, "rev-parse", f"{scientific_commit}:{relative_path}")
        if actual_blob != expected_blob:
            mismatches.append(
                f"{relative_path}: expected {expected_blob}, found {actual_blob}"
            )
    if mismatches:
        _fail("canonical file mismatch:\n" + "\n".join(mismatches))

    pyproject = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    project = pyproject["project"]
    package = manifest["package"]
    if project["name"] != package["name"]:
        _fail(f"package name drift: {project['name']} != {package['name']}")
    if project["version"] != package["version"]:
        _fail(f"package version drift: {project['version']} != {package['version']}")
    if project["requires-python"] != package["python"]:
        _fail(
            "Python requirement drift: "
            f"{project['requires-python']} != {package['python']}"
        )

    required_current_paths = (
        "README.md",
        "REPRODUCIBILITY.md",
        "docs/final_evidence_ledger.md",
        "manuscript/claim_evidence_map.md",
        "reproducibility/release_manifest.json",
    )
    missing = [path for path in required_current_paths if not (root / path).exists()]
    if missing:
        _fail("missing release paths: " + ", ".join(missing))

    ledger = (root / "docs/final_evidence_ledger.md").read_text(encoding="utf-8")
    for identifier in ("H1", "H3", "H2-R", "H2-A"):
        if identifier not in ledger:
            _fail(f"final evidence ledger no longer names {identifier}")

    print(
        json.dumps(
            {
                "status": "verified",
                "repository": manifest["repository"],
                "scientific_commit": scientific_commit,
                "canonical_file_count": len(manifest["canonical_files"]),
                "package": package,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AssertionError, FileNotFoundError, KeyError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"release-contract verification failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
