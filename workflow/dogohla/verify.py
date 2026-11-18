"""Integrity verifier for vendored frozen code in DogoHLA workflow.

Checks every file listed in MANIFEST.sha256 against its recorded SHA-256 hash.
Must be called at every start before importing frozen modules.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

FROZEN_DIR = Path(__file__).resolve().parent / "frozen"
MANIFEST_PATH = FROZEN_DIR / "MANIFEST.sha256"


def verify() -> dict[str, int]:
    """Check every vendored file in MANIFEST.sha256 against its sha256 checksum."""
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(f"Missing frozen manifest: {MANIFEST_PATH}")

    checked = 0
    bad = {}
    for line in MANIFEST_PATH.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        expected_sha, rel_path = parts[0], parts[1].strip()
        file_path = FROZEN_DIR / rel_path
        if not file_path.exists():
            bad[rel_path] = "MISSING"
            continue
        actual_sha = hashlib.sha256(file_path.read_bytes()).hexdigest()
        if actual_sha != expected_sha:
            bad[rel_path] = f"MISMATCH (expected {expected_sha}, got {actual_sha})"
        checked += 1

    if bad:
        raise AssertionError(f"Vendored frozen files failed integrity check ({len(bad)} failures):\n" +
                             "\n".join(f"  {k}: {v}" for k, v in sorted(bad.items())))

    return {"status": "ok", "files_checked": checked}


def use(stage: str | None = None) -> Path:
    """Verify manifest and add frozen directories to sys.path in priority order."""
    verify()
    # Add base package dir
    pkg_dir = Path(__file__).resolve().parent
    if str(pkg_dir) not in sys.path:
        sys.path.insert(0, str(pkg_dir))
    
    # Priority order for frozen modules:
    # 1. Specified stage (if any)
    # 2. mosaic-allgenes-20261114
    # 3. reads-repull-20261114
    # 4. dose-response-20261114
    # 5. mosaic-typer-20261113
    # 6. baselines
    # 7. confirmatory-test-20261113
    stages = [
        "mosaic-allgenes-20261114",
        "reads-repull-20261114",
        "dose-response-20261114",
        "mosaic-typer-20261113",
        "baselines",
        "confirmatory-test-20261113"
    ]
    if stage and stage in stages:
        stages.remove(stage)
        stages.insert(0, stage)

    for st in reversed(stages):
        p = str(FROZEN_DIR / st)
        if p not in sys.path:
            sys.path.insert(0, p)

    return FROZEN_DIR


if __name__ == "__main__":
    res = verify()
    print(f"Verified {res['files_checked']} vendored frozen files successfully.")
