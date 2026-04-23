#!/usr/bin/env python3
"""Run AlphaFold 3 using the container's `alphafold3` package but with a
*host-supplied* `alphafold3.model.confidence_types` module.

This lets you keep the AF3 container/install untouched while still emitting
high-precision (7-decimal) ipTM/PAE values in `*_summary_confidences.json`
for downstream ensemble ranking.

Must be run *inside* the AF3 container, with `run_alphafold` importable
on the normal Python path.

Configuration (env vars)
------------------------
    AF3_ROOT
        Filesystem path to your AF3 source clone, added to `sys.path` so
        `import run_alphafold` resolves. Required.

    COPAL_AF3_HOST_CONFIDENCE_TYPES
        Path to the replacement `confidence_types.py` (this directory's
        copy, made readable from inside the container).
        Defaults to the file next to this script.

The wrapper passes through any CLI args to `run_alphafold.main` via
`absl.app.run`.

This file is original COPAL code (MIT-licensed). The replacement
`confidence_types.py` it loads is CC-BY-NC-SA 4.0 (see LICENSE in this
directory).
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

DEFAULT_HOST_CONFIDENCE_TYPES = (
    Path(__file__).resolve().parent / "confidence_types.py"
)


def _resolve_paths() -> tuple[Path, Path]:
    af3_root = os.environ.get("AF3_ROOT")
    if not af3_root:
        raise SystemExit(
            "AF3_ROOT env var is not set. Point it at your AlphaFold 3 source "
            "clone (e.g. /opt/alphafold3 or wherever your install lives)."
        )
    host_conf = Path(
        os.environ.get(
            "COPAL_AF3_HOST_CONFIDENCE_TYPES",
            str(DEFAULT_HOST_CONFIDENCE_TYPES),
        )
    )
    if not host_conf.exists():
        raise SystemExit(
            f"Host confidence_types not found at {host_conf}. "
            "Make this file's directory visible to the container or set "
            "COPAL_AF3_HOST_CONFIDENCE_TYPES."
        )
    return Path(af3_root), host_conf


def main() -> None:
    af3_root, host_conf = _resolve_paths()
    if str(af3_root) not in sys.path:
        sys.path.insert(0, str(af3_root))

    spec = importlib.util.spec_from_file_location(
        "alphafold3.model.confidence_types", str(host_conf),
    )
    if spec is None or spec.loader is None:
        raise SystemExit(f"Cannot load {host_conf} as a module.")
    mod = importlib.util.module_from_spec(spec)
    sys.modules["alphafold3.model.confidence_types"] = mod
    spec.loader.exec_module(mod)

    import run_alphafold  # noqa: WPS433  — must come *after* injection
    from absl import app
    app.run(run_alphafold.main)


if __name__ == "__main__":
    main()
