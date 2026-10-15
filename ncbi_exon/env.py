"""Minimal .env loader (stdlib only).

Why not python-dotenv? We're stdlib-only for the validation phase. This loader
covers the common case: KEY=value lines, optional surrounding quotes, '#' comments,
no variable interpolation. Switch to python-dotenv later if we need anything fancier.
"""

from __future__ import annotations

import os
from pathlib import Path


def load_dotenv(path: str | Path | None = None, *, override: bool = False) -> dict[str, str]:
    """Load KEY=value pairs from a .env file into os.environ.

    - If `path` is None, walks up from cwd looking for a `.env` file.
    - Existing env vars are NOT overwritten unless `override=True`.
    - Returns the dict of keys that were applied (useful for logging).
    - Never prints values; caller is responsible for any logging.
    """
    env_path = Path(path) if path else _find_dotenv()
    if env_path is None or not env_path.exists():
        return {}

    applied: dict[str, str] = {}
    for raw in env_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip()
        # Strip a single layer of matching quotes.
        if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
            value = value[1:-1]
        if not key:
            continue
        if override or key not in os.environ:
            os.environ[key] = value
            applied[key] = value
    return applied


def _find_dotenv(start: Path | None = None) -> Path | None:
    here = (start or Path.cwd()).resolve()
    for d in (here, *here.parents):
        candidate = d / ".env"
        if candidate.exists():
            return candidate
    return None
