#!/usr/bin/env python3
"""
config_loader.py — Centralized config path resolution for BlackCotton.
"""

import os
from pathlib import Path
from typing import Any

import yaml

BASE_DIR = Path(__file__).resolve().parent.parent
CONFIG_DIR = BASE_DIR / "config"

PRODUCTION_CONFIG_NAME = "parameters_production.yaml"
DEVELOPMENT_CONFIG_NAME = "parameters.yaml"
ENV_CONFIG_PATH = "BLACKCOTTON_CONFIG_PATH"


def resolve_config_path(preferred_path: str | Path | None = None) -> Path:
    """Resolve active config path with explicit > env > production > development priority."""
    if preferred_path is not None:
        return Path(preferred_path).expanduser().resolve()

    env_path = os.environ.get(ENV_CONFIG_PATH, "").strip()
    if env_path:
        return Path(env_path).expanduser().resolve()

    production_path = CONFIG_DIR / PRODUCTION_CONFIG_NAME
    if production_path.exists():
        return production_path

    return CONFIG_DIR / DEVELOPMENT_CONFIG_NAME


def load_config(preferred_path: str | Path | None = None) -> dict[str, Any]:
    """Load YAML config using the standard resolution order."""
    config_path = resolve_config_path(preferred_path)
    with open(config_path, "r") as f:
        loaded = yaml.safe_load(f)
    if not isinstance(loaded, dict):
        raise ValueError(f"Config must deserialize to a mapping: {config_path}")
    return loaded

