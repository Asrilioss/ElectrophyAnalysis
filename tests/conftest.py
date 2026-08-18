"""Shared pytest fixtures and import helpers for the ElectrophyAnalysis suite.

The pipelines are top-level scripts rather than an installed package, so tests
load them by path. Modules whose heavy optional dependencies (pyabf, efel,
ipfx) are missing are skipped rather than failing, which lets the pure-numeric
tests run in a minimal environment and in CI.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def load_script(name):
    """Import a top-level script by filename, or skip if its deps are absent."""
    path = REPO_ROOT / f"{name}.py"
    if not path.is_file():
        pytest.skip(f"{name}.py not found at repository root")
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        sys.modules.pop(name, None)
        pytest.skip(f"{name}.py requires an unavailable dependency: {exc}")
    return module


@pytest.fixture(scope="session")
def slow_depol():
    return load_script("slow_depol")


@pytest.fixture(scope="session")
def single_ap():
    return load_script("single_ap")
