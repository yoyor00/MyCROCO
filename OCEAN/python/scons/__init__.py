"""SCons backend for the CROCO build.

The heavy SCons-dependent classes are imported lazily so that
``scons.backend`` (used by ``jobcomp --backend scons``) can be inspected and
dry-run without importing the SCons API.  The ``SConstruct`` in this directory
imports the public classes through the module ``__getattr__`` below.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .config import CrocoBuildConfig
    from .croco_bin import CrocoBin
    from .croco_build import CrocoBuild
    from .croco_env import CrocoEnv
    from .croco_tools import CrocoTools

__all__ = [
    "CrocoBin",
    "CrocoBuild",
    "CrocoBuildConfig",
    "CrocoEnv",
    "CrocoTools",
]

_LAZY = {
    "CrocoBin": (".croco_bin", "CrocoBin"),
    "CrocoBuild": (".croco_build", "CrocoBuild"),
    "CrocoBuildConfig": (".config", "CrocoBuildConfig"),
    "CrocoEnv": (".croco_env", "CrocoEnv"),
    "CrocoTools": (".croco_tools", "CrocoTools"),
}


def __getattr__(name: str):
    target = _LAZY.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from importlib import import_module

    module = import_module(target[0], __name__)
    value = getattr(module, target[1])
    globals()[name] = value
    return value


__version__ = "1.0"
