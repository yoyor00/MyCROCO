"""JSON build manifest shared between jobcomp and the SCons backend.

``Makedefs`` is make syntax; SCons used to parse it and the ``Makefile`` again.
This module serializes the already-resolved build context into a small,
versioned JSON document instead.  ``jobcomp`` is the only writer; SCons is a
read-only consumer.
"""

from __future__ import annotations

import json
import os
import shlex
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Mapping

from . import __version__
from .flags import classify_flags

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .context import BuildContext

SCHEMA_VERSION = 1
MANIFEST_NAME = "croco_build.json"

#: Informational list of helper programs the SCons backend knows how to build.
SCONS_TOOLS = (
    "mpc",
    "cross_matrix",
    "cppcheck",
    "srcscheck",
    "partit",
    "ncjoin",
    "ncrename",
)

_ENVIRONMENT_KEYS = (
    "PATH",
    "HOME",
    "LD_LIBRARY_PATH",
    "PKG_CONFIG_PATH",
    "C_INCLUDE_PATH",
)


def _environment() -> Dict[str, str]:
    """The subset of the OS environment the build may need."""
    values = {
        name: os.environ[name]
        for name in _ENVIRONMENT_KEYS
        if name in os.environ
    }
    values.update(
        {
            name: value
            for name, value in sorted(os.environ.items())
            if name.startswith("I_")
        }
    )
    return values


def _flags(value: str) -> List[str]:
    return shlex.split(value) if value else []


class BuildManifest:
    """Resolved build configuration, serialized by jobcomp."""

    def __init__(self, data: Mapping[str, Any]) -> None:
        self.data = dict(data)

    @classmethod
    def from_context(
        cls,
        ctx: "BuildContext",
        *,
        targets: List[str],
        sources: Mapping[str, List[str]],
    ) -> "BuildManifest":
        tools = ctx.tools
        cpp = _flags(tools.cpp)
        fflags = _flags(tools.fflags)
        ldflags = _flags(tools.ldflags)
        defines, includes, library_dirs, libraries = classify_flags(
            tools.cpp, tools.cppflags, tools.fflags, tools.ldflags
        )
        return cls(
            {
                "schema_version": SCHEMA_VERSION,
                "generator": {"name": "jobcomp", "version": __version__},
                "build": {
                    "backend": ctx.backend.value,
                    "source_dir": str(ctx.paths.source.resolve()),
                    "scratch_dir": str(ctx.paths.scratch.resolve()),
                    "run_dir": str(ctx.paths.run_dir.resolve()),
                    "agrif_dir": str(ctx.paths.agrif.resolve()),
                    "jobs": tools.nprocs,
                    "copy_mode": ctx.copy_mode.value,
                    "targets": list(targets),
                    "tools": list(SCONS_TOOLS),
                },
                "toolchain": {
                    "os": tools.os,
                    "fc": tools.fc,
                    "fc_invoked": tools.cft,
                    "ldr": tools.cft,
                    "mpif90": tools.mpif90,
                    "cc": tools.cc,
                    "make": tools.make,
                    "cpp": cpp,
                },
                "flags": {
                    "cppflags": _flags(tools.cppflags),
                    "fflags": fflags,
                    "fflags90": list(fflags),
                    "ldflags": ldflags,
                    "cflags": _flags(tools.cflags),
                    "defines": list(defines),
                    "include_dirs": list(includes),
                    "library_dirs": list(library_dirs),
                    "libraries": list(libraries),
                },
                "features": dict(ctx.features),
                "sources": {
                    "F": list(sources.get("F", [])),
                    "F90": list(sources.get("F90", [])),
                },
                "environment": _environment(),
            }
        )

    def to_dict(self) -> Dict[str, Any]:
        return dict(self.data)

    def write(self, path: Path, *, dry_run: bool = False) -> None:
        """Write the manifest atomically (or just log it in a dry run)."""
        path = Path(path)
        if dry_run:
            return
        path.parent.mkdir(parents=True, exist_ok=True)
        text = json.dumps(self.data, indent=2, sort_keys=True) + "\n"
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(path)
