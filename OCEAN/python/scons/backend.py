"""Glue between ``jobcomp`` and the SCons executable.

This module deliberately does not import the SCons Python API.  ``jobcomp``
can therefore build (and dry-run) the command even when SCons is not
installed in its own interpreter, as long as the ``scons`` executable is
available.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional

from jobcomp.errors import BuildError
from jobcomp.manifest import MANIFEST_NAME

SCONS_SCRIPT = Path(__file__).resolve().parent / "SConstruct"

#: Targets jobcomp installs next to the executable, like the make backend.
TARGETS = ("croco", "partit", "ncjoin")

#: Directory containing the local ``scons`` package (``OCEAN/python``).
PACKAGE_DIR = Path(__file__).resolve().parent.parent


def environment(ctx) -> Dict[str, str]:
    """Environment pointing SCons at the local package and build manifest."""
    env = dict(os.environ)
    env["CROCO_SCONS_PACKAGE"] = str(PACKAGE_DIR)
    env["CROCO_BUILD_CONFIG"] = str(
        ctx.paths.scratch.resolve() / MANIFEST_NAME
    )
    return env


def scons_argv() -> Optional[List[str]]:
    """How to invoke SCons, or ``None`` when it is not installed."""
    executable = shutil.which("scons")
    if executable:
        return [executable]
    if importlib.util.find_spec("SCons") is not None:
        return [sys.executable, "-m", "SCons"]
    return None


def command(ctx, *, require: bool = True) -> List[str]:
    """The SCons command line for a :class:`BuildContext`-like object."""
    argv = scons_argv()
    if argv is None:
        if require:
            raise BuildError(
                "SCons is not installed. Install it with "
                "`python3 -m pip install -r "
                "OCEAN/python/scons/requirements.txt` or use the default "
                "make backend."
            )
        # A dry run should show the plan even on a machine without SCons.
        argv = ["scons"]
    return [
        *argv,
        "-f",
        str(SCONS_SCRIPT),
        "-j",
        str(ctx.tools.nprocs),
        *TARGETS,
    ]
