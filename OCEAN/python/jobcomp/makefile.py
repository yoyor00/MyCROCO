"""Read simple variable assignments from a CROCO ``Makefile``.

The SCons backend used to parse the Makefile itself.  Jobcomp now reads the
source lists once and puts them in the build manifest.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List

_MISSING = object()


def read_variable(
    path: Path, name: str, default: object = _MISSING
) -> str:
    """Return the value of ``name = ...`` from a makefile.

    Backslash continuations are joined and comments are removed, matching the
    subset of GNU make syntax used by CROCO's ``Makefile``.
    """
    pattern = re.compile(rf"^\s*{re.escape(name)}\s*=")
    lines: List[str] = []
    started = False
    for raw_line in Path(path).read_text(encoding="utf-8").splitlines():
        line = re.sub(r"#.*", "", raw_line).rstrip()
        if not started:
            if not pattern.match(line):
                continue
            started = True
            line = line.split("=", 1)[1].strip()
        lines.append(line)
        if not line.endswith("\\"):
            break
    if not lines:
        if default is _MISSING:
            raise ValueError(f"No variable '{name}' found in {path}")
        return str(default)
    return " ".join(lines).replace("\\", " ").strip()


def read_source_lists(path: Path) -> Dict[str, List[str]]:
    """Return the fixed- and free-form Fortran source lists of a Makefile."""
    return {
        "F": read_variable(path, "SRCS").split(),
        "F90": read_variable(path, "SRCS90", "").split(),
    }
