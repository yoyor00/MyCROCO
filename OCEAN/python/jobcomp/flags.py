"""Helpers for turning resolved flag strings into structured lists."""

from __future__ import annotations

import shlex
from typing import List, Tuple


def classify_flags(*flag_groups: str) -> Tuple[Tuple[str, ...], ...]:
    """Split flags into ``(defines, includes, -L, libraries)``, in order."""
    defines: List[str] = []
    includes: List[str] = []
    library_dirs: List[str] = []
    libraries: List[str] = []
    for group in flag_groups:
        tokens = shlex.split(group)
        index = 0
        while index < len(tokens):
            token = tokens[index]
            following = tokens[index + 1] if index + 1 < len(tokens) else None
            index += 1
            if token in ("-I", "-L", "-l", "-D") and following is not None:
                token, index = token + following, index + 1
            if token.startswith("-I"):
                includes.append(token[2:])
            elif token.startswith("-L"):
                library_dirs.append(token[2:])
            elif token.startswith("-l"):
                libraries.append(token[2:])
            elif token.startswith("-D"):
                defines.append(token[2:])
            elif token.endswith(".a"):
                libraries.append(token)
    return (
        tuple(dict.fromkeys(defines)),
        tuple(dict.fromkeys(includes)),
        tuple(dict.fromkeys(library_dirs)),
        tuple(dict.fromkeys(libraries)),
    )
