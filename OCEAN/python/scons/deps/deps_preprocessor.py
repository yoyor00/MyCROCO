"""Find preprocessor ``#include`` dependencies for the SCons scanner."""

from __future__ import annotations

import re
import sys
from typing import Iterable, List, Optional

from SCons.Node.FS import FS

_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*["\']([^"\']+)["\']', re.M)

#: Compiler-provided headers that are not part of the CROCO source tree.
_COMPILER_INCLUDES = frozenset(
    {
        "mpif.h",
        "mpi_f08.h",
        "omp_lib.h",
        "openacc.h",
    }
)


class Preprocessor:
    """Extract existing ``#include`` dependencies from a source node."""

    def __init__(self, node: FS.File):
        self.node = node

    def _contents(self) -> str:
        if not hasattr(self.node.attributes, "contents"):
            self.node.attributes.contents = self.node.get_text_contents()
        return self.node.attributes.contents

    def gen_deps_exists(
        self,
        paths: List[FS.Dir],
        inc_files: Optional[Iterable[str]] = None,
    ) -> List[str]:
        """Return the include files that exist in ``paths``."""
        if hasattr(self.node.attributes, "cache_preprocess_inc_files"):
            return self.node.attributes.cache_preprocess_inc_files

        includes = list(inc_files) if inc_files is not None else None
        if includes is None:
            includes = _INCLUDE_RE.findall(self._contents())

        existing = []
        for include in includes:
            if include in _COMPILER_INCLUDES:
                continue
            for directory in paths:
                candidate = directory.File(include)
                if not candidate.exists():
                    continue
                path = str(candidate)
                if path.startswith("./"):
                    path = path[2:]
                existing.append(path)
                break
            else:
                print(
                    "WARNING: Preprocessor scanner - File "
                    f"'{include}' not found in paths "
                    f"'{[str(_) for _ in paths]}'",
                    file=sys.stderr,
                )

        self.node.attributes.cache_preprocess_inc_files = existing
        return existing
