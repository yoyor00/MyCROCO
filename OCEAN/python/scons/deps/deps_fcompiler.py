"""Extract Fortran module and ``include`` dependencies from source files."""

from __future__ import annotations

import re
from typing import Callable, List

from SCons import Util
from SCons.Node.FS import FS

from .deps_preprocessor import Preprocessor

_NO_DEPS_MODULES = frozenset(
    {
        "openacc",
        "openmp",
        "omp_lib",
        "mpi",
        "mpi_f08",
        "mpi_f77",
        "mpi_f90",
        "perf_regions_fortran",
        "iso_fortran_env",
    }
)

_USE_RE = re.compile(
    r"(?i)(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*"
    r"(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)",
    re.M,
)
_INCLUDE_RE = re.compile(
    r"""(?i)(?:^|['">]\s*;)\s*INCLUDE\s+"""
    r"""(?:\w+_)?[<"'](.+?)(?=["'>])""",
    re.M,
)
_MODULE_RE = re.compile(
    r"""(?i)^\s*MODULE\s+"""
    r"""(?!PROCEDURE|SUBROUTINE|FUNCTION|PURE|ELEMENTAL)(\w+)""",
    re.M,
)


class FCompiler:
    """Find ``use`` modules and Fortran ``include`` files in one source."""

    def __init__(self, node: FS.File):
        self.node = node

    def _contents(self) -> str:
        if not hasattr(self.node.attributes, "contents"):
            self.node.attributes.contents = self.node.get_text_contents()
        return self.node.attributes.contents

    def _cached(
        self, name: str, factory: Callable[[], List[str]]
    ) -> List[str]:
        if hasattr(self.node.attributes, name):
            return getattr(self.node.attributes, name)
        value = factory()
        setattr(self.node.attributes, name, value)
        return value

    def use_modules(self) -> List[str]:
        """Return the modules referenced by ``use`` statements."""
        return self._cached("cache_use_modules", self._find_use_modules)

    def _find_use_modules(self) -> List[str]:
        modules = [_.lower() for _ in _USE_RE.findall(self._contents())]
        modules = Util.unique(modules)
        modules = [m for m in modules if m not in _NO_DEPS_MODULES]
        defined = self.defined_modules()
        return [m for m in modules if m not in defined]

    def includes(self, include_paths: List[FS.Dir]) -> List[str]:
        """Return Fortran ``include`` files that exist in the search paths."""
        return self._cached(
            "cache_inc_dep_files",
            lambda: Preprocessor(self.node).gen_deps_exists(
                include_paths, inc_files=_INCLUDE_RE.findall(self._contents())
            ),
        )

    def defined_modules(self) -> List[str]:
        """Return the modules defined by this source file."""
        return self._cached("cache_def_modules", self._find_defined_modules)

    def _find_defined_modules(self) -> List[str]:
        modules = [_.lower() for _ in _MODULE_RE.findall(self._contents())]
        return [m for m in modules if m not in _NO_DEPS_MODULES]
