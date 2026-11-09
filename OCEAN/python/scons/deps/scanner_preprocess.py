"""SCons scanner for preprocessor ``#include`` dependencies."""

from __future__ import annotations

from .deps_preprocessor import Preprocessor
from .scanner_base import CrocoScanner


class _PreprocessScanner(CrocoScanner):
    name = "PreprocessScanner"
    suffixes = (".F", ".F90", ".h", ".inc")

    def scan(self, node, env, path):
        return env.File(Preprocessor(node).gen_deps_exists(path))


def ScannerPreprocess(path_function) -> CrocoScanner:
    return _PreprocessScanner(path_function)
