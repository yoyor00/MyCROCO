"""SCons scanner for Fortran ``use`` and ``include`` dependencies."""

from __future__ import annotations

from SCons.Environment import Environment
from SCons.Node import Node

from .deps_fcompiler import FCompiler
from .scanner_base import CrocoScanner


class _FortranScanner(CrocoScanner):
    name = "FortranScan"
    suffixes = (".f", ".f90")

    def scan(self, node: Node, env: Environment, path):
        compiler = FCompiler(node)
        includes = compiler.includes(path)
        module_suffix = env.subst("$FORTRANMODSUFFIX")
        modules = []
        for module in compiler.use_modules():
            name = f"{module}{module_suffix}"
            for directory in path:
                candidate = directory.File(name)
                if candidate.exists():
                    name = str(candidate)
                    break
            modules.append(name)
        return env.File(includes + modules)


def ScannerFortran(path_function) -> CrocoScanner:
    return _FortranScanner(path_function)
