
from __future__ import annotations
from SCons.Scanner import Classic, Current
from SCons.Node import Node
from SCons.Environment import Environment
from .deps_fcompiler import FCompiler


class _F90Scanner(Classic):
    """
    See the SCons Fortran scanner documentation.

    Modified F90 scanner from SCons documentation to get all 'use'
    and 'include' statements. Note that the 'include' statements are compiler
    extensions in Fortran and not precompilers.
    """

    def __init__(
        self,
        name,
        suffixes,
        path_function,
        *args,
        **kwargs,
    ) -> None:

        def _scan_function(node, env, path, self=self):
            # Return the nodes filename
            node = node.rfile()

            # If it doesn't exist, return empty list
            if not node.exists():
                return []

            return self.scan(node, env, path)

        kwargs["function"] = _scan_function
        kwargs["path_function"] = path_function
        kwargs["recursive"] = 1
        kwargs["skeys"] = suffixes
        kwargs["name"] = name

        # Bypasses the parent classic initializer
        Current.__init__(self, *args, **kwargs)

    def scan(self, node: Node, env: Environment, path):

        mod_suffix = env.subst("$FORTRANMODSUFFIX")

        fcompiler = FCompiler(node)

        # Load include files
        inc_files = fcompiler.getIncludes(path)

        # Load 'use [module]' modules
        use_modules = fcompiler.getUseModules()

        # Append .mod file suffix
        use_modules_files = [
            f"{module}{mod_suffix}" for module in use_modules
        ]

        return env.File(inc_files + use_modules_files)


def ScannerFortran(path_function):
    scanner = _F90Scanner(
        name="FortranScan",
        suffixes=[".f", ".f90"],
        path_function=path_function,
    )

    return scanner
