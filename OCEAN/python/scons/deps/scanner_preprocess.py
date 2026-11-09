
from __future__ import annotations
from SCons.Scanner import Classic, Current
from .deps_preprocessor import Preprocessor


class _PreprocessScanner(Classic):
    """
    https://scons.org/doc/HTML/scons-api/_modules/SCons/Scanner/Fortran/#F90Scanner

    Modified F90 scanner from SCons documentation to get all precompiler
    'include' statements
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

        # bypasses the parent Classic initializer
        Current.__init__(self, *args, **kwargs)

    def scan(self, node, env, path):
        """
        Simple preprocessor to determine requirements for Croco
        """

        incs = Preprocessor(node)
        preprocess_inc_files = incs.gen_deps_exists(path)

        return env.File(preprocess_inc_files)


def ScannerPreprocess(path_function):
    scanner = _PreprocessScanner(
        suffixes=[".F", ".F90", ".h", ".inc"],
        name="PreprocessScanner",
        path_function=path_function,
    )

    return scanner
