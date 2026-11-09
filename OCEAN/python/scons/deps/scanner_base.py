"""Shared SCons scanner setup for the CROCO dependency scanners."""

from __future__ import annotations

from SCons.Environment import Environment
from SCons.Node import Node
from SCons.Scanner import Classic, Current


class CrocoScanner(Classic):
    """Base scanner that skips missing files and delegates to ``scan``."""

    name = "CrocoScanner"
    suffixes = ()

    def __init__(self, path_function) -> None:
        def scan_function(node, env, path, self=self):
            node = node.rfile()
            if not node.exists():
                return []
            return self.scan(node, env, path)

        Current.__init__(
            self,
            function=scan_function,
            path_function=path_function,
            recursive=1,
            skeys=self.suffixes,
            name=self.name,
        )

    def scan(self, node: Node, env: Environment, path):
        raise NotImplementedError
