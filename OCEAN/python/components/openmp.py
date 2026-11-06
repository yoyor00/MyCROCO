"""OpenMP: the per-compiler ``-openmp`` spelling."""

from __future__ import annotations

import re

from jobcomp.context import BuildContext
from jobcomp.toolchain import CompilerSpec

from .base import Component


class OpenMpComponent(Component):
    name = "openmp"
    key = "openmp"
    label = "OMP"

    def configure(self, ctx: BuildContext) -> None:
        if not self.enabled(ctx):
            return
        if ctx.tools.os not in ("Linux", "Darwin"):
            return
        ctx.tools.add_fflags(self.flag_for(ctx, ctx.tools.fc))

    @staticmethod
    def flag_for(ctx: BuildContext, fc: str) -> str:
        """The ``-openmp`` spelling of ``fc`` (Intel classic: version)."""
        spec = CompilerSpec.find(fc)
        if spec is None:
            return "-openmp"
        if spec.name in ("ifort", "ifc"):
            proc = ctx.shell.probe(["ifort", "--version"])
            match = re.search(r"(\d+)", proc.stdout + proc.stderr)
            if match and int(match.group(1)) > 18:
                return "-qopenmp"
            return "-openmp"
        return spec.openmp or "-openmp"
