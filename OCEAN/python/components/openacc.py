"""OpenACC: the accelerator flags of the PGI/NVIDIA compilers."""

from __future__ import annotations

from jobcomp.context import BuildContext
from jobcomp.toolchain import CompilerSpec

from .base import Component


class OpenAccComponent(Component):
    name = "openacc"
    key = "openaccisdefined"
    label = "OPENACC"

    def configure(self, ctx: BuildContext) -> None:
        spec = CompilerSpec.find(ctx.tools.fc)
        if spec is not None and spec.name in ("pgfortran", "nvfortran"):
            ctx.tools.add_fflags("-acc -Minfo=accel")
