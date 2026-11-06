"""XIOS: I/O server sources, libraries and the ``xios_server.exe`` link."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec

from .base import Component


class XiosComponent(Component):
    name = "xios"
    key = "xiosisdefined"
    label = "XIOS"
    announce = True

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        # XIOS coupling sources (e.g. oxide.F) live in the XIOS directory, not
        # in OCEAN; the bash script copies them unconditionally when present.
        yield CopySpec(src=ctx.paths.root / "XIOS", dst=ctx.paths.scratch,
                       glob="*.F", optional=True)

    def configure(self, ctx: BuildContext) -> None:
        if not self.enabled(ctx):
            return
        tools = ctx.tools
        root = ctx.xios_dir
        tools.add_ldflags(f"{root}/lib/libxios.a -lstdc++ -lnetcdff -lnetcdf")
        tools.add_cppflags(f"-I{root}/inc")
        tools.add_fflags(f"-I{root}/inc")

    def prepare(self, ctx: BuildContext) -> None:
        if not ctx.feature(self.name):
            return
        self.symlink(
            ctx, Path(f"{ctx.xios_dir}/bin/xios_server.exe"), ctx.run_dir
        )
