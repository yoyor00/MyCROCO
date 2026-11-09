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
        yield CopySpec.glob(
            ctx.paths.root / "XIOS", ctx.paths.scratch, "*.F"
        )

    def configure(self, ctx: BuildContext) -> None:
        tools = ctx.tools
        root = ctx.xios_dir
        tools.add_ldflags(f"{root}/lib/libxios.a -lstdc++ -lnetcdff -lnetcdf")
        tools.add_cppflags(f"-I{root}/inc")
        tools.add_fflags(f"-I{root}/inc")

    def prepare(self, ctx: BuildContext) -> None:
        self.symlink(
            ctx, Path(f"{ctx.xios_dir}/bin/xios_server.exe"), ctx.run_dir
        )
