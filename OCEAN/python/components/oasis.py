"""OASIS-MCT (PRISM): the psmile/mct libraries and their includes."""

from __future__ import annotations

from jobcomp.context import BuildContext

from .base import Component


class OasisComponent(Component):
    name = "oasis"
    key = "oacplisdefined"
    label = "OASIS"
    announce = True

    def configure(self, ctx: BuildContext) -> None:
        tools = ctx.tools
        root = ctx.prism_dir
        channel = "MPI1"
        libraries = (
            f"{root}/lib/libpsmile.{channel}.a "
            f"{root}/lib/libmct.a "
            f"{root}/lib/libmpeu.a "
            f"{root}/lib/libscrip.a"
        )
        includes = (
            f"-I{root}/build/lib/psmile.{channel} "
            f"-I{root}/build/lib/mct"
        )
        tools.add_ldflags(libraries)
        tools.add_cppflags(includes)
        tools.add_fflags(includes)
