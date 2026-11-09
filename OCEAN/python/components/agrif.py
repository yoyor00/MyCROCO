"""AGRIF: the library build and the ``CROCOFILES`` header tree."""

from __future__ import annotations

import shutil
from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopyMode, CopySpec
from jobcomp.errors import LOG

from .base import Component
from .crocofiles import Crocofiles


class AgrifComponent(Component):
    """The AGRIF library and the ``CROCOFILES`` header tree."""

    name = "agrif"
    key = "agrifisdefined"
    label = "AGRIF"
    announce = True

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        yield CopySpec.tree(ctx.paths.agrif, ctx.paths.scratch)

    def configure(self, ctx: BuildContext) -> None:
        tools = ctx.tools
        tools.add_fflags("-IAGRIF")
        tools.ldflags = f"-LAGRIF -lagrif {tools.ldflags}".strip()

    def prepare(self, ctx: BuildContext) -> None:
        tools = ctx.tools
        # The AGRIF variant of Makedefs carries the extra rules.  It must win
        # whatever the copy mode says: a stale template means a broken build.
        ctx.copier.file(
            ctx.paths.scratch / "Makedefs.generic.AGRIF",
            ctx.paths.scratch / "Makedefs.generic",
            mode=CopyMode.ALWAYS,
        )
        tools.nprocs = 1  # compile the library on one process only
        self.resolve_c_compiler(ctx)
        argv = [
            tools.make,
            "-C",
            "AGRIF",
            f"FC={tools.cft}",
            f"CPP={tools.cpp}",
        ]
        cppflags = self.cppflags(ctx)
        if cppflags:
            argv.append(f"CPPFLAGS={cppflags}")
        argv += [f"FFLAGS={tools.fflags}", f"CFLAGS={tools.cflags}"]
        ctx.shell.run(argv, cwd=ctx.paths.scratch)
        if tools.os == "Darwin":
            ctx.shell.run(
                ["ranlib", "AGRIF/libagrif.a"], cwd=ctx.paths.scratch
            )
        Crocofiles(ctx).generate()

    def cppflags(self, ctx: BuildContext) -> str:
        """Extra preprocessor flags of the AGRIF library (MPI variant)."""
        if not ctx.feature("mpi"):
            return ""
        return f"-DAGRIF_MPI {ctx.tools.mpiinc}".strip()

    # -- the C compiler the AGRIF library is built with --------------------
    def resolve_c_compiler(self, ctx: BuildContext) -> None:
        """Fill ``ctx.tools.cc``/``cflags`` the way the bash script does."""
        tools = ctx.tools
        scratch = ctx.paths.scratch
        cwd = scratch if scratch.is_dir() else None
        cflags = tools.cflags or tools.make_variable(
            ctx.shell, "CFLAGS", cwd=cwd
        )
        cc = tools.cc or tools.make_variable(ctx.shell, "CC", cwd=cwd)
        if cc and shutil.which(cc) and tools.is_gnu_c(ctx.shell, cc):
            LOG.info("Using the GNU C compiler. Adding -fcommon to CFLAGS")
            cflags = f"{cflags} -fcommon -std=gnu11".strip()
        tools.cc = cc
        tools.cflags = cflags
