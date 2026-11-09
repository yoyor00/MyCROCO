"""MPI: the Fortran compiler becomes ``MPIF90``."""

from __future__ import annotations

from jobcomp.context import BuildContext

from .base import Component


class MpiComponent(Component):
    name = "mpi"
    key = "mpiisdefined"
    label = "MPI"
    announce = True

    def configure(self, ctx: BuildContext) -> None:
        tools = ctx.tools
        tools.add_ldflags(tools.mpilib)
        tools.add_cppflags(tools.mpiinc)
        tools.add_fflags(tools.mpiinc)
        tools.cft = tools.mpif90
