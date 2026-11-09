"""The order in which the components are asked to do their work.

Keeping the lists explicit makes the ordering requirements (AGRIF library
before OpenMP flags, PISCES namelists before AGRIF, ...) visible instead of
implicit in the script line order.  Adding a feature is a class plus one entry
in the list that describes *when* it runs.

The component instances are shared by all lists: a component is stateless, so
one instance per feature is enough.
"""

from __future__ import annotations

from jobcomp.context import BuildContext

from .agrif import AgrifComponent
from .mpi import MpiComponent
from .mustang import MustangComponent
from .oasis import OasisComponent
from .obstruction import ObstructionComponent
from .openacc import OpenAccComponent
from .openmp import OpenMpComponent
from .pisces import PiscesComponent
from .sources import ModelSources
from .stogen import StogenComponent
from .xios import XiosComponent

SOURCES = ModelSources()
AGRIF = AgrifComponent()
XIOS = XiosComponent()
PISCES = PiscesComponent()
MUSTANG = MustangComponent()
OBSTRUCTION = ObstructionComponent()
STOGEN = StogenComponent()
MPI = MpiComponent()
OASIS = OasisComponent()
OPENMP = OpenMpComponent()
OPENACC = OpenAccComponent()

#: Everything copied into the scratch directory (in this order).
COPY_COMPONENTS = (
    SOURCES,
    AGRIF,
    XIOS,
    PISCES,
    MUSTANG,
    OBSTRUCTION,
    STOGEN,
)

#: Components that drop a namelist in the run directory (never overwritten).
RUN_DIR_COMPONENTS = (STOGEN, PISCES)

#: Flags that must be known before the AGRIF library is built.
CONFIGURE_BEFORE = (AGRIF, MPI, XIOS, OASIS)

#: Pre-build work that touches the filesystem, in the order the bash script
#: did it: PISCES namelists, the AGRIF library and headers, the XIOS server.
PREPARE_COMPONENTS = (PISCES, AGRIF, XIOS)

#: Flags that are only needed for the main compile.
CONFIGURE_AFTER = (OPENMP, OPENACC)


def detect_features(ctx: BuildContext) -> None:
    """Detect every keyed feature once, before any phase uses it."""
    for component in CONFIGURE_BEFORE + CONFIGURE_AFTER:
        component.detect(ctx)


def copy_into_scratch(ctx: BuildContext) -> None:
    for component in COPY_COMPONENTS:
        ctx.copy_many(component.copy_specs(ctx))


def install_run_dir_files(ctx: BuildContext) -> None:
    for component in RUN_DIR_COMPONENTS:
        component.run_dir_files(ctx)


def configure_before(ctx: BuildContext) -> None:
    for component in CONFIGURE_BEFORE:
        if component.enabled(ctx):
            component.configure(ctx)


def prepare(ctx: BuildContext) -> None:
    for component in PREPARE_COMPONENTS:
        if component.enabled(ctx):
            component.prepare(ctx)


def configure_after(ctx: BuildContext) -> None:
    for component in CONFIGURE_AFTER:
        if component.enabled(ctx):
            component.configure(ctx)
