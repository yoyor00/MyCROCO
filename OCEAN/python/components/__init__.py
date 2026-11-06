"""The independent features of the CROCO build.

One class per feature (AGRIF, XIOS, PISCES, MUSTANG, OBSTRUCTION, STOGEN,
OASIS, MPI, OpenMP, OpenACC and the sources themselves), each in its own
module and each declaring only what it needs:

``copy_specs``      what goes into ``Compile/``
``run_dir_files``   namelists / static files installed next to the executable
``configure``       the flags it adds to the toolchain
``prepare``         pre-build work that is not a simple copy

:mod:`components.pipelines` holds the order in which they are asked, which is
the order the bash script used.
"""

from __future__ import annotations

from .agrif import AgrifComponent
from .base import CPP_OUTPUT_DROP, Component
from .mpi import MpiComponent
from .mustang import MustangComponent
from .oasis import OasisComponent
from .obstruction import ObstructionComponent
from .openacc import OpenAccComponent
from .openmp import OpenMpComponent
from .pipelines import (
    CONFIGURE_AFTER,
    CONFIGURE_BEFORE,
    COPY_COMPONENTS,
    PREPARE_COMPONENTS,
    RUN_DIR_COMPONENTS,
)
from .pisces import PiscesComponent
from .sources import ModelSources
from .stogen import StogenComponent
from .xios import XiosComponent

__all__ = [
    "CPP_OUTPUT_DROP",
    "Component",
    "AgrifComponent",
    "ModelSources",
    "MpiComponent",
    "MustangComponent",
    "OasisComponent",
    "ObstructionComponent",
    "OpenAccComponent",
    "OpenMpComponent",
    "PiscesComponent",
    "StogenComponent",
    "XiosComponent",
    "CONFIGURE_AFTER",
    "CONFIGURE_BEFORE",
    "COPY_COMPONENTS",
    "PREPARE_COMPONENTS",
    "RUN_DIR_COMPONENTS",
]
