"""The order in which the components are asked to do their work.

Keeping the lists explicit makes the ordering requirements (AGRIF library
before OpenMP flags, PISCES namelists before AGRIF, ...) visible instead of
implicit in the script line order.  Adding a feature is a class plus one entry
in the list that describes *when* it runs.
"""

from __future__ import annotations

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

#: Everything copied into the scratch directory (in this order).
COPY_COMPONENTS = (
    ModelSources(),
    AgrifComponent(),
    XiosComponent(),
    PiscesComponent(),
    MustangComponent(),
    ObstructionComponent(),
    StogenComponent(),
)

#: Components that drop a namelist in the run directory (never overwritten).
RUN_DIR_COMPONENTS = (StogenComponent(), PiscesComponent())

#: Flags that must be known before the AGRIF library is built.
CONFIGURE_BEFORE = (
    AgrifComponent(),
    MpiComponent(),
    XiosComponent(),
    OasisComponent(),
)

#: Pre-build work that touches the filesystem, in the order the bash script
#: did it: PISCES namelists, the AGRIF library and headers, the XIOS server.
PREPARE_COMPONENTS = (PiscesComponent(), AgrifComponent(), XiosComponent())

#: Flags that are only needed for the main compile.
CONFIGURE_AFTER = (OpenMpComponent(), OpenAccComponent())
