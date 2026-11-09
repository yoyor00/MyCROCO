"""The resolved build parameters, for build systems other than this one.

``jobcomp`` can hand its configuration to another build system (scons, cmake,
a Makefile of your own) without running anything::

    from jobcomp.jobcomp_driver import Jobcomp

    params = Jobcomp(["--fc", "gfortran", "--jobs", "8"]).parameters()
    print(params.cc, params.fflags)
    print(params.include_dirs, params.library_dirs, params.libraries)

Everything is resolved exactly as for a build — options, then ``CROCO_*``
environment variables, then ``jobcomp.conf``, then the defaults, plus the
feature-dependent flags of AGRIF, MPI, XIOS, OASIS, OpenMP and OpenACC — and
nothing is written: copying the sources, writing ``Makedefs`` and compiling
stay the business of ``jobcomp`` itself (``Build``).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Tuple

from components import (
    AGRIF,
    configure_after,
    configure_before,
    detect_features,
)

from .build import EXECUTABLES
from .context import BuildContext
from .flags import classify_flags as _classify
from .makedefs import Makedefs


@dataclass(frozen=True)
class BuildParameters:
    """What a build needs to know, as plain data (see the module docstring)."""

    source: Path
    scratch: Path
    run_dir: Path
    make: str
    nprocs: int
    cpp: str
    agrif_cppflags: str  # extra flags of the AGRIF library (MPI variant)
    fc: str  # Fortran compiler, already switched to MPIF90 when MPI is on
    cc: str  # C compiler of the AGRIF library ("" when not known / not used)
    cflags: str
    fflags: str
    cppflags: str
    ldflags: str
    features: Mapping[str, bool]
    defines: Tuple[str, ...]
    include_dirs: Tuple[str, ...]
    library_dirs: Tuple[str, ...]
    libraries: Tuple[str, ...]
    executables: Tuple[str, ...]
    makedefs: Mapping[str, str]

    @classmethod
    def resolve(cls, ctx: BuildContext) -> "BuildParameters":
        """Resolve the flags of ``ctx`` without copying or compiling."""
        detect_features(ctx)
        configure_before(ctx)
        if ctx.feature("agrif"):
            # only the AGRIF library is built from C
            AGRIF.resolve_c_compiler(ctx)
        configure_after(ctx)
        tools = ctx.tools
        defines, includes, library_dirs, libraries = _classify(
            tools.cpp, tools.cppflags, tools.fflags, tools.ldflags
        )
        return cls(
            source=ctx.paths.source,
            scratch=ctx.paths.scratch,
            run_dir=ctx.paths.run_dir,
            make=tools.make,
            nprocs=tools.nprocs,
            cpp=tools.cpp,
            agrif_cppflags=(
                AGRIF.cppflags(ctx) if ctx.feature("agrif") else ""
            ),
            fc=tools.cft,
            cc=tools.cc,
            cflags=tools.cflags,
            fflags=tools.fflags,
            cppflags=tools.cppflags,
            ldflags=tools.ldflags,
            features=dict(ctx.features),
            defines=defines,
            include_dirs=includes,
            library_dirs=library_dirs,
            libraries=libraries,
            executables=EXECUTABLES,
            makedefs=Makedefs(tools).variables(),
        )
