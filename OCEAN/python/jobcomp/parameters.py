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

import shlex
from dataclasses import dataclass
from pathlib import Path
from typing import List, Mapping, Tuple

from components import CONFIGURE_AFTER, CONFIGURE_BEFORE
from components.agrif import AgrifComponent

from .build import EXECUTABLES
from .context import BuildContext
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
        for component in CONFIGURE_BEFORE:
            component.configure(ctx)
        if ctx.feature("agrif"):
            # only the AGRIF library is built from C
            AgrifComponent().resolve_c_compiler(ctx)
        for component in CONFIGURE_AFTER:
            component.configure(ctx)
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
                AgrifComponent().cppflags(ctx)
                if ctx.feature("agrif")
                else ""
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


def _classify(*flag_groups: str) -> Tuple[Tuple[str, ...], ...]:
    """Split flags into ``(defines, includes, -L, libraries)``, in order."""
    defines: List[str] = []
    includes: List[str] = []
    library_dirs: List[str] = []
    libraries: List[str] = []
    for group in flag_groups:
        tokens = shlex.split(group)
        index = 0
        while index < len(tokens):
            token = tokens[index]
            following = tokens[index + 1] if index + 1 < len(tokens) else None
            index += 1
            if token in ("-I", "-L", "-l", "-D") and following is not None:
                token, index = token + following, index + 1
            if token.startswith("-I"):
                includes.append(token[2:])
            elif token.startswith("-L"):
                library_dirs.append(token[2:])
            elif token.startswith("-l"):
                libraries.append(token[2:])
            elif token.startswith("-D"):
                defines.append(token[2:])
            elif token.endswith(".a"):
                libraries.append(token)
    return (
        tuple(dict.fromkeys(defines)),
        tuple(dict.fromkeys(includes)),
        tuple(dict.fromkeys(library_dirs)),
        tuple(dict.fromkeys(libraries)),
    )
