"""Compiler table, NetCDF detection and the resolved toolchain.

The compiler ``if/elif`` chain of the bash script is a data table
(:data:`COMPILERS`), and the two NetCDF flavours (``nf-config`` /
``nc-config``) are probed once.
"""

from __future__ import annotations

import os
import platform
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional, Tuple

from .errors import LOG, BuildError
from .shell import Shell


@dataclass(frozen=True)
class CompilerSpec:
    """One row of the compiler table."""

    name: str
    cpp: str
    fflags: str
    openmp: str = ""
    cft: Optional[str] = None  # compiler actually invoked (None: same as name)
    #: Match any name starting with ``name``: ``gfortran`` also accepts
    #: ``gfortran-12`` (the bash chain tested ``$FC == gfortran*``).
    prefix_match: bool = False

    def matches(self, fc: str) -> bool:
        if self.prefix_match:
            return fc.startswith(self.name)
        return fc == self.name

    @classmethod
    def find(cls, fc: str) -> Optional["CompilerSpec"]:
        """The table row for ``fc``, or ``None`` when it is unknown."""
        for spec in COMPILERS:
            if spec.matches(fc):
                return spec
        return None


COMPILERS: Tuple[CompilerSpec, ...] = (
    CompilerSpec(
        "ifort",
        "cpp -traditional -DLinux -DIfort",
        "-O2 -mcmodel=medium -fno-alias -i4 -r8 -fp-model precise",
    ),
    CompilerSpec(
        "ifc",
        "cpp -traditional -DLinux -DIfort",
        "-O2 -mcmodel=medium -fno-alias -i4 -r8 -fp-model precise",
        cft="ifort",
    ),
    CompilerSpec(
        "ifx",
        "cpp -traditional -DLinux -DIfort",
        "-O2 -mcmodel=medium -fno-alias -i4 -r8 -fp-model precise",
        openmp="-qopenmp",
    ),
    CompilerSpec(
        "gfortran",
        "cpp -traditional -DLinux",
        "-O2 -mcmodel=medium -fdefault-real-8 -fdefault-double-8 -std=legacy",
        openmp="-fopenmp",
        prefix_match=True,
    ),
    CompilerSpec(
        "pgfortran",
        "cpp -traditional -DLinux -DXLF",
        "-g -fast -r8 -i4 -mcmodel=medium -Mbackslash",
        openmp="-openmp",
    ),
    CompilerSpec(
        "nvfortran",
        "cpp -traditional -DLinux -DXLF",
        "-g -fast -r8 -i4 -mcmodel=medium -Mbackslash",
        openmp="-openmp",
    ),
)


@dataclass
class Paths:
    """The directories one build works with."""

    run_dir: Path  # where croco/partit/ncjoin end up (``pwd``)
    source: Path  # CROCO sources
    root: Path  # ``source/..``: PISCES, AGRIF, MUSTANG, ...
    scratch: Path  # ``./Compile``
    agrif: Path  # ``root/AGRIF``

    @classmethod
    def create(cls, source: Path, scratch: Path) -> "Paths":
        source = Path(source)
        return cls(
            run_dir=Path.cwd(),
            source=source,
            root=source.parent,
            scratch=Path(scratch),
            agrif=source.parent / "AGRIF",
        )


@dataclass(frozen=True)
class ToolchainSettings:
    """The values :class:`Toolchain` needs, already resolved."""

    fc: str
    fc_from_env: bool
    mpif90: str
    fflags_override: Optional[str]
    nprocs: int
    netcdf_lib: str
    netcdf_inc: str
    netcdf_lib_default: str
    netcdf_inc_default: str
    cc_override: Optional[str] = None
    cflags_override: Optional[str] = None
    make_override: Optional[str] = None
    mpilib: str = ""
    mpiinc: str = ""
    environ: Optional[Mapping[str, str]] = None


@dataclass
class Toolchain:
    """Everything the build needs to know about compilers and flags.

    ``append`` adds a non-empty fragment, separated by single blanks, so the
    flags stay readable in ``Makedefs`` and in the ``make`` command line.
    """

    os: str  # "Linux" | "Darwin"
    fc: str  # compiler requested by the user
    cft: str  # compiler actually invoked (may become mpif90)
    mpif90: str
    mpilib: str = ""
    mpiinc: str = ""
    cpp: str = ""
    fflags: str = ""
    ldflags: str = ""
    cppflags: str = ""
    cc: str = ""
    cflags: str = ""
    make: str = "make"
    nprocs: int = 1

    # -- flag accumulation -------------------------------------------------
    def append(self, attribute: str, flags: str) -> None:
        """Append ``flags`` to ``attribute`` when it is not empty."""
        if flags:
            current = getattr(self, attribute)
            setattr(self, attribute, " ".join((current, flags)) if current
                    else flags)

    def add_fflags(self, flags: str) -> None:
        self.append("fflags", flags)

    def add_cppflags(self, flags: str) -> None:
        self.append("cppflags", flags)

    def add_ldflags(self, flags: str) -> None:
        self.append("ldflags", flags)

    # -- make probes -------------------------------------------------------
    def make_variable(self, shell: Shell, name: str, cwd=None) -> str:
        """Ask make for the value of a variable (``CC``, ``CFLAGS``, ...)."""
        snippet = f"dummy_target:\n\t@echo $({name})\n"
        return shell.capture(
            [self.make, "-f", "-", "dummy_target"],
            cwd=cwd,
            text_input=snippet,
        ).strip()

    def is_gnu_c(self, shell: Shell, cc: Optional[str] = None) -> bool:
        """Is ``cc`` (or the resolved C compiler) the GNU C compiler?"""
        compiler = cc or self.cc
        if not compiler:
            return False
        proc = shell.probe([compiler, "-v"])
        return "gcc version" in (proc.stdout + proc.stderr)

    # -- construction ------------------------------------------------------
    @classmethod
    def from_settings(
        cls, settings: ToolchainSettings, os_name: str
    ) -> "Toolchain":
        """Build the toolchain, reproducing the bash branch selection."""
        environ = os.environ if settings.environ is None else settings.environ
        if os_name not in ("Linux", "Darwin"):
            raise BuildError(
                f"Unknown Operating System, please add it to jobcomp "
                f"({os_name})"
            )
        cft, cpp, fflags = cls._compiler_flags(
            settings.fc, settings.fc_from_env
        )
        toolchain = cls(
            os=os_name,
            fc=settings.fc,
            cft=cft,
            mpif90=settings.mpif90,
            mpilib=environ.get("CROCO_MPILIB") or settings.mpilib,
            mpiinc=environ.get("CROCO_MPIINC") or settings.mpiinc,
            cpp=cpp,
            fflags=fflags,
            ldflags=settings.netcdf_lib,
            cppflags=f"{settings.netcdf_inc} -ICROCOFILES/AGRIF_INC".strip(),
            cc=settings.cc_override or "",
            cflags=settings.cflags_override or "",
            make=settings.make_override or cls.default_make(),
            nprocs=settings.nprocs,
        )
        if settings.fflags_override:
            toolchain.fflags = settings.fflags_override
        # NetCDF modules are needed by the Fortran compiler itself.  The bash
        # script uses the *detected* NETCDFINC here, not the user override;
        # kept for bug-for-bug compatibility.
        toolchain.add_fflags(settings.netcdf_inc_default)
        return toolchain

    @classmethod
    def _compiler_flags(
        cls, fc: str, fc_from_env: bool
    ) -> Tuple[str, str, str]:
        """Return ``(cft, cpp, fflags)`` for a requested compiler."""
        spec = CompilerSpec.find(fc)
        if spec is None:
            # Historical behaviour: an unknown compiler is only fatal when it
            # comes from the hard-coded default; a user-provided one is
            # ignored.  Keep the two-branch warning of the bash script.
            LOG.warning("Unknown Fortran compiler, please add it to jobcomp")
            if not fc_from_env:
                raise BuildError(
                    "Unknown Fortran compiler, please add it to jobcomp"
                )
            LOG.warning("(Ignoring unknown Fortran compiler)")
            return fc, "", ""
        cft = spec.cft or fc
        cpp = spec.cpp
        if spec.name in ("pgfortran", "nvfortran"):
            if cls.hostname().startswith("jean-zay"):
                cpp = f"{cpp} -DJEANZAY"
        return cft, cpp, spec.fflags

    @classmethod
    def resolve(cls, *, shell: Shell, os_name: str, **kwargs) -> "Toolchain":
        """Compatibility wrapper: resolve from keyword values."""
        return cls.from_settings(ToolchainSettings(**kwargs), os_name=os_name)

    @staticmethod
    def default_fflags(fc: str) -> str:
        """The compiler-table FFLAGS for ``fc`` (``""`` when unknown)."""
        spec = CompilerSpec.find(fc)
        return spec.fflags if spec is not None else ""

    # -- environment probes ------------------------------------------------
    @staticmethod
    def hostname() -> str:
        return os.environ.get("HOSTNAME") or platform.node()

    @staticmethod
    def default_make() -> str:
        return "gmake" if shutil.which("gmake") else "make"

    @staticmethod
    def detect_netcdf(shell: Shell) -> Tuple[str, str]:
        """Return ``(lib, inc)`` flags, like the nf-config/nc-config logic."""
        lib = ""
        inc = ""
        if shutil.which("nf-config"):
            lib = shell.capture(["nf-config", "--flibs"]).strip()
            inc = "-I" + shell.capture(["nf-config", "--includedir"]).strip()
        if shutil.which("nc-config"):
            libdir = shell.capture(["nc-config", "--libdir"]).strip()
            incdir = shell.capture(["nc-config", "--includedir"]).strip()
            lib = f"{lib} -L{libdir} -lnetcdf".strip()
            inc = f"{inc} -I{incdir}".strip()
        return lib, inc
