"""Command line, configuration and entry point of the CROCO build driver.

``Jobcomp`` is the whole program: it parses the options, resolves the
toolchain and the directories, and hands over to :class:`jobcomp.build.Build`.
``main`` is only the entry point used by ``OCEAN/jobcomp.py``.
"""

from __future__ import annotations

import argparse
import logging
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from .build import Build
from .config import Config
from .context import BuildContext
from .copy import CopyMode
from .errors import LOG, BuildError
from .features import FeatureDetector
from .parameters import BuildParameters
from .shell import Shell
from .toolchain import CompilerSpec, Paths, Toolchain

HELP = """-----------------------------------------------------------------
JOBCOMP script to compile CROCO on OS = Linux or Darwin (Mac OS)
-----------------------------------------------------------------
Usage: jobcomp [OPTIONS]
       jobcomp --sh [OPTIONS]    (run the historical bash script jobcomp.sh)

Options:
  -h, --help                  Show this help message and exit.
  --sh                        Hand over to the original bash implementation
                              (jobcomp.sh, shipped unmodified).  All other
                              options are passed through to it untouched.
  --src <path>                Specify the CROCO source directory
                              (default: ../croco/OCEAN).
  --fc <compiler>             Specify the Fortran compiler
                              (default: gfortran).  Available = gfortran,
                              ifort, ifc, ifx, nvfortran, pgfortran
  --mpif90 <compiler>         Specify the MPI Fortran compiler
                              (default: mpif90).
  --fflags <flags>            Specify the compilation flags
                              (default: if fc=ifort or ifx:
                               '-O2 -mcmodel=medium -fno-alias -i4 -r8
                                -fp-model precise' ).
                              (default: if fc=gfortran:
                               '-O2 -mcmodel=medium -fdefault-real-8
                                -fdefault-double-8 -std=legacy' ).
                              (default: if fc=nvfortran or pgfortran:
                               '-g -fast -r8 -i4 -mcmodel=medium
                                -Mbackslash' ).
  --jobs <n>                  Number of processes for compilation (default: 1)
  --netcdf-inc <path>         Specify the NetCDF include directory.
                              (default: nf-config --includedir)
  --netcdf-lib <path>         Specify the NetCDF library directory.
                              (default: nf-config --flibs)
  --xios-dir <path>           Specify the XIOS root directory.
                              (default: $HOME/xios)
  --prism-dir <path>          Specify the OASIS-MCT root directory.
                              (default: ../../../oasis3-mct/compile_oa3-mct)

Incremental build options:
  --copy-mode <mode>          newer (default) | always.  "newer" only copies
                              a file when it is missing or its timestamp is
                              newer than the copy in Compile/.
  --keep-scratch              Do not wipe ./Compile before building.
  --clean                     Wipe ./Compile before building (default).
  --delete-stale              With --keep-scratch, remove files that no longer
                              exist in the source tree.

Advanced options:
  --cc <compiler>             C compiler (default: value of $(CC) in make)
  --cflags <flags>            C flags (appended to the make value)
  --make <command>            Make command (default: gmake if available)
  --dry-run                   Print what would be done, change nothing.
  -v, --verbose               Verbose output (repeat for more).
  -q, --quiet                 Only report errors.
  --log-file <path>           Also write the log to this file.

Configuration file:
  A ``jobcomp.conf`` file next to this script is read before the environment
  variables are looked at, using the historical names: SOURCE1, FC, MPIF90,
  MPILIB, MPIINC, PRISM_ROOT_DIR, XIOS_ROOT_DIR (one KEY=VALUE per line).
  Precedence: command line > CROCO_* environment > jobcomp.conf > defaults.
"""


class _ArgumentParser(argparse.ArgumentParser):
    """argparse with the historical "Unknown option" message and exit 1."""

    def error(self, message: str):  # pragma: no cover - argparse plumbing
        text = message.replace("unrecognized arguments: ", "Unknown option: ")
        sys.stderr.write(f"{text}\n")
        sys.stderr.write(HELP)
        raise SystemExit(1)


class Jobcomp:
    """The driver: one instance per invocation."""

    #: Options that take a value (``--fflags -O2`` must work like in bash).
    VALUE_OPTIONS = (
        "--src",
        "--fc",
        "--mpif90",
        "--fflags",
        "--jobs",
        "--netcdf-inc",
        "--netcdf-lib",
        "--xios-dir",
        "--prism-dir",
        "--copy-mode",
        "--cc",
        "--cflags",
        "--make",
        "--log-file",
    )

    #: Options that take no value (never glued onto the previous option).
    FLAGS = (
        "-h",
        "--help",
        "--sh",
        "--print-fflags",
        "-v",
        "--verbose",
        "-q",
        "--quiet",
        "--keep-scratch",
        "--clean",
        "--delete-stale",
        "--dry-run",
    )

    def __init__(
        self, argv: Optional[Sequence[str]] = None, script_path=None
    ) -> None:
        self.argv = list(sys.argv[1:] if argv is None else argv)
        self.config = Config(script_path)

    # -- entry point -------------------------------------------------------
    def run(self) -> int:
        """Do whatever the command line asks for, and return an exit code."""
        passthrough = self.legacy_argv()
        if passthrough is not None:
            try:
                return self.run_legacy(passthrough)
            except BuildError as exc:
                sys.stderr.write(f"{exc}\n")
                return exc.returncode or 1
        args = self.parse_args()
        if args.help:
            sys.stdout.write(HELP)
            return 0
        if args.print_fflags:
            sys.stdout.write(self.default_fflags(args) + "\n")
            return 0
        self.setup_logging(args)
        try:
            Build(self.create_context(args)).run()
        except BuildError as exc:
            LOG.error("%s", exc)
            return exc.returncode or 1
        except OSError as exc:
            LOG.error("%s", exc)
            return 1
        except KeyboardInterrupt:  # pragma: no cover
            LOG.error("interrupted")
            return 130
        return 0

    # -- the historical bash script ---------------------------------------
    def legacy_argv(self) -> Optional[List[str]]:
        """The arguments to hand to the bash ``jobcomp``, or ``None``.

        ``--sh`` is looked for before ``argparse`` runs so that options only
        the bash script understands reach it untouched instead of being
        rejected.  A ``--sh`` where a value is expected (``--fflags --sh``)
        is a value, not a request to hand over.
        """
        passthrough: List[str] = []
        expect_value = False
        found = False
        for token in self.argv:
            if expect_value:
                passthrough.append(token)
                expect_value = False
                continue
            if token == "--sh":
                found = True
                continue
            passthrough.append(token)
            expect_value = token in self.VALUE_OPTIONS
        return passthrough if found else None

    def run_legacy(self, argv: Sequence[str]) -> int:
        """Exec the historical bash script, keeping its output and status."""
        command = ["bash", str(self.config.legacy_script()), *argv]
        if os.name != "posix":  # pragma: no cover - the bash script is POSIX
            return subprocess.call(command)
        os.execvp(command[0], command)
        return 127  # pragma: no cover - execvp only returns on failure

    # -- options -----------------------------------------------------------
    def parse_args(self) -> argparse.Namespace:
        argv = self.glue_flag_values(self.argv)
        return self.build_parser().parse_args(argv)

    def build_parser(self) -> _ArgumentParser:
        parser = _ArgumentParser(
            prog="jobcomp", add_help=False, allow_abbrev=False
        )
        parser.add_argument("-h", "--help", action="store_true", dest="help")
        parser.add_argument("--sh", action="store_true")
        # plumbing for the scripts: print the FFLAGS the toolchain would use
        parser.add_argument("--print-fflags", action="store_true")
        parser.add_argument("--src")
        parser.add_argument("--fc")
        parser.add_argument("--mpif90")
        parser.add_argument("--fflags")
        parser.add_argument("--jobs", type=int)
        parser.add_argument("--netcdf-inc")
        parser.add_argument("--netcdf-lib")
        parser.add_argument("--xios-dir")
        parser.add_argument("--prism-dir")
        parser.add_argument(
            "--copy-mode",
            choices=[mode.value for mode in CopyMode],
            default=None,
        )
        parser.add_argument("--keep-scratch", action="store_true")
        parser.add_argument("--clean", action="store_true")
        parser.add_argument("--delete-stale", action="store_true")
        parser.add_argument("--cc")
        parser.add_argument("--cflags")
        parser.add_argument("--make")
        parser.add_argument("--dry-run", action="store_true")
        parser.add_argument("-v", "--verbose", action="count", default=0)
        parser.add_argument("-q", "--quiet", action="store_true")
        parser.add_argument("--log-file")
        return parser

    def glue_flag_values(self, argv: Sequence[str]) -> List[str]:
        """Allow ``--fflags -O2 -g`` like the bash script did.

        ``argparse`` refuses values that start with a dash, while ``jobcomp``
        always took the next word as the value (``FFLAGS1="$2"``).  For the
        value-taking options we therefore rewrite ``--fflags -O2`` into
        ``--fflags=-O2`` unless the next word is another known option.
        """
        out: List[str] = []
        index = 0
        while index < len(argv):
            token = argv[index]
            nxt = argv[index + 1] if index + 1 < len(argv) else None
            if (
                token in self.VALUE_OPTIONS
                and nxt is not None
                and nxt not in self.VALUE_OPTIONS
                and nxt not in self.FLAGS
            ):
                out.append(f"{token}={nxt}")
                index += 2
                continue
            out.append(token)
            index += 1
        return out

    def setup_logging(self, args: argparse.Namespace) -> None:
        level = logging.INFO
        if args.quiet:
            level = logging.WARNING
        elif args.verbose:
            level = logging.DEBUG
        LOG.setLevel(level)
        LOG.handlers.clear()
        stream = logging.StreamHandler(sys.stdout)
        stream.setFormatter(logging.Formatter("%(message)s"))
        LOG.addHandler(stream)
        if args.log_file:
            handler = logging.FileHandler(args.log_file, mode="w")
            handler.setFormatter(
                logging.Formatter("%(asctime)s %(levelname)s %(message)s")
            )
            LOG.addHandler(handler)
        LOG.propagate = False

    # -- settings ----------------------------------------------------------
    def parameters(self) -> BuildParameters:
        """The resolved build parameters, without copying or compiling.

        This is the entry point for another build system: the same flags and
        paths a build would use, as data (see :mod:`jobcomp.parameters`).
        """
        return BuildParameters.resolve(self.create_context(self.parse_args()))

    def default_fflags(self, args: argparse.Namespace) -> str:
        """The ``FFLAGS`` a build with these options would start from.

        Used by the coupling scripts that used to run
        ``sed -e "s|-O3|-O2|g"`` on the bash script to weaken the optimisation
        before compiling.
        """
        fc = self.config.pick(args.fc, "CROCO_CFT1")
        spec = CompilerSpec.find(fc)
        fflags = spec.fflags if spec is not None else ""
        override = self.config.pick(args.fflags, "CROCO_FFLAGS1", "")
        return override or fflags

    def create_context(self, args: argparse.Namespace) -> BuildContext:
        """Turn the parsed options into everything a build needs."""
        os_name = platform.system()
        LOG.info("OPERATING SYSTEM IS: %s", os_name)
        shell = Shell(dry_run=args.dry_run)
        tools = self.toolchain(args, shell, os_name)
        paths = self.paths(args)
        return BuildContext(
            paths=paths,
            tools=tools,
            detector=FeatureDetector(
                shell, tools.cpp, paths.scratch, paths.source, paths.run_dir
            ),
            shell=shell,
            copy_mode=CopyMode(
                self.config.pick(args.copy_mode, "CROCO_COPY_MODE",
                                 CopyMode.NEWER.value)
            ),
            dry_run=args.dry_run,
            keep_scratch=args.keep_scratch and not args.clean,
            delete_stale=args.delete_stale,
            xios_dir=self.config.expand(
                self.config.pick(args.xios_dir, "CROCO_XIOS_ROOT_DIR")
            ),
            prism_dir=self.config.expand(
                self.config.pick(args.prism_dir, "CROCO_PRISM_ROOT_DIR")
            ),
        )

    def toolchain(
        self, args: argparse.Namespace, shell: Shell, os_name: str
    ) -> Toolchain:
        """Resolve the compilers, the flags and the NetCDF paths."""
        config = self.config
        detected_lib, detected_inc = Toolchain.detect_netcdf(shell)
        fc_from_env = args.fc is not None or bool(
            config.environ.get("CROCO_CFT1")
        )
        tools = Toolchain.resolve(
            shell=shell,
            os_name=os_name,
            fc=config.pick(args.fc, "CROCO_CFT1"),
            fc_from_env=fc_from_env,
            mpif90=config.pick(args.mpif90, "CROCO_MPIF90"),
            fflags_override=config.pick(args.fflags, "CROCO_FFLAGS1", "")
            or None,
            nprocs=int(config.pick(args.jobs, "NPROCS", "1")),
            netcdf_lib=config.pick(
                args.netcdf_lib, "CROCO_NETCDFLIB", detected_lib
            ),
            netcdf_inc=config.pick(
                args.netcdf_inc, "CROCO_NETCDFINC", detected_inc
            ),
            netcdf_lib_default=detected_lib,
            netcdf_inc_default=detected_inc,
            cc_override=args.cc,
            cflags_override=args.cflags,
            make_override=args.make,
        )
        if tools.nprocs < 1:
            raise BuildError(f"invalid number of jobs: {tools.nprocs}")
        return tools

    def paths(self, args: argparse.Namespace) -> Paths:
        """Where the sources are, and the ``./Compile`` scratch directory."""
        paths = Paths.create(
            self.config.expand(self.config.pick(args.src, "CROCO_SRC")),
            Path("./Compile"),
        )
        if not paths.source.is_dir():
            raise BuildError(
                f"CROCO source directory not found: {paths.source}"
            )
        return paths


def main(argv: Optional[Sequence[str]] = None, script_path=None) -> int:
    """Entry point used by ``OCEAN/jobcomp.py``."""
    return Jobcomp(argv, script_path).run()
