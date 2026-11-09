"""Command line, configuration and entry point of the CROCO build driver.

``Jobcomp`` parses the options, resolves the settings, and hands over to
:class:`jobcomp.build.Build`.  ``main`` is only the entry point used by
``OCEAN/jobcomp.py``.
"""

from __future__ import annotations

import argparse
import logging
import os
import platform
import subprocess
import sys
from dataclasses import dataclass
from typing import Any, Callable, List, Optional, Sequence, Tuple

from .build import Build
from .config import Config
from .context import BuildContext
from .copy import CopyMode
from .errors import LOG, BuildError
from .parameters import BuildParameters
from .settings import BuildSettings
from .shell import Shell
from .toolchain import Toolchain

USAGE = (
    "%(prog)s [OPTIONS]\n"
    "       %(prog)s --sh [OPTIONS]    "
    "(run the historical bash script jobcomp.sh)"
)

HELP_DESCRIPTION = """\
-----------------------------------------------------------------
JOBCOMP script to compile CROCO on OS = Linux or Darwin (Mac OS)
-----------------------------------------------------------------"""

HELP_EPILOG = """\
Configuration file:
  A jobcomp.conf file next to this script is read before the environment
  variables are looked at, using the historical names: SOURCE1, FC, MPIF90,
  MPILIB, MPIINC, PRISM_ROOT_DIR, XIOS_ROOT_DIR (one KEY=VALUE per line).
  Precedence: command line > CROCO_* environment > jobcomp.conf > defaults."""

_UNSET = object()


@dataclass(frozen=True)
class Option:
    """One command-line option, described once for parser and help."""

    flags: Tuple[str, ...]
    dest: str
    help: str
    section: str = "Options"
    action: Optional[str] = None
    takes_value: bool = False
    value_type: Optional[Callable[[str], Any]] = None
    choices: Optional[Tuple[str, ...]] = None
    default: Any = _UNSET


OPTIONS: Tuple[Option, ...] = (
    Option(
        ("-h", "--help"),
        "help",
        "Show this help message and exit.",
        action="store_true",
    ),
    Option(
        ("--sh",),
        "sh",
        "Hand over to the original bash implementation (jobcomp.sh, shipped "
        "unmodified).  All other options are passed through to it untouched.",
        action="store_true",
    ),
    Option(
        ("--src",),
        "src",
        "Specify the CROCO source directory (default: ../croco/OCEAN).",
        takes_value=True,
    ),
    Option(
        ("--fc",),
        "fc",
        "Specify the Fortran compiler (default: gfortran).  Available = "
        "gfortran, ifort, ifc, ifx, nvfortran, pgfortran",
        takes_value=True,
    ),
    Option(
        ("--mpif90",),
        "mpif90",
        "Specify the MPI Fortran compiler (default: mpif90).",
        takes_value=True,
    ),
    Option(
        ("--fflags",),
        "fflags",
        "Specify the compilation flags, overriding the compiler defaults.  "
        "A value may start with a dash.",
        takes_value=True,
    ),
    Option(
        ("--jobs",),
        "jobs",
        "Number of processes for compilation (default: 1).",
        takes_value=True,
        value_type=int,
    ),
    Option(
        ("--netcdf-inc",),
        "netcdf_inc",
        "Specify the NetCDF include directory (default: nf-config "
        "--includedir).",
        takes_value=True,
    ),
    Option(
        ("--netcdf-lib",),
        "netcdf_lib",
        "Specify the NetCDF library directory (default: nf-config --flibs).",
        takes_value=True,
    ),
    Option(
        ("--xios-dir",),
        "xios_dir",
        "Specify the XIOS root directory (default: $HOME/xios).",
        takes_value=True,
    ),
    Option(
        ("--prism-dir",),
        "prism_dir",
        "Specify the OASIS-MCT root directory (default: "
        "../../../oasis3-mct/compile_oa3-mct).",
        takes_value=True,
    ),
    Option(
        ("--copy-mode",),
        "copy_mode",
        "newer (default) | always.  \"newer\" only copies a file when it is "
        "missing or its timestamp is newer than the copy in Compile/.",
        section="Incremental build options",
        takes_value=True,
        choices=tuple(mode.value for mode in CopyMode),
        default=None,
    ),
    Option(
        ("--keep-scratch",),
        "keep_scratch",
        "Do not wipe ./Compile before building.",
        section="Incremental build options",
        action="store_true",
    ),
    Option(
        ("--clean",),
        "clean",
        "Wipe ./Compile before building (default).",
        section="Incremental build options",
        action="store_true",
    ),
    Option(
        ("--delete-stale",),
        "delete_stale",
        "With --keep-scratch, remove files that no longer exist in the "
        "source tree.",
        section="Incremental build options",
        action="store_true",
    ),
    Option(
        ("--cc",),
        "cc",
        "C compiler (default: value of $(CC) in make).",
        section="Advanced options",
        takes_value=True,
    ),
    Option(
        ("--cflags",),
        "cflags",
        "C flags (appended to the make value).",
        section="Advanced options",
        takes_value=True,
    ),
    Option(
        ("--make",),
        "make",
        "Make command (default: gmake if available).",
        section="Advanced options",
        takes_value=True,
    ),
    Option(
        ("--print-fflags",),
        "print_fflags",
        "Print the FFLAGS the toolchain would use, then exit.",
        section="Advanced options",
        action="store_true",
    ),
    Option(
        ("--dry-run",),
        "dry_run",
        "Print what would be done, change nothing.",
        section="Advanced options",
        action="store_true",
    ),
    Option(
        ("-v", "--verbose"),
        "verbose",
        "Verbose output (repeat for more).",
        section="Advanced options",
        action="count",
        default=0,
    ),
    Option(
        ("-q", "--quiet"),
        "quiet",
        "Only report errors.",
        section="Advanced options",
        action="store_true",
    ),
    Option(
        ("--log-file",),
        "log_file",
        "Also write the log to this file.",
        section="Advanced options",
        takes_value=True,
    ),
)

#: Options that take a value (``--fflags -O2`` must work like in bash).
VALUE_OPTIONS = tuple(
    flag
    for option in OPTIONS
    if option.takes_value
    for flag in option.flags
    if flag.startswith("--")
)

#: Options that take no value (never glued onto the previous option).
FLAGS = tuple(
    flag
    for option in OPTIONS
    if not option.takes_value
    for flag in option.flags
)


class _ArgumentParser(argparse.ArgumentParser):
    """argparse with the historical "Unknown option" message and exit 1."""

    def error(self, message: str):  # pragma: no cover - argparse plumbing
        text = message.replace("unrecognized arguments: ", "Unknown option: ")
        sys.stderr.write(f"{text}\n")
        sys.stderr.write(self.format_help())
        raise SystemExit(1)


class Jobcomp:
    """The driver: one instance per invocation."""

    VALUE_OPTIONS = VALUE_OPTIONS
    FLAGS = FLAGS

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
            sys.stdout.write(self.build_parser().format_help())
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
            expect_value = token in VALUE_OPTIONS
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
            prog="jobcomp",
            add_help=False,
            allow_abbrev=False,
            usage=USAGE,
            description=HELP_DESCRIPTION,
            epilog=HELP_EPILOG,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        groups = {}
        for option in OPTIONS:
            group = groups.get(option.section)
            if group is None:
                group = parser.add_argument_group(option.section)
                groups[option.section] = group
            kwargs = {"dest": option.dest, "help": option.help}
            if option.action is not None:
                kwargs["action"] = option.action
            if option.value_type is not None:
                kwargs["type"] = option.value_type
            if option.choices is not None:
                kwargs["choices"] = option.choices
            if option.default is not _UNSET:
                kwargs["default"] = option.default
            group.add_argument(*option.flags, **kwargs)
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
                token in VALUE_OPTIONS
                and nxt is not None
                and nxt not in VALUE_OPTIONS
                and nxt not in FLAGS
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
        override = self.config.pick(args.fflags, "CROCO_FFLAGS1", "")
        return override or Toolchain.default_fflags(fc)

    def create_context(self, args: argparse.Namespace) -> BuildContext:
        """Turn the parsed options into everything a build needs."""
        os_name = platform.system()
        LOG.info("OPERATING SYSTEM IS: %s", os_name)
        shell = Shell(dry_run=args.dry_run)
        settings = BuildSettings.resolve(self.config, args, shell, os_name)
        if not settings.source.is_dir():
            raise BuildError(
                f"CROCO source directory not found: {settings.source}"
            )
        return BuildContext.from_settings(settings, shell)


def main(argv: Optional[Sequence[str]] = None, script_path=None) -> int:
    """Entry point used by ``OCEAN/jobcomp.py``."""
    return Jobcomp(argv, script_path).run()
