"""Resolved build settings: one object instead of a dozen arguments.

The command line, the ``CROCO_*`` environment, ``jobcomp.conf`` and the
defaults are combined here, once.  ``Toolchain``, ``BuildContext`` and
``BuildParameters`` then consume a :class:`BuildSettings` instead of each
repeating the precedence rules.
"""

from __future__ import annotations

import platform
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .config import Config
from .copy import CopyMode
from .errors import BuildError
from .shell import Shell
from .toolchain import Toolchain, ToolchainSettings


@dataclass(frozen=True)
class BuildSettings:
    """Everything that is resolved before a build context is created."""

    os_name: str
    source: Path
    scratch: Path
    tools: ToolchainSettings
    copy_mode: CopyMode
    keep_scratch: bool
    delete_stale: bool
    dry_run: bool
    xios_dir: Path
    prism_dir: Path

    @classmethod
    def resolve(
        cls,
        config: Config,
        args,
        shell: Shell,
        os_name: Optional[str] = None,
    ) -> "BuildSettings":
        """Apply ``CLI > environment > jobcomp.conf > defaults`` once."""
        os_name = platform.system() if os_name is None else os_name
        netcdf_lib, netcdf_inc = Toolchain.detect_netcdf(shell)

        nprocs_value = config.pick(args.jobs, "NPROCS", "1")
        try:
            nprocs = int(nprocs_value)
        except (TypeError, ValueError):
            raise BuildError(f"invalid number of jobs: {nprocs_value}")
        if nprocs < 1:
            raise BuildError(f"invalid number of jobs: {nprocs}")

        mode_value = config.pick(
            args.copy_mode, "CROCO_COPY_MODE", CopyMode.NEWER.value
        )
        try:
            copy_mode = CopyMode(mode_value)
        except ValueError:
            raise BuildError(f"invalid copy mode: {mode_value}")

        tools = ToolchainSettings(
            fc=config.pick(args.fc, "CROCO_CFT1"),
            fc_from_env=args.fc is not None
            or bool(config.environ.get("CROCO_CFT1")),
            mpif90=config.pick(args.mpif90, "CROCO_MPIF90"),
            fflags_override=config.pick(
                args.fflags, "CROCO_FFLAGS1", ""
            )
            or None,
            nprocs=nprocs,
            netcdf_lib=config.pick(
                args.netcdf_lib, "CROCO_NETCDFLIB", netcdf_lib
            ),
            netcdf_inc=config.pick(
                args.netcdf_inc, "CROCO_NETCDFINC", netcdf_inc
            ),
            netcdf_lib_default=netcdf_lib,
            netcdf_inc_default=netcdf_inc,
            cc_override=args.cc,
            cflags_override=args.cflags,
            make_override=args.make,
        )
        return cls(
            os_name=os_name,
            source=config.expand(config.pick(args.src, "CROCO_SRC")),
            scratch=Path("./Compile"),
            tools=tools,
            copy_mode=copy_mode,
            keep_scratch=args.keep_scratch and not args.clean,
            delete_stale=args.delete_stale,
            dry_run=args.dry_run,
            xios_dir=config.expand(
                config.pick(args.xios_dir, "CROCO_XIOS_ROOT_DIR")
            ),
            prism_dir=config.expand(
                config.pick(args.prism_dir, "CROCO_PRISM_ROOT_DIR")
            ),
        )
