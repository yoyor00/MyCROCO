"""The build itself: clean, copy, configure, write ``Makedefs``, compile.

``Build`` holds the phase order in one place and takes the components from
:mod:`components`, so the sequence can be read without following a 500-line
script.  The two copy rules that are not features (the namelists dropped in
the run directory, and the run-directory files that win over the sources) live
here, next to the phase they belong to.
"""

from __future__ import annotations

import shutil
import sys

from components import (
    CONFIGURE_AFTER,
    CONFIGURE_BEFORE,
    COPY_COMPONENTS,
    PREPARE_COMPONENTS,
    RUN_DIR_COMPONENTS,
)

from .context import BuildContext
from .copy import CopyMode, CopySpec
from .errors import LOG
from .makedefs import Makedefs

#: Files the makefile leaves in the scratch directory for the run directory.
EXECUTABLES = ("croco", "partit", "ncjoin")

#: Run-directory files that win over the shared sources (``cppdefs.h``,
#: ``param.h``, a local ``Makefile``, ...).  Always copied: a local file wins
#: even when its timestamp is older than the source it replaces.
OVERRIDES = ("*.F90", "*.F", "*.h", "*.h90", "Make*")


class Build:
    """Runs one complete build described by a :class:`BuildContext`."""

    def __init__(self, ctx: BuildContext) -> None:
        self.ctx = ctx

    # -- phases ------------------------------------------------------------
    def run(self) -> None:
        self.clean_scratch()
        self.copy_sources()
        self.configure()
        self.prepare()
        self.write_makedefs()
        self.compile()
        self.finish()

    def clean_scratch(self) -> None:
        """Wipe ``./Compile`` unless the caller asked to keep it."""
        ctx = self.ctx
        scratch = ctx.paths.scratch
        if ctx.keep_scratch:
            LOG.info("Keeping scratch directory %s", scratch)
            if not ctx.dry_run:
                scratch.mkdir(parents=True, exist_ok=True)
            return
        if ctx.dry_run:
            LOG.info("would remove %s", scratch)
            return
        LOG.debug("cleaning %s", scratch)
        if scratch.is_dir() and not scratch.is_symlink():
            shutil.rmtree(scratch)
        elif scratch.exists():
            scratch.unlink()
        scratch.mkdir(parents=True)

    def copy_sources(self) -> None:
        """Copy every component in, then the run-directory overrides."""
        ctx = self.ctx
        for component in COPY_COMPONENTS:
            ctx.copy_many(component.copy_specs(ctx))
        for component in RUN_DIR_COMPONENTS:
            component.run_dir_files(ctx)
        for pattern in OVERRIDES:
            ctx.copier.apply(
                CopySpec(
                    src=ctx.run_dir,
                    dst=ctx.paths.scratch,
                    glob=pattern,
                    mode=CopyMode.ALWAYS,
                )
            )
        ctx.copier.finish()

    def configure(self) -> None:
        """Feature-dependent flags, before the AGRIF library is built."""
        for component in CONFIGURE_BEFORE:
            component.configure(self.ctx)

    def prepare(self) -> None:
        """Pre-build work (PISCES namelists, AGRIF library and headers)."""
        for component in PREPARE_COMPONENTS:
            component.prepare(self.ctx)
        for component in CONFIGURE_AFTER:
            component.configure(self.ctx)

    def write_makedefs(self) -> None:
        Makedefs(self.ctx.tools).write(
            self.ctx.paths.scratch, dry_run=self.ctx.dry_run
        )

    def compile(self) -> None:
        """``make depend``, ``make -j``, then move the results out."""
        ctx = self.ctx
        scratch = ctx.paths.scratch
        make = ctx.tools.make
        ctx.shell.run([make, "depend"], cwd=scratch)
        ctx.shell.run([make, "-j", str(ctx.tools.nprocs)], cwd=scratch)
        for executable in EXECUTABLES:
            self.move_result(executable)

    def finish(self) -> None:
        self._print_banner()
        LOG.info("CROCO is OK")

    # -- helpers -----------------------------------------------------------
    def move_result(self, executable: str) -> None:
        """Move ``Compile/<executable>`` to the run directory."""
        ctx = self.ctx
        path = ctx.paths.scratch / executable
        if ctx.dry_run:
            if path.exists():
                LOG.debug("would move %s to %s", path, ctx.run_dir)
            return
        if path.exists():
            shutil.move(str(path), str(ctx.run_dir / executable))

    def _print_banner(self) -> None:
        """Print ``croco_ascii.txt``, as the bash script does at the end."""
        banner = self.ctx.paths.scratch / "croco_ascii.txt"
        if banner.exists() and not self.ctx.dry_run:
            sys.stdout.write(banner.read_text(errors="replace"))
