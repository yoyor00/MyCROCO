"""The build itself: clean, copy, configure, write build metadata, compile.

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
    configure_after,
    configure_before,
    copy_into_scratch,
    detect_features,
    install_run_dir_files,
    prepare as prepare_components,
)

from .context import BuildContext
from .copy import CopyMode, CopySpec
from .errors import LOG
from .makefile import read_source_lists
from .makedefs import Makedefs
from .manifest import MANIFEST_NAME, BuildManifest
from .settings import BuildBackend

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
        if self.ctx.backend is BuildBackend.SCONS:
            self.write_manifest()
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
        copy_into_scratch(ctx)
        install_run_dir_files(ctx)
        for pattern in OVERRIDES:
            ctx.copier.apply(
                CopySpec.glob(
                    ctx.run_dir,
                    ctx.paths.scratch,
                    pattern,
                    mode=CopyMode.ALWAYS,
                )
            )
        ctx.copier.finish()

    def configure(self) -> None:
        """Feature-dependent flags, before the AGRIF library is built."""
        detect_features(self.ctx)
        configure_before(self.ctx)

    def prepare(self) -> None:
        """Pre-build work (PISCES namelists, AGRIF library and headers)."""
        prepare_components(self.ctx)
        configure_after(self.ctx)

    def write_makedefs(self) -> None:
        Makedefs(self.ctx.tools).write(
            self.ctx.paths.scratch, dry_run=self.ctx.dry_run
        )

    def write_manifest(self) -> None:
        """Write the JSON configuration consumed by the SCons backend."""
        ctx = self.ctx
        path = ctx.paths.scratch / MANIFEST_NAME
        if ctx.dry_run:
            LOG.info("would write %s", path)
            return
        sources = read_source_lists(ctx.paths.scratch / "Makefile")
        manifest = BuildManifest.from_context(
            ctx, targets=list(EXECUTABLES), sources=sources
        )
        manifest.write(path)
        LOG.debug("wrote %s", path)

    def compile(self) -> None:
        """Run the selected compilation backend, then move the results out."""
        if self.ctx.backend is BuildBackend.SCONS:
            self.compile_scons()
        else:
            self.compile_make()
        for executable in EXECUTABLES:
            self.move_result(executable)

    def compile_make(self) -> None:
        """``make depend`` followed by ``make -j``."""
        ctx = self.ctx
        scratch = ctx.paths.scratch
        make = ctx.tools.make
        ctx.shell.run([make, "depend"], cwd=scratch)
        ctx.shell.run([make, "-j", str(ctx.tools.nprocs)], cwd=scratch)

    def compile_scons(self) -> None:
        """Run the packaged SCons build from the scratch directory."""
        from scons.backend import command, environment

        ctx = self.ctx
        argv = command(ctx, require=not ctx.dry_run)
        ctx.shell.run(
            argv, cwd=ctx.paths.scratch, env=environment(ctx)
        )

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
