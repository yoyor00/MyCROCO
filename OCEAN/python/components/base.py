"""``Component``: the protocol every build feature implements.

A component declares what it copies (``copy_specs``), which files it drops in
the run directory (``run_dir_files``), which flags it adds (``configure``) and
what has to happen before the build (``prepare``).  The order in which the
components are asked is in :mod:`components.pipelines`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Optional

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec
from jobcomp.errors import LOG


class Component:
    """One independent feature of the build.

    Subclasses override the hooks they need; ``detect`` performs the
    historical "Checking COMPILE<label>..." dance and caches the answer in
    the context.
    """

    name = "component"
    key: Optional[str] = None  # testkeys.F marker (lower case)
    label: Optional[str] = None
    announce = False  # does the bash script print " => <label> activated"?

    def detect(self, ctx: BuildContext) -> bool:
        """Detect this feature once and cache the answer in the context."""
        if self.key is None:
            ctx.features[self.name] = True
            return True
        if self.label:
            LOG.info("Checking COMPILE%s...", self.label)
        active = ctx.detector.has(self.key)
        if active and self.announce and self.label:
            LOG.info(" => %s activated", self.label)
        ctx.features[self.name] = active
        return active

    def enabled(self, ctx: BuildContext) -> bool:
        """Was this feature detected?  Detection is lazy if not done yet."""
        if self.name in ctx.features:
            return ctx.feature(self.name)
        return self.detect(ctx)

    # -- hooks -------------------------------------------------------------
    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        return ()

    def run_dir_files(self, ctx: BuildContext) -> None:
        """Namelists / static files to install next to the executable."""

    def configure(self, ctx: BuildContext) -> None:
        """Add flags, based on the detected features."""

    def prepare(self, ctx: BuildContext) -> None:
        """Filesystem / pre-build work that is not a simple copy."""

    # -- helpers -----------------------------------------------------------
    def symlink(
        self, ctx: BuildContext, target: Path, directory: Path
    ) -> None:
        """Force ``directory/<name>`` to point at ``target``."""
        link = Path(directory) / Path(target).name
        if ctx.dry_run:
            LOG.info("would symlink %s -> %s", link, target)
            return
        try:
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(target)
        except OSError as exc:
            # Missing XIOS install: the bash script keeps going as well.
            LOG.warning("cannot create symlink %s: %s", link, exc)
