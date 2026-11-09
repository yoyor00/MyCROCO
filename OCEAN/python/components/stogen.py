"""STOGEN: the stochastic generator sources and its namelist."""

from __future__ import annotations

from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec
from jobcomp.errors import LOG

from .base import Component


class StogenComponent(Component):
    """STOGEN (stochastic generator) sources and its namelist."""

    name = "stogen"

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        yield CopySpec.contents(ctx.paths.root / "STOGEN", ctx.paths.scratch)

    def run_dir_files(self, ctx: BuildContext) -> None:
        if (ctx.run_dir / "namelist_sto").exists():
            LOG.info("  file namelist_sto exist in Run directory")
            return
        source = ctx.paths.root / "STOGEN" / "namelist_sto"
        if not source.exists():
            if ctx.dry_run:
                LOG.info("would copy %s -> %s (missing in dry run)", source,
                         ctx.run_dir)
                return
            # The bash script's ``cp`` fails and the build goes on.
            LOG.warning("%s not found, namelist_sto left unchanged", source)
            return
        ctx.copier.file(source, ctx.run_dir)
        LOG.info("  file namelist_sto copied from source directory")
