"""PISCES: biogeochemistry sources, namelists and the p2z/p5z choice."""

from __future__ import annotations

from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopyMode, CopySpec
from jobcomp.errors import LOG

from .base import Component


class PiscesComponent(Component):
    """PISCES (biogeochemistry) sources and its namelists."""

    name = "pisces"

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        root = ctx.paths.root
        scratch = ctx.paths.scratch
        yield CopySpec.contents(root / "PISCES", scratch)
        yield CopySpec.contents(root / "PISCES" / "LIB", scratch)
        yield CopySpec.contents(root / "PISCES" / "SED", scratch)
        yield CopySpec.glob(root / "PISCES", ctx.run_dir, "kRGB61*")

    def run_dir_files(self, ctx: BuildContext) -> None:
        if (ctx.run_dir / "namelist_pisces_ref").exists():
            LOG.info("  file namelist_pisces exists in Run directory")
            return
        pisces = ctx.paths.root / "PISCES"
        if not pisces.is_dir():
            return
        for path in sorted(pisces.glob("namelist_pisces*")):
            ctx.copier.file(path, ctx.run_dir)
        LOG.info("  file namelist_pisces copied from source directory")
        for path in sorted((pisces / "SED").glob("namelist_sediment*")):
            ctx.copier.file(path, ctx.run_dir)
        LOG.info("  file namelist_sediment copied from source directory")

    def prepare(self, ctx: BuildContext) -> None:
        for kind, text in (
            ("p2z", "Simple version of PISCES is defined"),
            ("p5z", "Quota version of PISCES is defined"),
        ):
            if ctx.detector.has(f"{kind}isdefined"):
                LOG.info("%s: copy the appropriate namelist", text)
                self.install_namelist(ctx, kind)

    @staticmethod
    def install_namelist(ctx: BuildContext, kind: str) -> None:
        """Copy ``namelist_pisces_cfg_<kind>`` to ``namelist_pisces_cfg``."""
        source = ctx.run_dir / f"namelist_pisces_cfg_{kind}"
        if not source.exists() and not ctx.dry_run:
            LOG.warning("%s not found, PISCES namelist left unchanged", source)
            return
        for name in ("namelist_pisces_cfg", "namelist_pisces_cfg.1"):
            # The chosen version wins over an older file of the other version.
            ctx.copier.file(
                source,
                ctx.run_dir / name,
                overwrite=True,
                mode=CopyMode.ALWAYS,
            )
