"""MUSTANG: the nesting tools copied into the scratch directory."""

from __future__ import annotations

from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec

from .base import Component


class MustangComponent(Component):
    name = "mustang"

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        yield CopySpec.contents(ctx.paths.root / "MUSTANG", ctx.paths.scratch)
