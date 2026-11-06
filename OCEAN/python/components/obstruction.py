"""OBSTRUCTION: the immersed-boundary sources."""

from __future__ import annotations

from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec

from .base import Component


class ObstructionComponent(Component):
    name = "obstruction"

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        yield CopySpec(src=ctx.paths.root / "OBSTRUCTION",
                       dst=ctx.paths.scratch, contents=True)
