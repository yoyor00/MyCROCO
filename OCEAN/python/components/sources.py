"""The CROCO sources themselves."""

from __future__ import annotations

from typing import Iterable

from jobcomp.context import BuildContext
from jobcomp.copy import CopySpec

from .base import Component


class ModelSources(Component):
    """The CROCO sources themselves plus their local includes/Makefiles."""

    name = "sources"

    SOURCE_PATTERNS = ("*.F", "*.F90", "*.h", "*.py", "Make*", "amr.in",
                       "croco_ascii.txt")

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        for pattern in self.SOURCE_PATTERNS:
            yield CopySpec.glob(ctx.paths.source, ctx.paths.scratch, pattern)
