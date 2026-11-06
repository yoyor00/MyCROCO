"""``BuildContext``: everything a build step may need."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Optional

from .copy import Copier, CopyMode, CopySpec
from .features import FeatureDetector
from .shell import Shell
from .toolchain import Paths, Toolchain


@dataclass
class BuildContext:
    """One build: where, with which tools, and how."""

    paths: Paths
    tools: Toolchain
    detector: FeatureDetector
    shell: Shell
    copy_mode: CopyMode = CopyMode.NEWER
    dry_run: bool = False
    keep_scratch: bool = False
    delete_stale: bool = False
    xios_dir: Path = Path(".")
    prism_dir: Path = Path(".")
    copier: Optional[Copier] = None
    features: Dict[str, bool] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.copier is None:
            self.copier = Copier(
                mode=self.copy_mode,
                dry_run=self.dry_run,
                delete_stale=self.delete_stale,
            )

    def feature(self, name: str) -> bool:
        """Was feature ``name`` detected?  (Cached by the components.)"""
        return bool(self.features.get(name, False))

    # -- convenience accessors used by the components ----------------------
    @property
    def scratch(self) -> Path:
        return self.paths.scratch

    @property
    def run_dir(self) -> Path:
        return self.paths.run_dir

    def copy_many(self, specs: Iterable[CopySpec]) -> None:
        for spec in specs:
            self.copier.apply(spec)
