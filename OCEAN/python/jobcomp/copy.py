"""Copying files and directory trees into the scratch directory.

The bash script copied everything unconditionally (``cp`` / ``cp -r``); here a
copy is described by a :class:`CopySpec` and applied by a :class:`Copier`,
which owns the "copy always / only if newer" rule and the counters.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from .errors import LOG, BuildError


class CopyMode(str, Enum):
    """When should a copy actually happen?"""

    #: Only when the destination is missing, older or of another size.
    NEWER = "newer"
    #: Unconditional overwrite (the bash script's behaviour).
    ALWAYS = "always"


@dataclass
class CopyStats:
    """Accounting for a copy phase."""

    copied: int = 0
    skipped: int = 0
    removed: int = 0
    bytes_copied: int = 0

    def __str__(self) -> str:
        return (
            f"{self.copied} copied ({self.bytes_copied} bytes), "
            f"{self.skipped} unchanged, {self.removed} removed"
        )


@dataclass(frozen=True)
class CopySpec:
    """Declarative description of "what goes where, under which rule".

    ``src`` is a directory, a file, or the base directory of ``glob``.
    ``tree`` copies ``src`` as a whole into ``dst/<name>`` (``cp -r src dst``).
    ``contents`` copies the entries of ``src`` into ``dst``
    (``cp -r src/* dst``), recursing into sub-directories.
    """

    src: Path
    dst: Path
    glob: Optional[str] = None
    tree: bool = False
    contents: bool = False
    overwrite: bool = True
    optional: bool = True
    mode: Optional[CopyMode] = None  # None: the mode of the Copier

    def replace(self, *, src=None, dst=None) -> "CopySpec":
        """The same rule applied to another source and/or destination."""
        return CopySpec(
            src=self.src if src is None else src,
            dst=self.dst if dst is None else dst,
            glob=self.glob,
            tree=self.tree,
            contents=self.contents,
            overwrite=self.overwrite,
            optional=self.optional,
            mode=self.mode,
        )


class Copier:
    """Applies :class:`CopySpec` objects with one shared mode and counter."""

    def __init__(
        self,
        mode: CopyMode = CopyMode.NEWER,
        dry_run: bool = False,
        delete_stale: bool = False,
        log: logging.Logger = LOG,
    ) -> None:
        self.mode = mode
        self.dry_run = dry_run
        self.delete_stale = delete_stale
        self.log = log
        self.stats = CopyStats()

    # -- single files ------------------------------------------------------
    def file(
        self,
        src,
        dst,
        *,
        overwrite: bool = True,
        mode: Optional[CopyMode] = None,
    ) -> bool:
        """Copy ``src`` to ``dst``; returns True when something was copied.

        ``dst`` may be a directory, in which case the file keeps its name.
        ``overwrite=False`` keeps an existing destination untouched (this is
        how "do not overwrite the run directory namelists" is expressed).
        ``mode`` overrides the mode of the copier for this one copy.
        """
        mode = self.mode if mode is None else mode
        src = Path(src)
        dst = Path(dst)
        if not src.exists():
            if self.dry_run:
                self.log.info(
                    "would copy %s -> %s (missing in dry run)", src, dst
                )
                return True
            raise FileNotFoundError(str(src))
        if dst.is_dir():
            dst = dst / src.name
        if src.resolve() == dst.resolve():
            self.stats.skipped += 1
            return False

        reason = "new"
        if dst.exists():
            if not overwrite:
                self.stats.skipped += 1
                self.log.debug("keep %s (already present)", dst)
                return False
            if mode is CopyMode.NEWER:
                src_stat = src.stat()
                dst_stat = dst.stat()
                if (
                    src_stat.st_size == dst_stat.st_size
                    and src_stat.st_mtime_ns <= dst_stat.st_mtime_ns
                ):
                    self.stats.skipped += 1
                    self.log.debug("keep %s (up to date)", dst)
                    return False
                reason = "mtime"
        size = src.stat().st_size
        if self.dry_run:
            self.log.info("would copy %s -> %s (%s)", src, dst, reason)
            self.stats.copied += 1
            return True
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)  # copy2 keeps mtime: make stays incremental
        self.stats.copied += 1
        self.stats.bytes_copied += size
        self.log.debug("copied %s -> %s (%s)", src, dst, reason)
        return True

    # -- directory trees ---------------------------------------------------
    def tree(
        self, src, dst, *, mode: Optional[CopyMode] = None
    ) -> CopyStats:
        """Mirror the directory ``src`` into ``dst``, honouring the mode."""
        src = Path(src)
        dst = Path(dst)
        if not src.is_dir():
            raise NotADirectoryError(str(src))
        seen = set()
        for path in sorted(src.rglob("*")):
            rel = path.relative_to(src)
            target = dst / rel
            if path.is_dir():
                if not self.dry_run:
                    target.mkdir(parents=True, exist_ok=True)
                continue
            seen.add(rel)
            self.file(path, target, mode=mode)
        if self.delete_stale and dst.is_dir():
            for path in sorted(dst.rglob("*")):
                rel = path.relative_to(dst)
                if not path.is_file() or rel in seen:
                    continue
                if self.dry_run:
                    self.log.info("would remove %s", path)
                else:
                    path.unlink()
                self.stats.removed += 1
        return self.stats

    # -- specs -------------------------------------------------------------
    def apply(self, spec: CopySpec) -> None:
        """Apply one :class:`CopySpec` (glob, tree, contents or plain file)."""
        src = spec.src
        mode = self.mode if spec.mode is None else spec.mode
        if spec.glob is not None:
            if not src.is_dir():
                if spec.optional:
                    return
                raise BuildError(f"missing source directory: {src}")
            matches = sorted(src.glob(spec.glob))
            if not matches and not spec.optional:
                raise BuildError(f"no {spec.glob} in {src}")
            for match in matches:
                self._entry(match, spec)
        elif spec.tree or spec.contents:
            if not src.exists():
                if spec.optional:
                    return
                raise BuildError(f"missing source directory: {src}")
            if spec.tree:
                self.tree(spec.src, spec.dst / src.name, mode=mode)
            else:
                for entry in sorted(src.iterdir()):
                    self._entry(entry, spec)
        else:
            if not src.exists():
                if spec.optional:
                    return
                raise BuildError(f"missing source file: {src}")
            self.file(
                src, spec.dst, overwrite=spec.overwrite, mode=mode
            )

    def _entry(self, path: Path, spec: CopySpec) -> None:
        if path.is_dir():
            # ``cp -r src/<dir> dst/<dir>``: the child becomes the new source
            self.tree(path, spec.dst / path.name, mode=spec.mode)
        else:
            self.file(
                path,
                spec.dst,
                overwrite=spec.overwrite,
                mode=spec.mode,
            )

    def finish(self) -> None:
        self.log.info("Copy phase: %s", self.stats)
