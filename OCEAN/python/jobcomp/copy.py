"""Copying files and directory trees into the scratch directory.

The bash script copied everything unconditionally (``cp`` / ``cp -r``); here a
copy is described by a :class:`CopySpec` and applied by a :class:`Copier`,
which owns the "copy always / only if newer" rule and the counters.

The four kinds of copy are explicit (:class:`CopyKind`), so a call site reads
``CopySpec.tree(...)`` or ``CopySpec.contents(...)`` instead of a combination
of boolean flags.
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


class CopyKind(str, Enum):
    """What a :class:`CopySpec` describes."""

    FILE = "file"
    GLOB = "glob"
    TREE = "tree"
    CONTENTS = "contents"


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

    ``src`` is a directory, a file, or the base directory of ``pattern``.
    ``TREE`` copies ``src`` as a whole into ``dst/<name>`` (``cp -r src dst``).
    ``CONTENTS`` copies the entries of ``src`` into ``dst``
    (``cp -r src/* dst``), recursing into sub-directories.
    """

    src: Path
    dst: Path
    kind: CopyKind = CopyKind.FILE
    pattern: Optional[str] = None
    overwrite: bool = True
    optional: bool = True
    mode: Optional[CopyMode] = None  # None: the mode of the Copier

    @classmethod
    def file(cls, src, dst, **kwargs) -> "CopySpec":
        """A single-file copy (``cp src dst``)."""
        return cls(Path(src), Path(dst), CopyKind.FILE, **kwargs)

    @classmethod
    def glob(cls, src, dst, pattern: str, **kwargs) -> "CopySpec":
        """Files matching ``pattern`` in ``src`` (``cp src/<glob> dst``)."""
        return cls(
            Path(src), Path(dst), CopyKind.GLOB, pattern=pattern, **kwargs
        )

    @classmethod
    def tree(cls, src, dst, **kwargs) -> "CopySpec":
        """A whole directory tree (``cp -r src dst/<name>``)."""
        return cls(Path(src), Path(dst), CopyKind.TREE, **kwargs)

    @classmethod
    def contents(cls, src, dst, **kwargs) -> "CopySpec":
        """The entries of one directory (``cp -r src/* dst``)."""
        return cls(Path(src), Path(dst), CopyKind.CONTENTS, **kwargs)


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
        """Apply one :class:`CopySpec`."""
        if spec.kind is CopyKind.GLOB:
            self._apply_glob(spec)
        elif spec.kind is CopyKind.TREE:
            self._apply_tree(spec)
        elif spec.kind is CopyKind.CONTENTS:
            self._apply_contents(spec)
        else:
            self._apply_file(spec)

    def _apply_glob(self, spec: CopySpec) -> None:
        if not spec.src.is_dir():
            self._missing(spec, f"missing source directory: {spec.src}")
            return
        matches = sorted(spec.src.glob(spec.pattern))
        if not matches:
            self._missing(spec, f"no {spec.pattern} in {spec.src}")
            return
        for match in matches:
            self._entry(match, spec, spec.dst)

    def _apply_tree(self, spec: CopySpec) -> None:
        if not spec.src.exists():
            self._missing(spec, f"missing source directory: {spec.src}")
            return
        self.tree(spec.src, spec.dst / spec.src.name, mode=spec.mode)

    def _apply_contents(self, spec: CopySpec) -> None:
        if not spec.src.exists():
            self._missing(spec, f"missing source directory: {spec.src}")
            return
        for entry in sorted(spec.src.iterdir()):
            self._entry(entry, spec, spec.dst)

    def _apply_file(self, spec: CopySpec) -> None:
        if not spec.src.exists():
            self._missing(spec, f"missing source file: {spec.src}")
            return
        self.file(
            spec.src,
            spec.dst,
            overwrite=spec.overwrite,
            mode=spec.mode,
        )

    def _missing(self, spec: CopySpec, message: str) -> None:
        if not spec.optional:
            raise BuildError(message)

    def _entry(self, path: Path, spec: CopySpec, dst: Path) -> None:
        if path.is_dir():
            # ``cp -r src/<dir> dst/<dir>``: the child becomes the new source
            self.tree(path, dst / path.name, mode=spec.mode)
        else:
            self.file(
                path,
                dst,
                overwrite=spec.overwrite,
                mode=spec.mode,
            )

    def finish(self) -> None:
        self.log.info("Copy phase: %s", self.stats)
