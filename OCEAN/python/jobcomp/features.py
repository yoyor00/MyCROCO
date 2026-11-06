"""Feature detection: one preprocessor run instead of one ``grep`` per key.

The bash script preprocessed ``testkeys.F`` once per feature
(``$CPP1 testkeys.F | grep -i -q <key>``); the answer only depends on the
``cppdefs.h`` of the run directory, so it is computed once and cached.
"""

from __future__ import annotations

import shlex
from pathlib import Path
from typing import Optional

from .errors import LOG
from .shell import Shell


class FeatureDetector:
    """Answers ``has("<key>isdefined")`` from a single preprocessor run."""

    def __init__(
        self,
        shell: Shell,
        cpp: str,
        scratch: Path,
        source: Path,
        run_dir: Optional[Path] = None,
    ) -> None:
        self.shell = shell
        self.cpp = cpp
        self.scratch = Path(scratch)
        self.source = Path(source)
        self.run_dir = run_dir if run_dir is not None else Path.cwd()
        self._output: Optional[str] = None

    def has(self, key: str) -> bool:
        if not key:
            return True
        return key.lower() in self._preprocess().lower()

    def _preprocess(self) -> str:
        """Preprocess ``cppdefs.h`` + ``testkeys.F``, once."""
        if self._output is not None:
            return self._output
        self._output = ""
        testkeys = self._find("testkeys.F")
        if testkeys is None or not self.cpp:
            LOG.debug("no testkeys.F or no cpp: cannot detect features")
            return self._output
        # The include is resolved by hand so that the same code works before
        # (dry run) and after the copy phase: the scratch copy wins, then the
        # run directory, then the sources.
        header = self._find("cppdefs.h")
        text = "\n".join(
            line
            for line in testkeys.read_text(errors="replace").splitlines()
            if line.strip().replace(" ", "") != '#include"cppdefs.h"'
        )
        if header is not None:
            text = header.read_text(errors="replace") + "\n" + text
        argv = shlex.split(self.cpp) + ["-P"]
        self._output = self.shell.capture(
            argv, cwd=testkeys.parent, text_input=text
        )
        LOG.debug("feature detection: cpp on %s + %s", header, testkeys)
        return self._output

    def _find(self, name: str) -> Optional[Path]:
        """First ``name`` in the scratch tree, the run dir or the sources."""
        for directory in (self.scratch, self.run_dir, self.source):
            candidate = directory / name
            if candidate.exists():
                return candidate
        return None
