"""Running external commands, with the driver's logging and error style."""

from __future__ import annotations

import logging
import shlex
import subprocess
from pathlib import Path
from typing import Optional, Sequence

from .errors import LOG, BuildError


class Shell:
    """The one place where the driver talks to the outside world.

    ``dry_run`` turns every command into a log line, so the same code produces
    the build plan and the build itself.
    """

    def __init__(
        self, dry_run: bool = False, log: logging.Logger = LOG
    ) -> None:
        self.dry_run = dry_run
        self.log = log

    def run(
        self,
        argv: Sequence[str],
        *,
        cwd: Optional[Path] = None,
        text_input: Optional[str] = None,
        capture: bool = False,
        check: bool = True,
    ) -> subprocess.CompletedProcess:
        """Run a command, logging it.  Raises :class:`BuildError` on failure.

        Without ``capture`` the child inherits our streams, so ``make`` output
        is visible exactly like it was with the bash script.
        """
        pretty = self._pretty(argv)
        if cwd is not None:
            self.log.debug("run (in %s): %s", cwd, pretty)
        else:
            self.log.debug("run: %s", pretty)
        if self.dry_run:
            self.log.info("would run: %s", pretty)
            return subprocess.CompletedProcess(argv, 0, "", "")
        proc = self._call(
            argv, cwd=cwd, text_input=text_input, capture=capture
        )
        if check and proc.returncode != 0:
            raise BuildError(
                f"command failed (exit {proc.returncode}): {pretty}",
                returncode=proc.returncode,
            )
        return proc

    def probe(
        self,
        argv: Sequence[str],
        *,
        cwd: Optional[Path] = None,
        text_input: Optional[str] = None,
    ) -> subprocess.CompletedProcess:
        """Run a read-only query, even in a dry run.

        Feature detection, compiler probing and ``make`` variable lookups must
        answer the same way whether the plan is executed or only printed.
        """
        self.log.debug("probe: %s", self._pretty(argv))
        try:
            return self._call(
                argv, cwd=cwd, text_input=text_input, capture=True
            )
        except OSError:  # command not installed: empty answer, never fatal
            return subprocess.CompletedProcess(argv, 127, "", "")

    def capture(
        self,
        argv: Sequence[str],
        *,
        cwd: Optional[Path] = None,
        text_input: Optional[str] = None,
    ) -> str:
        """Standard output of a read-only query (``""`` if it cannot run)."""
        return self.probe(argv, cwd=cwd, text_input=text_input).stdout or ""

    # -- plumbing ----------------------------------------------------------
    @staticmethod
    def _pretty(argv: Sequence[str]) -> str:
        return " ".join(shlex.quote(str(part)) for part in argv)

    @staticmethod
    def _call(
        argv: Sequence[str],
        *,
        cwd: Optional[Path],
        text_input: Optional[str],
        capture: bool,
    ) -> subprocess.CompletedProcess:
        return subprocess.run(
            list(argv),
            cwd=str(cwd) if cwd is not None else None,
            input=text_input,
            text=True,
            capture_output=capture,
        )
