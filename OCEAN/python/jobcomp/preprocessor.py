"""Running ``cpp`` and filtering its output, in one place.

Both feature detection and the AGRIF ``CROCOFILES`` generation preprocess
small files.  They share the shell invocation, the debug logging and the
historical ``grep -v`` filter here.
"""

from __future__ import annotations

import shlex
from typing import List, Optional, Sequence

from .errors import LOG
from .shell import Shell

#: Lines dropped from preprocessor output, like the bash ``grep -v`` filter.
CPP_OUTPUT_DROP = ("!", "#", "%", "*")


class Preprocessor:
    """The preprocessor command of a build, with its output filter."""

    def __init__(self, shell: Shell, cpp: str = "") -> None:
        self.shell = shell
        self.cpp = cpp

    def argv(self, *extra: str) -> List[str]:
        """The configured ``cpp`` command followed by ``extra`` arguments."""
        return shlex.split(self.cpp) + [str(argument) for argument in extra]

    def run(
        self,
        argv: Sequence[str],
        *,
        text_input: Optional[str] = None,
        cwd=None,
    ) -> str:
        """Run a preprocessor command; never fatal, like ``cpp ... | grep``."""
        proc = self.shell.probe(argv, cwd=cwd, text_input=text_input)
        if proc.returncode not in (0, 127):
            LOG.debug(
                "cpp returned %s: %s", proc.returncode, proc.stderr.strip()
            )
        return proc.stdout or ""

    @staticmethod
    def filter(text: str) -> str:
        """Drop the lines the bash script filters out with ``grep -v``."""
        return "".join(
            line
            for line in text.splitlines(True)
            if not any(char in line for char in CPP_OUTPUT_DROP)
        )
