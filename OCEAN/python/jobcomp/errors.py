"""The single failure type of the driver, and the shared logger."""

from __future__ import annotations

import logging

#: One logger for the whole driver; the CLI installs the handlers.
LOG = logging.getLogger("jobcomp")


class BuildError(RuntimeError):
    """A build step failed; ``returncode`` is the shell exit status to use."""

    def __init__(self, message: str, returncode: int = 1) -> None:
        super().__init__(message)
        self.returncode = returncode
