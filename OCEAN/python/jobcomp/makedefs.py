"""The ``Makedefs`` file: five substitutions in ``Makedefs.generic``.

The bash script built a ``sed`` script from the flags
(``s?$(FFLAGS1)?$FFLAGS1?g`` ...) and ran it on the template; here the
template is filled in directly, which yields the same file with tidy blanks.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

from .errors import LOG, BuildError
from .toolchain import Toolchain


class Makedefs:
    """Renders ``Makedefs.generic`` with the flags of a toolchain."""

    #: Token in the template, and the ``Toolchain`` attribute filling it.
    TOKENS: Tuple[Tuple[str, str], ...] = (
        ("$(FFLAGS1)", "fflags"),
        ("$(LDFLAGS1)", "ldflags"),
        ("$(CPP1)", "cpp"),
        ("$(CFT1)", "cft"),
        ("$(CPPFLAGS1)", "cppflags"),
    )

    def __init__(self, tools: Toolchain) -> None:
        self.tools = tools

    def items(self):
        """Yield ``(token, make variable, value)`` for every substitution."""
        for token, attribute in self.TOKENS:
            yield token, token[2:-1], getattr(self.tools, attribute)

    def variables(self) -> Dict[str, str]:
        """The make variables it defines (``FFLAGS1``, ``CFT1``, ...)."""
        return {name: value for _, name, value in self.items()}

    def render(self, template: str) -> str:
        """The ``Makedefs`` text for ``template`` (``Makedefs.generic``)."""
        text = template
        for token, _, value in self.items():
            text = text.replace(token, value)
        return text

    def write(self, directory: Path, *, dry_run: bool = False) -> None:
        """Write ``Makedefs`` next to ``Makedefs.generic`` in ``directory``."""
        template = Path(directory) / "Makedefs.generic"
        target = Path(directory) / "Makedefs"
        if not template.exists():
            if dry_run:
                LOG.info("would write %s", target)
                return
            raise BuildError(f"missing {template}")
        text = self.render(template.read_text())
        if dry_run:
            LOG.info("would write %s", target)
            return
        target.write_text(text)
