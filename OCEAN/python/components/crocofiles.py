"""The ``CROCOFILES`` header tree consumed by the AGRIF ``conv`` tool."""

from __future__ import annotations

import shutil

from jobcomp.context import BuildContext
from jobcomp.errors import LOG
from jobcomp.preprocessor import Preprocessor


class Crocofiles:
    """Preprocess the AGRIF headers into ``Compile/CROCOFILES``."""

    def __init__(self, ctx: BuildContext) -> None:
        self.ctx = ctx

    def generate(self) -> None:
        """Build the ``CROCOFILES`` header tree."""
        scratch = self.ctx.paths.scratch
        crocofiles = scratch / "CROCOFILES"
        headers = sorted(scratch.glob("*.h")) + sorted(scratch.glob("*.h90"))
        if self.ctx.dry_run:
            LOG.info(
                "would create %s (+%d preprocessed headers)",
                crocofiles,
                len(headers),
            )
            return
        (crocofiles / "AGRIF_MODELFILES").mkdir(parents=True, exist_ok=True)
        (crocofiles / "AGRIF_INC").mkdir(parents=True, exist_ok=True)

        amr_cpp = Preprocessor(self.ctx.shell, self.ctx.tools.cpp)
        amr = amr_cpp.run(amr_cpp.argv("amr.in"), cwd=scratch)
        (crocofiles / "amr.scrum").write_text(Preprocessor.filter(amr))

        conv = scratch / "AGRIF" / "conv"
        if conv.exists():
            shutil.move(str(conv), str(crocofiles / conv.name))

        cppdefs = scratch / "cppdefs.h"
        prefix = ""
        if cppdefs.exists():
            prefix = cppdefs.read_text(errors="replace")
        # The bash script preprocesses the headers with plain ``cpp -P``,
        # independently of the compiler-specific CPP1 setting.
        header_cpp = Preprocessor(self.ctx.shell, "cpp")
        for path in headers:
            LOG.info("%s", path.name)
            data = prefix + path.read_text(errors="replace")
            output = header_cpp.run(
                header_cpp.argv("-P"), text_input=data, cwd=scratch
            )
            (crocofiles / path.name).write_text(Preprocessor.filter(output))

        amr_header = crocofiles / "private_scratch_AMR.h"
        if amr_header.exists():
            amr_header.replace(crocofiles / "private_scratch.h")
