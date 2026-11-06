"""AGRIF: the library build and the ``CROCOFILES`` header tree."""

from __future__ import annotations

import shlex
import shutil
from pathlib import Path
from typing import Iterable, List, Optional

from jobcomp.context import BuildContext
from jobcomp.copy import CopyMode, CopySpec
from jobcomp.errors import LOG

from .base import CPP_OUTPUT_DROP, Component


class AgrifComponent(Component):
    """The AGRIF library and the ``CROCOFILES`` header tree."""

    name = "agrif"
    key = "agrifisdefined"
    label = "AGRIF"
    announce = True

    def copy_specs(self, ctx: BuildContext) -> Iterable[CopySpec]:
        yield CopySpec(src=ctx.paths.agrif, dst=ctx.paths.scratch, tree=True)

    def configure(self, ctx: BuildContext) -> None:
        if not self.enabled(ctx):
            return
        tools = ctx.tools
        tools.add_fflags("-IAGRIF")
        tools.ldflags = f"-LAGRIF -lagrif {tools.ldflags}".strip()

    def prepare(self, ctx: BuildContext) -> None:
        if not ctx.feature(self.name):
            return
        tools = ctx.tools
        # The AGRIF variant of Makedefs carries the extra rules.  It must win
        # whatever the copy mode says: a stale template means a broken build.
        ctx.copier.file(
            ctx.paths.scratch / "Makedefs.generic.AGRIF",
            ctx.paths.scratch / "Makedefs.generic",
            mode=CopyMode.ALWAYS,
        )
        tools.nprocs = 1  # compile the library on one process only
        self.resolve_c_compiler(ctx)
        argv = [
            tools.make,
            "-C",
            "AGRIF",
            f"FC={tools.cft}",
            f"CPP={tools.cpp}",
        ]
        cppflags = self.cppflags(ctx)
        if cppflags:
            argv.append(f"CPPFLAGS={cppflags}")
        argv += [f"FFLAGS={tools.fflags}", f"CFLAGS={tools.cflags}"]
        ctx.shell.run(argv, cwd=ctx.paths.scratch)
        if tools.os == "Darwin":
            ctx.shell.run(["ranlib", "AGRIF/libagrif.a"],
                          cwd=ctx.paths.scratch)
        self._prepare_crocofiles(ctx)

    def cppflags(self, ctx: BuildContext) -> str:
        """Extra preprocessor flags of the AGRIF library (MPI variant)."""
        if not ctx.feature("mpi"):
            return ""
        return f"-DAGRIF_MPI {ctx.tools.mpiinc}".strip()

    # -- the C compiler the AGRIF library is built with --------------------
    def resolve_c_compiler(self, ctx: BuildContext) -> None:
        """Fill ``ctx.tools.cc``/``cflags`` the way the bash script does."""
        tools = ctx.tools
        cflags = tools.cflags or self._make_variable(ctx, "CFLAGS")
        cc = tools.cc or self._make_variable(ctx, "CC")
        if cc and shutil.which(cc) and self._is_gnu_c(ctx, cc):
            LOG.info("Using the GNU C compiler. Adding -fcommon to CFLAGS")
            cflags = f"{cflags} -fcommon -std=gnu11".strip()
        tools.cc = cc
        tools.cflags = cflags

    @staticmethod
    def _make_variable(ctx: BuildContext, name: str) -> str:
        """Ask make for the value of a variable (``CC``, ``CFLAGS``, ...)."""
        snippet = f"dummy_target:\n\t@echo $({name})\n"
        scratch = ctx.paths.scratch
        cwd = scratch if scratch.is_dir() else None
        argv = [ctx.tools.make, "-f", "-", "dummy_target"]
        return ctx.shell.capture(argv, cwd=cwd, text_input=snippet).strip()

    @staticmethod
    def _is_gnu_c(ctx: BuildContext, cc: str) -> bool:
        proc = ctx.shell.probe([cc, "-v"])
        return "gcc version" in (proc.stdout + proc.stderr)

    # -- the CROCOFILES header tree ---------------------------------------
    def _prepare_crocofiles(self, ctx: BuildContext) -> None:
        """Build the ``CROCOFILES`` header tree consumed by ``conv``."""
        scratch = ctx.paths.scratch
        crocofiles = scratch / "CROCOFILES"
        headers = sorted(scratch.glob("*.h")) + sorted(scratch.glob("*.h90"))
        if ctx.dry_run:
            LOG.info("would create %s (+%d preprocessed headers)",
                     crocofiles, len(headers))
            return
        (crocofiles / "AGRIF_MODELFILES").mkdir(parents=True, exist_ok=True)
        (crocofiles / "AGRIF_INC").mkdir(parents=True, exist_ok=True)

        amr = self._cpp(
            ctx, shlex.split(ctx.tools.cpp) + ["amr.in"], cwd=scratch
        )
        (crocofiles / "amr.scrum").write_text(self._filter(amr))

        conv = scratch / "AGRIF" / "conv"
        if conv.exists():
            shutil.move(str(conv), str(crocofiles / conv.name))

        cppdefs = scratch / "cppdefs.h"
        prefix = ""
        if cppdefs.exists():
            prefix = cppdefs.read_text(errors="replace")
        for path in headers:
            LOG.info("%s", path.name)
            data = prefix + path.read_text(errors="replace")
            output = self._cpp(
                ctx, ["cpp", "-P"], text_input=data, cwd=scratch
            )
            (crocofiles / path.name).write_text(self._filter(output))

        amr_header = crocofiles / "private_scratch_AMR.h"
        if amr_header.exists():
            amr_header.replace(crocofiles / "private_scratch.h")

    @staticmethod
    def _filter(text: str) -> str:
        """Drop the lines the bash script filters out with ``grep -v``."""
        return "".join(
            line
            for line in text.splitlines(True)
            if not any(char in line for char in CPP_OUTPUT_DROP)
        )

    @staticmethod
    def _cpp(
        ctx: BuildContext,
        argv: List[str],
        *,
        text_input: Optional[str] = None,
        cwd: Optional[Path] = None,
    ) -> str:
        """Run the preprocessor; never fatal, like ``cpp ... | grep``."""
        proc = ctx.shell.probe(argv, cwd=cwd, text_input=text_input)
        if proc.returncode not in (0, 127):
            LOG.debug(
                "cpp returned %s: %s", proc.returncode, proc.stderr.strip()
            )
        return proc.stdout or ""
