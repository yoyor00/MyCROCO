"""Build the default CROCO executable through the compile-pass pipeline."""

from __future__ import annotations

from typing import List, Optional

from .compile_passes.compile_pass_base import CompilePassBase
from .compile_passes.croco_change_loops import CrocoChangeLoops
from .compile_passes.croco_mpc import CrocoMPC
from .compile_passes.croco_preprocessor import CrocoPreprocessor
from .croco_build import CrocoBuild
from .croco_tools import CrocoTools


class CrocoBin:
    """Class to build the default CROCO binary."""

    def __init__(
        self,
        croco_build: CrocoBuild,
        compile_passes_pre: Optional[List[CompilePassBase]] = None,
        compile_passes_default: Optional[List[CompilePassBase]] = None,
        compile_passes_post: Optional[List[CompilePassBase]] = None,
        croco_tools: Optional[CrocoTools] = None,
    ):
        self.croco_build = croco_build
        self.env = croco_build.env

        # Tool builders must exist before the generated check_* sources are
        # compiled.  Reuse the caller's instance when one is provided.
        self.croco_tools = croco_tools or CrocoTools(self.croco_build)

        # Source lists come from the jobcomp manifest, not from a second
        # parse of the Makefile.
        sources = self.croco_build.croco_env.config.sources
        pipe = [self.env.File(name) for name in sources.fixed]
        pipe += [self.env.File(name) for name in sources.free]

        if compile_passes_default is None:
            compile_passes_default = [
                CrocoPreprocessor(self.croco_build, skip_f90_files=False),
                CrocoChangeLoops(self.croco_build, skip_f90_files=True),
                CrocoMPC(
                    self.croco_build,
                    self.croco_build.dep_mpc,
                    skip_f90_files=True,
                ),
            ]

        for compile_pass in compile_passes_pre or []:
            pipe = compile_pass.apply(pipe)
        for compile_pass in compile_passes_default:
            pipe = compile_pass.apply(pipe)
        for compile_pass in compile_passes_post or []:
            pipe = compile_pass.apply(pipe)

        objects = self.croco_build.compile_fort(pipe)
        self.program = self.croco_build.link(objects, output_file="croco")
