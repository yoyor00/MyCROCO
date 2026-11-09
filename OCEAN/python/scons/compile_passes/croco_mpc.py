"""Run the CROCO ``mpc`` source pass."""

from __future__ import annotations

from typing import List

from SCons.Node.FS import File

from ..croco_env import CrocoEnv
from .compile_pass_base import CommandPass


class CrocoMPC(CommandPass):
    """Apply ``mpc`` to fixed-form sources."""

    stage = "mpc"
    action = './mpc < "$SOURCE" > "$TARGET"'
    skip_f90 = True

    def __init__(
        self,
        croco_env: CrocoEnv,
        dep_mpc: File,
        skip_f90_files: bool = True,
    ):
        super().__init__(croco_env, skip_f90_files=skip_f90_files)
        self.dep_mpc = dep_mpc

    def sources(self, source: File) -> List[File]:
        return [source, self.dep_mpc]
