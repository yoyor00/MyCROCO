
from __future__ import annotations
from typing import Union, List
from SCons.Environment import Environment
from SCons.Node.FS import File
from SCons.Node import NodeList
from SCons.Scanner import Scanner

from ..croco_env import CrocoEnv
from ..helpers import Helpers


class CrocoMPC:

    def __init__(
        self,
        croco_env: CrocoEnv,
        dep_mpc: File,
        skip_f90_files: bool = True,
    ):
        self.croco_env: CrocoEnv = croco_env
        self.env: Environment = croco_env.env
        self.skip_f90_files: bool = skip_f90_files

        # Dependency to the mpc executable used by this pass.
        self.dep_mpc: File = dep_mpc

        self.helpers: Helpers = Helpers(croco_env)

        self.dummy_scanner = Scanner(
            function=lambda node, env, path: []
        )

    def apply(self, input_files: Union[List, File, str], **kwargs) -> NodeList:
        input_files = self.helpers.convert_to_nodelist_file(input_files)
        output_deps: list[File] = NodeList()

        for input_file in input_files:
            input_file: File

            if self.skip_f90_files and self.helpers.is_f90(input_file):
                # Skip F90 files
                output_deps += [input_file]
                continue

            output_file = self.helpers.file_add_stage(input_file, "mpc")

            output_deps += self.env.Command(
                target=output_file,
                source=[input_file, self.dep_mpc],
                action="""./mpc < "$SOURCE" > "$TARGET" """,
                source_scanner=self.dummy_scanner,  # Disable scanner
                **kwargs,
            )

        return output_deps
