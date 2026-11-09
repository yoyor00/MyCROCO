
from __future__ import annotations
from typing import Union, List
from SCons.Environment import Environment
from SCons.Node.FS import File
from SCons.Node import NodeList
from SCons.Scanner import Scanner

from ..croco_env import CrocoEnv
from ..helpers import Helpers


class CrocoSed:
    """
    This class provides the basic build calls to hook in `sed`
    """

    def __init__(self, croco_env: CrocoEnv):
        self.croco_env: CrocoEnv = croco_env
        self.env: Environment = croco_env.env

        self.helpers: Helpers = Helpers(croco_env)

        self.dummy_scanner = Scanner(
            function=lambda node, env, path: []
        )

        self.enabled = True
        self.sed_bin: str = "sed"
        self.sed_options: str = ""
        self.files_included: list[str] = []
        self.files_excluded: list[str] = []

    def apply(self, input_files: Union[List, File, str], **kwargs) -> NodeList:
        """
        Apply `sed` on a file-by-file basis
        """

        if not self.enabled:
            return input_files

        input_files = self.helpers.convert_to_nodelist_file(input_files)

        output_deps: list[File] = NodeList()

        for input_file in input_files:
            input_file: File

            if 1:
                # Only run `sed` on files as specified
                input_file_base = input_file.name.split(".")[0]

                if len(self.files_included) > 0:
                    # By default, sed runs on all files unless an include
                    # list is provided.
                    if input_file_base not in self.files_included:
                        output_deps += [input_file]
                        continue

                if input_file_base in self.files_excluded:
                    output_deps += [input_file]
                    continue

            output_file = self.helpers.file_add_stage(input_file, "sed")

            output_deps += self.env.Command(
                target=output_file,
                source=input_file,
                action=(
                    f'{self.sed_bin} {self.sed_options} "'
                    f'{input_file}" > "{output_file}"'
                ),
                source_scanner=self.dummy_scanner,  # Disable scanner
                **kwargs,
            )

        return output_deps
