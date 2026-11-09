"""Common protocol and scaffolding for source-to-source passes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional, Union

from SCons.Environment import Environment
from SCons.Node import NodeList
from SCons.Node.FS import File
from SCons.Scanner import Scanner

from ..croco_env import CrocoEnv
from ..helpers import Helpers


class CompilePassBase(ABC):
    """One step in the Fortran source pipeline."""

    @abstractmethod
    def apply(
        self, input_files: Union[List, File, str], **kwargs
    ) -> NodeList:
        """Apply the pass and return the new pipeline files."""
        raise NotImplementedError


class CommandPass(CompilePassBase):
    """A pass that produces one staged file per input source."""

    stage = ""
    action = ""
    skip_f90 = False

    def __init__(
        self,
        croco_env: CrocoEnv,
        skip_f90_files: Optional[bool] = None,
    ) -> None:
        self.croco_env = croco_env
        self.env: Environment = croco_env.env
        self.helpers = Helpers(croco_env)
        self.dummy_scanner = Scanner(function=lambda node, env, path: [])
        self.skip_f90_files = (
            self.skip_f90 if skip_f90_files is None else skip_f90_files
        )

    def target(self, source: File) -> File:
        return self.helpers.file_add_stage(source, self.stage)

    def sources(self, source: File) -> List[File]:
        return [source]

    def apply(
        self, input_files: Union[List, File, str], **kwargs
    ) -> NodeList:
        files = self.helpers.convert_to_nodelist_file(input_files)
        output = NodeList()
        for source in files:
            if self.skip_f90_files and self.helpers.is_f90(source):
                output += [source]
                continue
            output += self.env.Command(
                target=self.target(source),
                source=self.sources(source),
                action=self.action,
                source_scanner=self.dummy_scanner,
                **kwargs,
            )
        return output
