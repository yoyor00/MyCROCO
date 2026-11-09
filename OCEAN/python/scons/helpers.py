"""Small helpers for manipulating SCons ``File`` nodes and names."""

from __future__ import annotations

import os
from typing import List, Union

from SCons.Environment import Environment
from SCons.Node import Node, NodeList
from SCons.Node.FS import File

from .croco_env import CrocoEnv


class Helpers:
    """Name and node helpers shared by the SCons build classes."""

    def __init__(self, croco_env: CrocoEnv):
        self.croco_env = croco_env
        self.env: Environment = croco_env.env

    def _derive_new_file(self, original_file: File, name: str) -> File:
        """Create ``name`` and propagate the original source attribute."""
        file = self.env.File(name)
        file.attributes.original_file = getattr(
            original_file.attributes, "original_file", original_file
        )
        return file

    @staticmethod
    def _split_name(file: File) -> tuple[str, str]:
        base, extension = os.path.splitext(file.name)
        if not extension:
            raise ValueError(f"file has no extension: {file.name}")
        return base, extension[1:]

    def file_add_stage(self, file: File, stage_name: str) -> File:
        """Add a stage suffix: ``foo.f90`` -> ``foo.stage.f90``."""
        base, extension = self._split_name(file)
        return self._derive_new_file(
            file, f"{base}.{stage_name}.{extension}"
        )

    def file_change_ext(self, file: File, new_ext: str) -> File:
        """Change the extension: ``foo.f90`` -> ``foo.o``."""
        base, _ = self._split_name(file)
        return self._derive_new_file(file, f"{base}.{new_ext}")

    def file_remove_stages(self, file: File) -> File:
        """Drop intermediate stages: ``foo.pre.cl.f90`` -> ``foo.f90``."""
        base, extension = self._split_name(file)
        stem = base.split(".", 1)[0]
        return self._derive_new_file(file, f"{stem}.{extension}")

    def convert_to_nodelist_file(
        self, input_files: Union[List, File, str, Node]
    ) -> NodeList:
        """Convert strings, nodes or lists into a ``NodeList`` of files."""
        if isinstance(input_files, (str, File)):
            items = [input_files]
        elif isinstance(input_files, NodeList):
            items = list(input_files)
        elif isinstance(input_files, Node):
            items = [input_files]
        elif isinstance(input_files, list):
            items = input_files
        else:
            raise TypeError(
                f"don't know how to handle type {type(input_files)!r}"
            )

        nodes = []
        for item in items:
            if isinstance(item, str):
                nodes.append(self.env.File(item))
            elif isinstance(item, Node):
                nodes.append(item)
            else:
                raise TypeError(
                    f"don't know how to handle type {type(item)!r}"
                )
        return NodeList(nodes)

    def is_f90(self, file: File) -> bool:
        """Return True for free-form Fortran files."""
        _, extension = self._split_name(file)
        return extension.lower() == "f90"
