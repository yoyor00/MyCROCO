
from __future__ import annotations
import os
from typing import Union, List

from SCons.Node.FS import File
from SCons.Node import Node, NodeList
from SCons.Environment import Environment

from .croco_env import CrocoEnv


class Helpers:
    """
    Some little helper routines
    """

    def __init__(self, croco_env: CrocoEnv):
        self.croco_env: CrocoEnv = croco_env
        self.env: Environment = croco_env.env

    def _derive_new_file(self, original_file: File, name: str):
        """
        Add a new file with given name.
        Also inherit original_file attribute if it exists or use current one.
        """
        file: File = self.env.File(name)
        if hasattr(original_file.attributes, "original_file"):
            file.attributes.original_file = (
                original_file.attributes.original_file
            )
        else:
            file.attributes.original_file = original_file

        return file

    def file_add_stage(self, file: File, stage_name: str) -> File:
        """
        Add a stage to an scons `File`.

        E.g.
        `foo.f90` => `foo.stage_new.f90`
        """

        assert isinstance(file, File)
        base, ext = os.path.splitext(file.name)
        assert ext[0] == "."
        ext = ext[1:]

        return self._derive_new_file(file, f"{base}.{stage_name}.{ext}")

    def file_change_ext(self, file: File, new_ext: str) -> File:
        """
        Swap the extension of an scons `File`.

        E.g.
        `foo.f90` => `foo.o`
        """

        assert isinstance(file, File)

        base, ext = os.path.splitext(file.name)
        assert ext[0] == "."

        return self._derive_new_file(file, f"{base}.{new_ext}")

    def file_remove_stages(self, file: File) -> File:
        """
        Remove stages from an scons `File`.

        This is required to obtain the original object file output name
        to create dependencies to it.

        E.g.
        `foo.stage1.stage2.stage3.f90` => `foo.f90`
        """

        assert isinstance(file, File)

        s = file.name.split(".")
        assert len(s) >= 2

        return self._derive_new_file(file, f"{s[0]}.{s[-1]}")

    def convert_to_nodelist_file(
        self, input_files: Union[List, File, str]
    ) -> NodeList:
        """
        Convert input files to be a `NodeList` of scons `File`.
        """

        if type(input_files) is str:
            ret_files = NodeList([self.env.File(input_files)])
        elif isinstance(input_files, Node):
            ret_files = NodeList([input_files])
        elif isinstance(input_files, NodeList):
            ret_files = input_files
        elif isinstance(input_files, list):
            ret_files = NodeList(input_files)
        else:
            raise TypeError(
                "Don't know how to handle type "
                f"'{type(input_files)}: {input_files}'"
            )

        for idx, file in enumerate(ret_files):
            file: File

            if isinstance(file, Node):
                continue

            if isinstance(file, str):
                ret_files[idx] = self.env.File(file)
                continue

            raise Exception(f"Don't know how to handle type '{type(file)}'")

        return ret_files

    def is_f90(self, file: File) -> File:
        """
        return True if the file has a .f90 extension, False otherwise.
        """

        assert isinstance(file, File)

        base, ext = os.path.splitext(file.name)
        assert ext[0] == "."

        return ext[1:].lower() == "f90"
