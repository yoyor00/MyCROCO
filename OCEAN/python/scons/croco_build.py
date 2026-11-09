"""Basic SCons build operations: preprocess, compile, link and copy stages."""

from __future__ import annotations

import os
from typing import List, Optional, Union

from SCons.Environment import Environment
from SCons.Node import NodeList
from SCons.Node.FS import File
from SCons.Scanner import Scanner

from .compile_passes.croco_mpc import CrocoMPC
from .compile_passes.croco_preprocessor import CrocoPreprocessor
from .croco_env import CrocoEnv
from .deps.deps_fcompiler import FCompiler
from .helpers import Helpers


class CrocoBuild:
    """Provide the elementary build calls used by the compile pipeline."""

    def __init__(self, croco_env: CrocoEnv):
        self.croco_env = croco_env
        self.env: Environment = croco_env.env
        self.helpers = Helpers(croco_env)
        self.dummy_scanner = Scanner(function=lambda node, env, path: [])
        self.preprocessor = CrocoPreprocessor(croco_env)

        # ``mpc`` is required by the preprocessing passes.
        self.dep_mpc = self.compile_link(
            input_files=self.preprocessor.apply("mpc.F"),
            output_file="mpc",
            LIBS=[],
        )
        self.mpc = CrocoMPC(croco_env, self.dep_mpc)

    # -- Fortran -----------------------------------------------------------
    def compile_fort(
        self,
        input_files: Union[List, File, str],
        input_additional_deps: Optional[NodeList] = None,
        rename_to_orig_filename: bool = True,
        fortran_input: Optional[str] = None,
        **kwargs,
    ) -> NodeList:
        """Compile one or more Fortran files."""
        inputs = self.helpers.convert_to_nodelist_file(input_files)
        additional = (
            self.helpers.convert_to_nodelist_file(input_additional_deps)
            if input_additional_deps is not None
            else None
        )

        output_deps = NodeList()
        for input_file in inputs:
            output_file = self.helpers.file_change_ext(input_file, "o")
            if rename_to_orig_filename:
                output_file = self.helpers.file_remove_stages(output_file)

            sources = NodeList([input_file])
            if additional is not None:
                sources += additional

            version = self._fortran_version(input_file, fortran_input)
            action = self._compile_action(version)
            output_files = [output_file]
            self.env.Command(
                target=output_files, source=sources, action=action, **kwargs
            )
            self._declare_module_side_effects(input_file, output_files)
            output_deps += NodeList(output_files)

        return output_deps

    @staticmethod
    def _fortran_version(source: File, requested: Optional[str]) -> str:
        if requested is not None:
            version = requested
        elif source.name.endswith(".f"):
            version = "f"
        elif source.name.endswith(".f90"):
            version = "f90"
        else:
            raise ValueError(f"Unknown file ending for '{source}'")
        if version not in ("f", "f90"):
            raise ValueError(f"Unknown Fortran input version: {version}")
        return version

    @staticmethod
    def _compile_action(version: str) -> str:
        if version == "f":
            return "$CFT $FFLAGS $_CPPINCFLAGS -c $SOURCE -o $TARGET"
        return "$CFT $FFLAGS90 $_CPPINCFLAGS -c $SOURCE -o $TARGET"

    def _declare_module_side_effects(
        self, source: File, output_files: List[File]
    ) -> None:
        """Declare the .mod files produced by ``source`` as side effects."""
        original = getattr(source.attributes, "original_file", None)
        if original is None and os.path.exists(source.name):
            original = source
        if original is None:
            return

        modules = FCompiler(original).defined_modules()
        suffix = self.env.subst("$FORTRANMODSUFFIX")
        self.env.SideEffect(
            [f"{module}{suffix}" for module in modules], output_files
        )

    # -- linking -----------------------------------------------------------
    def link(
        self,
        input_objects: Union[List, File, str],
        output_file: Union[File, str],
        **kwargs,
    ) -> File:
        """Link the object files into one executable."""
        inputs = self.helpers.convert_to_nodelist_file(input_objects)
        if not isinstance(output_file, File):
            output_file = self.env.File(output_file)
        # $FFLAGS is needed at link time too for options such as -acc=gpu.
        action = (
            "$LDR $FFLAGS $SOURCES -o $TARGET "
            "$_LIBDIRFLAGS $_LIBFLAGS $LDFLAGS"
        )
        out_deps = self.env.Command(
            target=output_file, source=inputs, action=action, **kwargs
        )
        return out_deps[0]

    # -- staging helpers ---------------------------------------------------
    def copy_stage(
        self, input_files: Union[List, File, str], stage_name: str
    ) -> NodeList:
        """Copy a file to a new stage name (``foo.F`` -> ``foo.tools.F``)."""
        inputs = self.helpers.convert_to_nodelist_file(input_files)
        output_deps = NodeList()
        for input_file in inputs:
            output_file = self.helpers.file_add_stage(input_file, stage_name)
            output_deps += self.env.Command(
                target=output_file,
                source=input_file,
                action='cp "$SOURCE" "$TARGET"',
                source_scanner=self.dummy_scanner,
            )
        return output_deps

    # -- common pipelines --------------------------------------------------
    def preproc_mpc(self, input_files: Union[List, File, str]) -> NodeList:
        return self.mpc.apply(self.preprocessor.apply(input_files))

    def compile_link(
        self,
        input_files: Union[List, File, str],
        output_file: str,
        **kwargs,
    ) -> File:
        return self.link(
            self.compile_fort(input_files), output_file=output_file, **kwargs
        )
