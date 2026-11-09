"""Basic SCons build operations: preprocess, compile, link and copy stages."""

from __future__ import annotations

import os
from typing import List, Union

from SCons.Environment import Environment
from SCons.Node import NodeList
from SCons.Node.FS import File
from SCons.Scanner import Scanner

from .compile_passes.croco_change_loops import CrocoChangeLoops
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

        self.changeloops = CrocoChangeLoops(croco_env)
        self.mpc = CrocoMPC(croco_env, self.dep_mpc)

    # -- Fortran / C -------------------------------------------------------
    def compile_fort(
        self,
        input_files: Union[List, File, str],
        input_additional_deps: NodeList = None,
        rename_to_orig_filename: bool = True,
        fortran_input: str = None,
        **kwargs,
    ) -> NodeList:
        """Compile one or more Fortran files."""
        inputs = self.helpers.convert_to_nodelist_file(input_files)
        additional = None
        if input_additional_deps is not None:
            additional = self.helpers.convert_to_nodelist_file(
                input_additional_deps
            )

        output_deps = NodeList()
        for input_file in inputs:
            assert isinstance(input_file, File)
            output_file = self.helpers.file_change_ext(input_file, "o")
            if rename_to_orig_filename:
                output_file = self.helpers.file_remove_stages(output_file)

            sources = NodeList([input_file])
            if additional is not None:
                sources += additional

            version = fortran_input
            if version is None:
                if input_file.name.endswith(".f"):
                    version = "f"
                elif input_file.name.endswith(".f90"):
                    version = "f90"
                else:
                    raise ValueError(f"Unknown file ending for '{input_file}'")
            if version == "f":
                action = "$CFT $FFLAGS $_CPPINCFLAGS -c $SOURCE -o $TARGET"
            elif version == "f90":
                action = "$CFT $FFLAGS90 $_CPPINCFLAGS -c $SOURCE -o $TARGET"
            else:
                raise ValueError(f"Unknown Fortran input version: {version}")

            output_files = [output_file]
            self.env.Command(
                target=output_files, source=sources, action=action, **kwargs
            )

            # Make the .mod files produced by this source known as side
            # effects of the object file.  The Fortran scanner then only has
            # to depend on the .mod file names.
            original = None
            if hasattr(input_file.attributes, "original_file"):
                original = input_file.attributes.original_file
            elif os.path.exists(input_file.name):
                original = input_file
            if original is not None:
                modules = FCompiler(original).getDefModules()
                suffix = self.env.subst("$FORTRANMODSUFFIX")
                self.env.SideEffect(
                    [f"{module}{suffix}" for module in modules],
                    output_files,
                )

            output_deps += self.helpers.convert_to_nodelist_file(output_files)

        return output_deps

    def compile_c(
        self,
        input_files: Union[List, File, str],
        rename_to_orig_filename: bool = True,
        **kwargs,
    ) -> NodeList:
        """Compile one or more C files."""
        inputs = self.helpers.convert_to_nodelist_file(input_files)
        output_deps = NodeList()
        for input_file in inputs:
            assert isinstance(input_file, File)
            output_file = self.helpers.file_change_ext(input_file, "o")
            if rename_to_orig_filename:
                output_file = self.helpers.file_remove_stages(output_file)
            output_deps += self.env.Command(
                target=output_file,
                source=input_file,
                action="$CC $CFLAGS $_CPPINCFLAGS -c $SOURCE -o $TARGET",
                **kwargs,
            )
        return output_deps

    # -- linking -----------------------------------------------------------
    def link(
        self,
        input_objects: Union[List, File, str],
        output_file: File,
        **kwargs,
    ) -> File:
        """Link the object files into one executable."""
        inputs = self.helpers.convert_to_nodelist_file(input_objects)
        if not isinstance(output_file, File):
            output_file = self.env.File(output_file)
        objects = " ".join(node.name for node in inputs)
        # $FFLAGS is needed at link time too for options such as -acc=gpu.
        action = (
            f"$LDR $FFLAGS {objects} -o {output_file} "
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
            assert isinstance(input_file, File)
            output_file = self.helpers.file_add_stage(input_file, stage_name)
            output_deps += self.env.Command(
                target=output_file,
                source=input_file,
                action='cp "$SOURCE" "$TARGET"',
                source_scanner=self.dummy_scanner,
            )
        return output_deps

    def copy_ext(
        self, input_files: Union[List, File, str], new_ext: str
    ) -> NodeList:
        """Copy a file to a new extension."""
        inputs = self.helpers.convert_to_nodelist_file(input_files)
        output_deps = NodeList()
        for input_file in inputs:
            assert isinstance(input_file, File)
            output_file = self.helpers.file_change_ext(input_file, new_ext)
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

    def preproc_changeloops_mpc(
        self, input_files: Union[List, File, str]
    ) -> NodeList:
        return self.mpc.apply(
            self.changeloops.apply(self.preprocessor.apply(input_files))
        )

    def compile_link(
        self,
        input_files: Union[List, File, str],
        output_file: str,
        **kwargs,
    ):
        return self.link(
            self.compile_fort(input_files), output_file=output_file, **kwargs
        )
