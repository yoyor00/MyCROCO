"""Build the CROCO helper programs and the executables next to ``croco``."""

from __future__ import annotations

from SCons.Environment import Environment
from SCons.Node import Node

from .croco_build import CrocoBuild


class CrocoTools:
    """Build helper tools, ``partit``, ``ncjoin`` and ``ncrename``."""

    def __init__(self, croco_build: CrocoBuild):
        self.croco_build = croco_build
        self.env: Environment = croco_build.env
        cb = croco_build

        self.dep_cross_matrix = self._tool("cross_matrix.F", "cross_matrix")
        self.dep_cppcheck = self._tool("cppcheck.F", "cppcheck")
        self.dep_srcscheck = self._tool("srcscheck.F", "srcscheck")

        # The model source list contains these generated files.  Their
        # builders must exist before CrocoBin asks SCons to compile them.
        self.dep_check_switches1: Node = self.env.Command(
            source=[self.dep_cppcheck, "cppdefs.h", "cppdefs_dev.h"],
            target="check_switches1.F",
            action=[
                "cat cppdefs.h cppdefs_dev.h > mergcpp.txt",
                "./cppcheck",
            ],
        )
        self.dep_check_srcs: Node = self.env.Command(
            source=[self.dep_srcscheck, "Makefile"],
            target="check_srcs.F",
            action=["./srcscheck"],
        )

        # Staged names avoid clashes with the main CROCO objects.
        insert_node_o = cb.compile_fort(
            cb.preprocessor.apply(cb.copy_stage("insert_node.F", "tools")),
            rename_to_orig_filename=False,
        )
        lenstr_o = cb.compile_fort(
            cb.preprocessor.apply(cb.copy_stage("lenstr.F", "tools")),
            rename_to_orig_filename=False,
        )

        # partit.F is only preprocessed; ncjoin.F and ncrename.F also use mpc.
        self.dep_partit: Node = cb.link(
            input_objects=cb.compile_fort(cb.preprocessor.apply("partit.F"))
            + insert_node_o
            + lenstr_o,
            output_file="partit",
        )
        self.dep_ncjoin: Node = cb.link(
            input_objects=cb.compile_fort(cb.preproc_mpc("ncjoin.F"))
            + lenstr_o,
            output_file="ncjoin",
        )
        self.dep_ncrename: Node = cb.link(
            input_objects=cb.compile_fort(cb.preproc_mpc("ncrename.F"))
            + lenstr_o,
            output_file="ncrename",
        )

    def _tool(self, source: str, output: str) -> Node:
        return self.croco_build.compile_link(
            self.croco_build.preproc_mpc(source),
            output_file=output,
            LIBS=[],
        )
