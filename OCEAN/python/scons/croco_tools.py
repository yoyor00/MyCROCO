"""Build the CROCO helper programs and the executables next to ``croco``."""

from __future__ import annotations

from SCons.Environment import Environment
from SCons.Node import Node

from .croco_build import CrocoBuild


class CrocoTools:
    """Build the helper programs, ``partit``, ``ncjoin`` and ``ncrename``."""

    def __init__(self, croco_build: CrocoBuild):
        self.croco_build = croco_build
        self.env: Environment = croco_build.env
        cb = croco_build

        # Current Makefile tools.  checkkwds/setup_kwds/checkdefs were removed
        # from CROCO; the SCons scanners find the dependencies directly.
        self.dep_cross_matrix: Node = cb.compile_link(
            cb.preproc_mpc("cross_matrix.F"),
            output_file="cross_matrix",
            LIBS=[],
        )
        self.dep_cppcheck: Node = cb.compile_link(
            cb.preproc_mpc("cppcheck.F"),
            output_file="cppcheck",
            LIBS=[],
        )
        self.dep_srcscheck: Node = cb.compile_link(
            cb.preproc_mpc("srcscheck.F"),
            output_file="srcscheck",
            LIBS=[],
        )

        # The model source list contains generated files check_srcs.F and
        # check_switches1.F.  Their builders must exist before CrocoBin asks
        # SCons to compile those sources.
        self.dep_check_switches1: Node = self.env.Command(
            source=[
                self.dep_cppcheck,
                "cppdefs.h",
                "cppdefs_dev.h",
            ],
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

        # Build these two objects under a staged name so they do not clash
        # with the main CROCO objects.
        insert_node_o = cb.compile_fort(
            cb.preprocessor.apply(cb.copy_stage("insert_node.F", "tools")),
            rename_to_orig_filename=False,
        )
        lenstr_o = cb.compile_fort(
            cb.preprocessor.apply(cb.copy_stage("lenstr.F", "tools")),
            rename_to_orig_filename=False,
        )

        # partit.F is only preprocessed in the current Makefile (no mpc),
        # whereas ncjoin.F and ncrename.F also go through mpc.
        partit_o = cb.compile_fort(cb.preprocessor.apply("partit.F"))
        self.dep_partit: Node = cb.link(
            input_objects=partit_o + insert_node_o + lenstr_o,
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

        self.tools: Node = self.env.Command(
            source=[
                cb.dep_mpc,
                self.dep_cross_matrix,
                self.dep_cppcheck,
                self.dep_srcscheck,
                self.dep_check_switches1,
                self.dep_check_srcs,
                self.dep_partit,
                self.dep_ncjoin,
                self.dep_ncrename,
            ],
            target="tools",
            action="",
        )
