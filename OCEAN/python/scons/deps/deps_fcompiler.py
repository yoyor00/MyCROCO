
from __future__ import annotations
import re
from typing import List
from SCons.Node.FS import FS
from SCons import Util
from .scanner_preprocess import Preprocessor


class FCompiler:
    """
    Find all `use` statements in Fortran files
    """

    NO_DEPS_MODS = [
        "openacc",
        "openmp",
        "omp_lib",
        "mpi",
        "mpi_f08",
        "mpi_f77",
        "mpi_f90",
        "perf_regions_fortran",
        "iso_fortran_env",
    ]

    def __init__(self, node: FS.File):
        self.node = node

    def getContents(self):
        if not hasattr(self.node.attributes, "contents"):
            self.node.attributes.contents = self.node.get_text_contents()
        return self.node.attributes.contents

    def getUseModules(self):
        """
        Return list of modules used via `use` statements.
        """

        if hasattr(self.node.attributes, "cache_use_modules"):
            a = self.node.attributes
            return a.cache_use_modules

        contents = self.getContents()

        #
        # Then, find all `use` statements
        #
        use_regex = (
            r"(?i)(?:^|;)\s*USE(?:\s+|(?:(?:\s*,\s*"
            r"(?:NON_)?INTRINSIC)?\s*::))\s*(\w+)"
        )
        cre_use = re.compile(use_regex, re.M)
        use_modules = cre_use.findall(contents)

        # Make all module names lower-case for case-insensitive comparison
        use_modules = [_.lower() for _ in use_modules]

        # Remove duplicates
        use_modules = Util.unique(use_modules)

        # Remove dependency for openacc, openmp, mpi, etc.
        use_modules = [m for m in use_modules if m not in self.NO_DEPS_MODS]

        # Remove all USE'd module names that are defined in the same file
        fcompiler = FCompiler(self.node)
        defined = fcompiler.getDefModules()
        use_modules = [m for m in use_modules if m.lower() not in defined]

        # We return module names as dependencies.
        # The caller can then append the `.mod` file extension
        # since one fortran file can output multiple .mod files and we like
        # to create a dependency to exactly this .mod file

        # Cache it
        a = self.node.attributes
        a.cache_use_modules = use_modules

        return use_modules

    def getIncludes(self, inc_paths: List[FS.Dir]):
        """
        Generate found dependencies which are also in `paths`

        Return tuple of
        'include files', 'used modules', 'defined modules'
        """

        # Did we cache it before?
        if hasattr(self.node.attributes, "cache_inc_dep_files"):
            a = self.node.attributes
            return a.cache_inc_dep_files

        contents = self.getContents()

        # Find Fortran ``include`` statements (not preprocessor includes).
        incl_regex = (
            r"""(?i)(?:^|['">]\s*;)\s*INCLUDE\s+"""
            r"""(?:\w+_)?[<"'](.+?)(?=["'>])"""
        )
        cre_incl = re.compile(incl_regex, re.M)
        inc_files = cre_incl.findall(contents)

        # We use the functionality of the precompiler includes
        # to sort out non-existing files.
        preproc_includes = Preprocessor(self.node)
        inc_dep_files = preproc_includes.gen_deps_exists(
            inc_paths, inc_files=inc_files
        )

        # Cache it
        a = self.node.attributes
        a.cache_inc_dep_files = inc_dep_files

        return inc_dep_files

    def getDefModules(self):
        """
        Generate found dependencies which are also in `paths`
        """

        # Did we cache it before?
        if hasattr(self.node.attributes, "cache_def_modules"):
            a = self.node.attributes
            return a.cache_def_modules

        contents = self.getContents()

        #
        # Find all `module` statements
        #
        def_regex = (
            r"""(?i)^\s*MODULE\s+"""
            r"""(?!PROCEDURE|SUBROUTINE|FUNCTION|PURE|ELEMENTAL)(\w+)"""
        )
        cre_def = re.compile(def_regex, re.M)
        def_modules = cre_def.findall(contents)

        def_modules = [_.lower() for _ in def_modules]

        # Remove dependency for openacc, openmp, mpi, etc.
        def_modules = [m for m in def_modules if m not in self.NO_DEPS_MODS]

        self.node.attributes.cache_def_modules = def_modules
        return def_modules
