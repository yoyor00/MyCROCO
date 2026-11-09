"""SCons environment for CROCO, loaded from the JSON build manifest.

``jobcomp`` writes ``Compile/croco_build.json`` before it invokes SCons, so the
SCons backend uses exactly the compiler, flags, libraries and source lists
resolved by the regular build driver.  SCons does not inherit the process
environment: every variable that must reach the compiler is set explicitly.
"""

from __future__ import annotations

import os
from typing import List, Sequence, Tuple

from SCons.Environment import Environment
from SCons.Script import Dir

from .config import ConfigError, CrocoBuildConfig
from .deps.scanner_fortran import ScannerFortran
from .deps.scanner_preprocess import ScannerPreprocess


def _without(
    values: Sequence[str],
    prefixes: Tuple[str, ...],
    suffixes: Tuple[str, ...] = (),
) -> List[str]:
    """Drop one-token and two-token options such as ``-I``/``-L``/``-l``."""
    result = []
    index = 0
    while index < len(values):
        token = values[index]
        if token in prefixes and index + 1 < len(values):
            index += 2
            continue
        if token.startswith(prefixes) or token.endswith(suffixes):
            index += 1
            continue
        result.append(token)
        index += 1
    return result


class CrocoEnv:
    """Set up the SCons environment used to build CROCO."""

    def __init__(self, config: CrocoBuildConfig):
        self.config = config
        if config.features.get("agrif"):
            raise ConfigError(
                "The SCons backend does not yet support AGRIF builds; "
                "use the default make backend."
            )
        self.env = Environment(tools=["default"])
        self.env["CPPPATH"] = ["."]
        self.load_basic_os_env()
        self.setup_env_vars()
        self.setup_workarounds()
        self.setup_scanners()
        if os.environ.get("CROCO_SCONS_VERBOSE") == "1":
            self.print_config()

        if self.env["CPPPATH"][0] != ".":
            raise RuntimeError(
                "the first entry in CPPPATH must be the current directory"
            )

    # -- manifest values ---------------------------------------------------
    def setup_env_vars(self) -> None:
        """Translate the manifest into SCons construction variables."""
        tool = self.config.toolchain
        flags = self.config.flags

        cpp = list(tool.cpp)
        self.env["CPP"] = cpp[0] if cpp else "cpp"
        self.env["CPPFLAGS"] = _without(
            cpp[1:] + list(flags.cppflags), ("-I", "-D")
        )
        self.env["CPPDEFINES"] = list(flags.defines)
        self.env["CFT"] = tool.fc_invoked
        self.env["FFLAGS"] = _without(list(flags.fflags), ("-I",))
        self.env["FFLAGS90"] = _without(list(flags.fflags90), ("-I",))
        self.env["LDR"] = tool.ldr
        self.env["CC"] = tool.cc
        self.env["CFLAGS"] = list(flags.cflags)

        self.env["CPPPATH"] += list(flags.include_dirs)
        self.env["LIBPATH"] = list(flags.library_dirs)
        self.env["LIBS"] = list(flags.libraries)
        self.env["LDFLAGS"] = _without(
            list(flags.ldflags), ("-L", "-l"), suffixes=(".a",)
        )

    def setup_workarounds(self) -> None:
        """Compiler-specific workarounds kept from the original example."""
        # jobcomp recommends -std=legacy, but flang-new does not support it.
        if "flang-new" in self.env["CFT"]:
            for name in ("FFLAGS", "FFLAGS90"):
                if "-std=legacy" in self.env[name]:
                    self.env[name].remove("-std=legacy")

    # -- OS environment ----------------------------------------------------
    def load_basic_os_env(self) -> None:
        """Copy the environment variables selected by jobcomp."""
        values = self.config.environment
        for name, value in values.items():
            self.env["ENV"][name] = value

        for libdir in values.get("LD_LIBRARY_PATH", "").split(os.pathsep):
            if os.path.isdir(libdir):
                self.env.Append(LIBPATH=libdir)

        for incdir in values.get("C_INCLUDE_PATH", "").split(os.pathsep):
            if os.path.isdir(incdir):
                self.env.Append(CPPPATH=incdir)

    # -- dependency scanners ----------------------------------------------
    def setup_scanners(self) -> None:
        """Add the CROCO preprocessor and Fortran dependency scanners."""
        src_paths = tuple(Dir(path) for path in self.env.get("CPPPATH", []))

        def path_function(env, dir, args, all_paths):
            return src_paths

        self.env["SCANNERS"] = [
            ScannerPreprocess(path_function),
            ScannerFortran(path_function),
        ] + self.env.get("SCANNERS", [])

    # -- diagnostics -------------------------------------------------------
    def print_config(self) -> None:
        """Print the resolved environment once, like the original example."""
        print("=" * 80)
        print("| Build environment variables:")
        print("=" * 80)
        for name in ("CPP", "CFT", "LDR", "LDFLAGS", "CC", "CFLAGS"):
            print(f"{name:>7}: '{self.env.get(name, '')}'")
        print("-" * 80)
        print(f"CPPPATH: '{self.env['CPPPATH']}'")
        print(f"   LIBS: '{self.env['LIBS']}'")
        print(f"LIBPATH: '{self.env['LIBPATH']}'")
        print("=" * 80)
