#!/usr/bin/env python3
"""Tests for the jobcomp -> SCons JSON build manifest."""

import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parents[1]  # ``OCEAN/python``

sys.path.insert(0, str(HERE))

from jobcomp import makefile, manifest  # noqa: E402
from jobcomp.copy import CopyMode  # noqa: E402
from jobcomp.flags import classify_flags  # noqa: E402
from jobcomp.settings import BuildBackend  # noqa: E402
from scons.config import ConfigError, CrocoBuildConfig  # noqa: E402


class TempDirTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)


class TestMakefileParser(TempDirTest):
    def test_reads_continued_variables_and_sources(self) -> None:
        makefile_path = self.dir / "Makefile"
        makefile_path.write_text(
            "# comment\n"
            " SRCS = main.F step.F\\\n"
            "        rhs3d.F\\\n"
            "\n"
            " SRCS90 = croco_namelist.F90\\\n"
            "          buffer.F90\n"
        )
        self.assertEqual(
            makefile.read_variable(makefile_path, "SRCS").split(),
            ["main.F", "step.F", "rhs3d.F"],
        )
        self.assertEqual(
            makefile.read_source_lists(makefile_path),
            {
                "F": ["main.F", "step.F", "rhs3d.F"],
                "F90": ["croco_namelist.F90", "buffer.F90"],
            },
        )

    def test_missing_variable_raises(self) -> None:
        makefile_path = self.dir / "Makefile"
        makefile_path.write_text("SRCS = main.F\n")
        with self.assertRaises(ValueError):
            makefile.read_variable(makefile_path, "SRCS90")


class TestBuildManifest(TempDirTest):
    def make_context(self) -> SimpleNamespace:
        source = self.dir / "OCEAN"
        scratch = self.dir / "Compile"
        run_dir = self.dir / "run"
        for path in (source, scratch, run_dir):
            path.mkdir()
        tools = SimpleNamespace(
            os="Linux",
            fc="gfortran",
            cft="gfortran",
            mpif90="mpif90",
            cc="cc",
            make="make",
            cpp="cpp -traditional -DLinux",
            cppflags="-I/usr/include -ICROCOFILES/AGRIF_INC",
            fflags="-O2 -I/usr/include",
            cflags="-O2",
            ldflags="-LAGRIF -lagrif -L/usr/lib -lnetcdff -lnetcdf",
            nprocs=4,
        )
        paths = SimpleNamespace(
            source=source,
            scratch=scratch,
            run_dir=run_dir,
            agrif=source.parent / "AGRIF",
        )
        return SimpleNamespace(
            tools=tools,
            paths=paths,
            features={"agrif": True, "mpi": False},
            backend=BuildBackend.SCONS,
            copy_mode=CopyMode.NEWER,
        )

    def test_manifest_round_trip(self) -> None:
        ctx = self.make_context()
        data = manifest.BuildManifest.from_context(
            ctx,
            targets=["croco", "partit", "ncjoin"],
            sources={"F": ["main.F"], "F90": ["buffer.F90"]},
        ).to_dict()
        self.assertEqual(data["schema_version"], manifest.SCHEMA_VERSION)
        self.assertEqual(data["build"]["backend"], "scons")
        self.assertEqual(data["build"]["jobs"], 4)
        self.assertEqual(
            data["toolchain"]["cpp"], ["cpp", "-traditional", "-DLinux"]
        )
        self.assertIn("-lnetcdff", data["flags"]["ldflags"])
        self.assertIn("netcdff", data["flags"]["libraries"])
        self.assertEqual(data["sources"]["F"], ["main.F"])

    def test_written_manifest_is_readable_by_scons(self) -> None:
        ctx = self.make_context()
        path = ctx.paths.scratch / manifest.MANIFEST_NAME
        manifest.BuildManifest.from_context(
            ctx,
            targets=["croco", "partit", "ncjoin"],
            sources={"F": ["main.F"], "F90": ["buffer.F90"]},
        ).write(path)
        config = CrocoBuildConfig.load(path)
        self.assertEqual(config.build.backend, "scons")
        self.assertEqual(config.toolchain.fc_invoked, "gfortran")
        self.assertEqual(config.flags.cppflags[0], "-I/usr/include")
        self.assertEqual(config.sources.fixed, ("main.F",))
        self.assertEqual(config.sources.free, ("buffer.F90",))

    def test_invalid_schema_version_is_rejected(self) -> None:
        path = self.dir / "croco_build.json"
        path.write_text(json.dumps({"schema_version": 999}))
        with self.assertRaises(ConfigError):
            CrocoBuildConfig.load(path)

    def test_flag_classification(self) -> None:
        defines, includes, libdirs, libs = classify_flags(
            "cpp -traditional -DLinux",
            "-I/usr/include -ICROCOFILES/AGRIF_INC",
            "-O2",
            "-LAGRIF -lagrif -L/usr/lib -lnetcdff -lnetcdf",
        )
        self.assertEqual(defines, ("Linux",))
        self.assertEqual(includes, ("/usr/include", "CROCOFILES/AGRIF_INC"))
        self.assertEqual(libdirs, ("AGRIF", "/usr/lib"))
        self.assertEqual(libs, ("agrif", "netcdff", "netcdf"))


if __name__ == "__main__":
    unittest.main()
