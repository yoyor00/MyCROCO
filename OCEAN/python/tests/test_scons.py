#!/usr/bin/env python3
"""Tests for the optional SCons backend and its jobcomp integration."""

import contextlib
import importlib.util
import io
import os
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

HERE = Path(__file__).resolve().parents[1]  # ``OCEAN/python``
OCEAN = HERE.parent
DRIVER = OCEAN / "jobcomp.py"

sys.path.insert(0, str(HERE))

import scons.backend as backend  # noqa: E402
from jobcomp import errors, jobcomp_driver, settings  # noqa: E402


def make_driver(argv=()):
    return jobcomp_driver.Jobcomp(list(argv), script_path=str(DRIVER))


def write(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


class TestSconsBackend(unittest.TestCase):
    def test_command_uses_the_packaged_sconstruct(self) -> None:
        ctx = SimpleNamespace(tools=SimpleNamespace(nprocs=4))
        with mock.patch.object(backend, "scons_argv", return_value=["scons"]):
            argv = backend.command(ctx)
        self.assertEqual(
            argv,
            [
                "scons",
                "-f",
                str(backend.SCONS_SCRIPT),
                "-j",
                "4",
                "croco",
                "partit",
                "ncjoin",
            ],
        )
        self.assertTrue(backend.SCONS_SCRIPT.is_file())

    def test_dry_run_works_without_scons(self) -> None:
        ctx = SimpleNamespace(tools=SimpleNamespace(nprocs=1))
        with mock.patch.object(backend, "scons_argv", return_value=None):
            argv = backend.command(ctx, require=False)
        self.assertEqual(argv[0], "scons")

    def test_missing_scons_is_fatal(self) -> None:
        ctx = SimpleNamespace(tools=SimpleNamespace(nprocs=1))
        with mock.patch.object(backend, "scons_argv", return_value=None):
            with self.assertRaises(errors.BuildError):
                backend.command(ctx, require=True)

    def test_environment_points_to_package_and_manifest(self) -> None:
        scratch = Path("/tmp/Compile")
        ctx = SimpleNamespace(paths=SimpleNamespace(scratch=scratch))
        env = backend.environment(ctx)
        self.assertEqual(env["CROCO_SCONS_PACKAGE"], str(backend.PACKAGE_DIR))
        self.assertEqual(
            env["CROCO_BUILD_CONFIG"], str(scratch / "croco_build.json")
        )


@unittest.skipUnless(
    importlib.util.find_spec("SCons") is not None, "SCons is not installed"
)
class TestSconsPackage(unittest.TestCase):
    def test_public_classes_are_importable(self) -> None:
        import scons

        for name in (
            "CrocoEnv",
            "CrocoBuild",
            "CrocoTools",
            "CrocoBin",
            "CrocoBuildConfig",
        ):
            self.assertTrue(callable(getattr(scons, name)), name)


class TestSconsJobcomp(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.source = self.dir / "OCEAN"
        write(self.source / "testkeys.F", "")
        write(self.source / "cppdefs.h", "")

    def context(self, argv):
        driver = make_driver(["--src", str(self.source), *argv])
        return driver.create_context(driver.parse_args())

    def test_scons_flag_selects_the_backend(self) -> None:
        ctx = self.context(["--scons"])
        self.assertIs(ctx.backend, settings.BuildBackend.SCONS)

    def test_backend_option_selects_the_backend(self) -> None:
        ctx = self.context(["--backend", "scons"])
        self.assertIs(ctx.backend, settings.BuildBackend.SCONS)

    def test_default_backend_is_make(self) -> None:
        ctx = self.context([])
        self.assertIs(ctx.backend, settings.BuildBackend.MAKE)

    def test_environment_selects_the_backend(self) -> None:
        os.environ["CROCO_BUILDER"] = "scons"
        self.addCleanup(os.environ.pop, "CROCO_BUILDER", None)
        ctx = self.context([])
        self.assertIs(ctx.backend, settings.BuildBackend.SCONS)

    def test_jobcomp_conf_selects_the_backend(self) -> None:
        script = write(self.dir / "jobcomp.py", "")
        write(self.dir / "jobcomp.conf", "BUILDER=scons\n")
        driver = jobcomp_driver.Jobcomp(
            ["--src", str(self.source)], script_path=str(script)
        )
        ctx = driver.create_context(driver.parse_args())
        self.assertIs(ctx.backend, settings.BuildBackend.SCONS)

    def test_dry_run_prints_the_scons_command(self) -> None:
        out = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            with contextlib.redirect_stdout(out):
                code = jobcomp_driver.main(
                    [
                        "--src",
                        str(self.source),
                        "--scons",
                        "--dry-run",
                        "--jobs",
                        "3",
                        "--netcdf-lib",
                        "-L/nc/lib -lnetcdf",
                        "--netcdf-inc",
                        "-I/nc/include",
                    ],
                    script_path=str(DRIVER),
                )
        finally:
            os.chdir(cwd)
        self.assertEqual(code, 0)
        text = out.getvalue()
        self.assertIn("would run: ", text)
        self.assertIn(f"-f {backend.SCONS_SCRIPT}", text)
        self.assertIn("-j 3 croco partit ncjoin", text)


if __name__ == "__main__":
    unittest.main()
