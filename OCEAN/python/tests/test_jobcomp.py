#!/usr/bin/env python3
"""Unit tests for the CROCO build driver (``OCEAN/jobcomp.py``).

Standard library only, so they run with ``python -m unittest`` as well as with
``pytest`` (the rest of the repository uses pytest).

    python -m unittest discover -s OCEAN/python/tests -v
    python -m pytest OCEAN/python/tests
"""

import contextlib
import io
import os
import shutil
import subprocess
import stat
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parents[1]  # ``OCEAN/python``
OCEAN = HERE.parent
DRIVER = OCEAN / "jobcomp.py"
LEGACY = OCEAN / "jobcomp.sh"  # the historical bash script, untouched
LINK = OCEAN / "jobcomp"  # what people type: a symlink to jobcomp.py

sys.path.insert(0, str(HERE))

import components  # noqa: E402  (sibling package: the build features)
from jobcomp import config, context, copy, jobcomp_driver  # noqa: E402
from jobcomp import errors, features, makedefs, parameters, shell  # noqa: E402
from jobcomp import toolchain  # noqa: E402


def make_driver(argv=()):
    """A driver whose ``jobcomp.conf``/``jobcomp`` resolve to the checkout."""
    return jobcomp_driver.Jobcomp(list(argv), script_path=str(DRIVER))


def make_shell(dry_run: bool = False) -> shell.Shell:
    return shell.Shell(dry_run=dry_run)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def write(path: Path, text: str = "x\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def make_executable(path: Path, text: str) -> Path:
    write(path, text)
    mode = path.stat().st_mode
    path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def touch_time(path: Path, seconds: float) -> None:
    os.utime(path, (seconds, seconds))


class TempDirTest(unittest.TestCase):
    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)

    def path(self, *parts: str) -> Path:
        return self.dir.joinpath(*parts)


# ---------------------------------------------------------------------------
# copy helpers
# ---------------------------------------------------------------------------


class TestCopierFile(TempDirTest):
    def test_always_copies(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = self.path("dst", "a.F")
        copier = copy.Copier()
        self.assertTrue(copier.file(src, dst))
        self.assertEqual(dst.read_text(), "a\n")
        self.assertEqual(copier.stats.copied, 1)
        self.assertEqual(copier.stats.bytes_copied, 2)

    def test_directory_destination_keeps_name(self) -> None:
        src = write(self.path("a.F"), "a\n")
        out = self.path("out")
        out.mkdir()
        copy.Copier().file(src, out)
        self.assertTrue((out / "a.F").exists())

    def test_missing_source_raises(self) -> None:
        with self.assertRaises(FileNotFoundError):
            copy.Copier().file(self.path("nope"), self.path("out"))

    def test_identical_path_is_skipped(self) -> None:
        src = write(self.path("a.F"), "a\n")
        copier = copy.Copier()
        self.assertFalse(copier.file(src, src))
        self.assertEqual(copier.stats.skipped, 1)

    def test_overwrite_false_keeps_existing_file(self) -> None:
        src = write(self.path("src", "namelist"), "new\n")
        dst = write(self.path("dst", "namelist"), "old\n")
        copier = copy.Copier()
        self.assertFalse(copier.file(src, dst, overwrite=False))
        self.assertEqual(dst.read_text(), "old\n")
        self.assertEqual(copier.stats.skipped, 1)

    def test_newer_mode_skips_up_to_date_file(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = self.path("dst", "a.F")
        copy.Copier().file(src, dst)
        touch_time(dst, time.time() + 100)
        copier = copy.Copier(mode=copy.CopyMode.NEWER)
        self.assertFalse(copier.file(src, dst))
        self.assertEqual(copier.stats.skipped, 1)

    def test_newer_mode_copies_when_source_is_newer(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = write(self.path("dst", "a.F"), "a\n")
        touch_time(dst, time.time() - 1000)
        copier = copy.Copier(mode=copy.CopyMode.NEWER)
        self.assertTrue(copier.file(src, dst))
        self.assertEqual(copier.stats.copied, 1)

    def test_dry_run_changes_nothing(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = self.path("dst", "a.F")
        copier = copy.Copier(dry_run=True)
        self.assertTrue(copier.file(src, dst))
        self.assertFalse(dst.exists())
        self.assertEqual(copier.stats.copied, 1)

    def test_newer_is_the_default_mode(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = write(self.path("dst", "a.F"), "a\n")
        touch_time(src, time.time() - 1000)
        copier = copy.Copier()
        self.assertIs(copier.mode, copy.CopyMode.NEWER)
        self.assertFalse(copier.file(src, dst))
        self.assertEqual(copier.stats.skipped, 1)

    def test_always_mode_can_still_be_asked_for(self) -> None:
        src = write(self.path("src", "a.F"), "a\n")
        dst = write(self.path("dst", "a.F"), "a\n")
        touch_time(src, time.time() - 1000)
        copier = copy.Copier(mode=copy.CopyMode.ALWAYS)
        self.assertTrue(copier.file(src, dst))
        self.assertEqual(copier.stats.copied, 1)

    def test_copy_preserves_mtime(self) -> None:
        src = write(self.path("a.F"), "a\n")
        touch_time(src, time.time() - 4242)
        dst = self.path("dst", "a.F")
        copy.Copier().file(src, dst)
        self.assertAlmostEqual(
            dst.stat().st_mtime, src.stat().st_mtime, delta=1
        )


class TestCopierTree(TempDirTest):
    def test_mirrors_tree(self) -> None:
        write(self.path("src", "a.F"), "a\n")
        write(self.path("src", "sub", "b.F90"), "b\n")
        copier = copy.Copier()
        stats = copier.tree(self.path("src"), self.path("dst"))
        self.assertEqual(stats.copied, 2)
        self.assertTrue(self.path("dst", "sub", "b.F90").exists())

    def test_always_spec_wins_over_the_default_mode(self) -> None:
        src = write(self.path("src", "cppdefs.h"), "#define LOCAL 1\n")
        dst = write(self.path("dst", "cppdefs.h"), "#define SOURCE 1\n")
        touch_time(src, time.time() - 1000)
        ctx = context.BuildContext(
            paths=toolchain.Paths(
                run_dir=self.dir,
                source=self.path("src"),
                root=self.dir,
                scratch=self.path("dst"),
                agrif=self.path("AGRIF"),
            ),
            tools=None,
            detector=None,
            shell=make_shell(),
        )
        ctx.copier.apply(
            copy.CopySpec.glob(
                self.path("src"),
                self.path("dst"),
                "cppdefs.h",
                mode=copy.CopyMode.ALWAYS,
            )
        )
        self.assertEqual(dst.read_text(), "#define LOCAL 1\n")

    def test_delete_stale_removes_extra_files(self) -> None:
        write(self.path("src", "a.F"), "a\n")
        write(self.path("dst", "gone.F"), "old\n")
        copier = copy.Copier(delete_stale=True)
        stats = copier.tree(self.path("src"), self.path("dst"))
        self.assertFalse(self.path("dst", "gone.F").exists())
        self.assertEqual(stats.removed, 1)


class TestCopySpecs(TempDirTest):
    def make_context(self, **kwargs) -> context.BuildContext:
        scratch = self.path("Compile")
        scratch.mkdir(exist_ok=True)
        return context.BuildContext(
            paths=toolchain.Paths(
                run_dir=self.dir,
                source=self.path("src"),
                root=self.dir,
                scratch=scratch,
                agrif=self.path("AGRIF"),
            ),
            tools=None,
            detector=None,
            shell=make_shell(),
            **kwargs,
        )

    def test_contents_merges_subdirectories(self) -> None:
        # cp -r PISCES/* Compile/  then  cp -r PISCES/LIB/* Compile/
        write(self.path("PISCES", "rootfile.F90"), "r\n")
        write(self.path("PISCES", "LIB", "lib.h90"), "l\n")
        write(self.path("PISCES", "SED", "sed.F90"), "s\n")
        ctx = self.make_context()
        ctx.copier.apply(
            copy.CopySpec.contents(self.path("PISCES"), ctx.paths.scratch)
        )
        ctx.copier.apply(
            copy.CopySpec.contents(
                self.path("PISCES", "LIB"), ctx.paths.scratch
            )
        )
        scratch = ctx.paths.scratch
        self.assertTrue((scratch / "LIB" / "lib.h90").exists())
        self.assertTrue((scratch / "SED" / "sed.F90").exists())
        self.assertTrue((scratch / "rootfile.F90").exists())
        self.assertTrue((scratch / "lib.h90").exists())  # flattened copy
        self.assertFalse((scratch / "LIB" / "LIB").exists())
        self.assertFalse((scratch / "LIB" / "rootfile.F90").exists())

    def test_tree_keeps_directory_name(self) -> None:
        write(self.path("AGRIF", "LIB", "fortran.c"), "c\n")
        ctx = self.make_context()
        ctx.copier.apply(
            copy.CopySpec.tree(self.path("AGRIF"), ctx.paths.scratch)
        )
        self.assertTrue(
            (ctx.paths.scratch / "AGRIF" / "LIB" / "fortran.c").exists()
        )

    def test_glob(self) -> None:
        write(self.path("XIOS", "oxide.F"), "f\n")
        write(self.path("XIOS", "iodef.xml"), "x\n")
        ctx = self.make_context()
        ctx.copier.apply(
            copy.CopySpec.glob(
                self.path("XIOS"), ctx.paths.scratch, "*.F"
            )
        )
        self.assertTrue((ctx.paths.scratch / "oxide.F").exists())
        self.assertFalse((ctx.paths.scratch / "iodef.xml").exists())

    def test_optional_missing_source_is_ignored(self) -> None:
        ctx = self.make_context()
        ctx.copier.apply(
            copy.CopySpec.file(self.path("nope"), ctx.paths.scratch)
        )

    def test_required_missing_source_raises(self) -> None:
        ctx = self.make_context()
        with self.assertRaises(errors.BuildError):
            ctx.copier.apply(
                copy.CopySpec.file(
                    self.path("nope"), ctx.paths.scratch, optional=False
                )
            )


# ---------------------------------------------------------------------------
# toolchain
# ---------------------------------------------------------------------------


def make_toolchain(**overrides):
    kwargs = dict(
        shell=make_shell(),
        os_name="Linux",
        fc="gfortran",
        fc_from_env=False,
        mpif90="mpif90",
        fflags_override=None,
        nprocs=4,
        netcdf_lib="-L/usr/lib -lnetcdff -lnetcdf",
        netcdf_inc="-I/usr/include",
        netcdf_lib_default="-L/usr/lib -lnetcdff -lnetcdf",
        netcdf_inc_default="-I/usr/include",
        cc_override=None,
        cflags_override=None,
        make_override="make",
    )
    kwargs.update(overrides)
    return toolchain.Toolchain.resolve(**kwargs)


class TestToolchain(unittest.TestCase):
    def test_gfortran_defaults(self) -> None:
        tools = make_toolchain()
        self.assertEqual(tools.cft, "gfortran")
        self.assertIn("-fdefault-real-8", tools.fflags)
        self.assertEqual(tools.cpp, "cpp -traditional -DLinux")
        self.assertEqual(tools.ldflags, "-L/usr/lib -lnetcdff -lnetcdf")
        self.assertEqual(
            tools.cppflags, "-I/usr/include -ICROCOFILES/AGRIF_INC"
        )
        self.assertEqual(tools.nprocs, 4)

    def test_gfortran_versioned_name(self) -> None:
        tools = make_toolchain(fc="gfortran-12")
        self.assertIn("-std=legacy", tools.fflags)
        self.assertEqual(tools.cft, "gfortran-12")

    def test_ifc_uses_ifort(self) -> None:
        self.assertEqual(make_toolchain(fc="ifc").cft, "ifort")

    def test_fflags_override_wins(self) -> None:
        tools = make_toolchain(fflags_override="-O0 -g")
        self.assertTrue(tools.fflags.startswith("-O0 -g"))

    def test_unknown_compiler_from_default_is_fatal(self) -> None:
        with self.assertRaises(errors.BuildError):
            make_toolchain(fc="f95", fc_from_env=False)

    def test_unknown_compiler_from_user_is_ignored(self) -> None:
        with self.assertLogs(errors.LOG, level="WARNING"):
            tools = make_toolchain(fc="f95", fc_from_env=True)
        self.assertEqual(tools.cft, "f95")
        self.assertEqual(tools.cpp, "")

    def test_unsupported_os_is_fatal(self) -> None:
        with self.assertRaises(errors.BuildError):
            make_toolchain(os_name="Windows")

    def test_detected_netcdf_inc_is_used_for_fflags(self) -> None:
        # Bug-for-bug: the bash script uses $NETCDFINC (detected) in FFLAGS1
        # even when the user overrides CROCO_NETCDFINC.
        tools = make_toolchain(netcdf_inc="-I/custom/include")
        self.assertIn("-I/custom/include", tools.cppflags)
        self.assertNotIn("-I/custom/include", tools.fflags)
        self.assertIn("-I/usr/include", tools.fflags)

    def test_jean_zay_gets_define(self) -> None:
        os.environ["HOSTNAME"] = "jean-zay-1"
        self.addCleanup(os.environ.pop, "HOSTNAME", None)
        tools = make_toolchain(fc="nvfortran")
        self.assertIn("-DJEANZAY", tools.cpp)


# ---------------------------------------------------------------------------
# feature detection
# ---------------------------------------------------------------------------

FAKE_CPP = """#!/bin/sh
printf 'run\\n' >> "$FAKE_CPP_LOG"
for arg in "$@"; do
    case "$arg" in
        -*) ;;
        *) if [ -f "$arg" ]; then cat "$arg"; exit 0; fi ;;
    esac
done
cat | tee "${FAKE_CPP_STDIN:-/dev/null}"
"""


class TestFeatureDetector(TempDirTest):
    def setUp(self) -> None:
        super().setUp()
        self.bin = self.path("bin")
        make_executable(self.bin / "cpp", FAKE_CPP)
        self.log = self.path("cpp.log")
        os.environ["PATH"] = f"{self.bin}{os.pathsep}{os.environ['PATH']}"
        os.environ["FAKE_CPP_LOG"] = str(self.log)
        self.addCleanup(os.environ.pop, "FAKE_CPP_LOG", None)

    def test_single_preprocessor_run(self) -> None:
        source = self.path("src")
        write(source / "testkeys.F", "#ifdef AGRIF\nAGRIFISDEFINED\n#endif\n")
        detector = features.FeatureDetector(
            make_shell(), "cpp -P", self.path("Compile"), source
        )
        self.assertTrue(detector.has("agrifisdefined"))
        self.assertFalse(detector.has("xiosisdefined"))
        self.assertTrue(detector.has("agrifisdefined"))
        self.assertEqual(self.log.read_text().count("run"), 1)

    def test_missing_testkeys_is_not_fatal(self) -> None:
        detector = features.FeatureDetector(
            make_shell(), "cpp -P", self.path("Compile"), self.path("src")
        )
        self.assertFalse(detector.has("agrifisdefined"))

    def test_run_directory_header_wins_in_dry_run(self) -> None:
        source = self.path("src")
        write(source / "testkeys.F", '#include "cppdefs.h"\nMARKER\n')
        write(source / "cppdefs.h", "#define FROM_SOURCE 1\n")
        run_dir = self.path("run")
        write(run_dir / "cppdefs.h", "#define FROM_RUN_DIR 1\n")
        fed = self.path("stdin.txt")
        os.environ["FAKE_CPP_STDIN"] = str(fed)
        self.addCleanup(os.environ.pop, "FAKE_CPP_STDIN", None)
        detector = features.FeatureDetector(
            make_shell(), "cpp -P", self.path("Compile"), source, run_dir
        )
        self.assertTrue(detector.has("marker"))
        text = fed.read_text()
        self.assertIn("#define FROM_RUN_DIR 1", text)
        self.assertNotIn("FROM_SOURCE", text)
        self.assertNotIn('#include "cppdefs.h"', text)
        self.assertIn("MARKER", text)

    def test_no_preprocessor(self) -> None:
        detector = features.FeatureDetector(
            make_shell(), "", self.path("Compile"), self.path("src")
        )
        self.assertFalse(detector.has("agrifisdefined"))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


class TestCli(TempDirTest):
    def test_help(self) -> None:
        args = make_driver(["--help"]).parse_args()
        self.assertTrue(args.help)

    def test_copy_mode_values(self) -> None:
        args = make_driver(["--copy-mode", "newer"]).parse_args()
        self.assertEqual(args.copy_mode, "newer")
        with self.assertRaises(SystemExit):
            with contextlib.redirect_stderr(io.StringIO()):
                make_driver(["--copy-mode", "hash"]).parse_args()

    def test_jobs_is_an_int(self) -> None:
        self.assertEqual(make_driver(["--jobs", "8"]).parse_args().jobs, 8)

    def test_newer_is_the_default_copy_mode(self) -> None:
        source = self.path("src")
        source.mkdir()
        cfg = make_driver(["--src", str(source)])
        args = cfg.parse_args()
        self.assertIsNone(args.copy_mode)  # resolved by create_context
        ctx = cfg.create_context(args)
        self.assertIs(ctx.copy_mode, copy.CopyMode.NEWER)

    def test_unknown_option_exits(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                make_driver(["--nope"]).parse_args()
        self.assertEqual(caught.exception.code, 1)

    def test_pick_precedence(self) -> None:
        cfg = config.Config(str(DRIVER))
        os.environ["JOBC0MP_TEST"] = "env"
        self.addCleanup(os.environ.pop, "JOBC0MP_TEST", None)
        self.assertEqual(cfg.pick("cli", "JOBC0MP_TEST", "default"), "cli")
        self.assertEqual(cfg.pick(None, "JOBC0MP_TEST", "default"), "env")
        self.assertEqual(cfg.pick(None, "MISSING_VAR", "default"), "default")

    def test_environment_overrides_defaults(self) -> None:
        source = self.path("src")
        write(source / "testkeys.F", "")
        write(source / "cppdefs.h", "")
        os.environ.update(
            CROCO_SRC=str(source), CROCO_CFT1="ifort", NPROCS="7",
            CROCO_COPY_MODE="newer",
        )
        for name in ("CROCO_SRC", "CROCO_CFT1", "NPROCS", "CROCO_COPY_MODE"):
            self.addCleanup(os.environ.pop, name, None)
        args = make_driver().parse_args()
        ctx = make_driver().create_context(args)
        self.assertEqual(ctx.tools.fc, "ifort")
        self.assertEqual(ctx.tools.nprocs, 7)
        self.assertIs(ctx.copy_mode, copy.CopyMode.NEWER)

    def test_cli_overrides_environment(self) -> None:
        source = self.path("src")
        write(source / "testkeys.F", "")
        os.environ["CROCO_CFT1"] = "ifort"
        self.addCleanup(os.environ.pop, "CROCO_CFT1", None)
        args = make_driver(
            ["--src", str(source), "--fc", "gfortran", "--jobs", "3"]
        ).parse_args()
        ctx = make_driver().create_context(args)
        self.assertEqual(ctx.tools.fc, "gfortran")
        self.assertEqual(ctx.tools.nprocs, 3)

    def test_missing_source_directory(self) -> None:
        args = make_driver(["--src", str(self.path("nope"))]).parse_args()
        with self.assertRaises(errors.BuildError):
            make_driver().create_context(args)

    def test_fflags_accepts_a_value_starting_with_a_dash(self) -> None:
        # The bash script took the next word verbatim, so ``--fflags -O2``
        # must not be mistaken for an option.
        self.assertEqual(
            make_driver(["--fflags", "-O2"]).parse_args().fflags, "-O2"
        )
        args = make_driver(["--fflags", "-O2 -g", "--jobs", "2"]).parse_args()
        self.assertEqual(args.fflags, "-O2 -g")
        self.assertEqual(args.jobs, 2)

    def test_option_values_are_not_glued_onto_following_options(self) -> None:
        args = make_driver(["--fc", "gfortran", "--jobs", "3"]).parse_args()
        self.assertEqual((args.fc, args.jobs), ("gfortran", 3))

    def test_the_option_table_builds_the_parser_and_the_help(self) -> None:
        driver = make_driver()
        parser = driver.build_parser()
        destinations = {action.dest for action in parser._actions}
        self.assertEqual(
            destinations, {option.dest for option in jobcomp_driver.OPTIONS}
        )
        self.assertTrue(set(driver.VALUE_OPTIONS).isdisjoint(driver.FLAGS))
        help_text = parser.format_help()
        for option in jobcomp_driver.OPTIONS:
            for flag in option.flags:
                self.assertIn(flag, help_text)

    def test_invalid_job_count_is_a_build_error(self) -> None:
        source = self.path("src")
        source.mkdir()
        args = make_driver(["--src", str(source), "--jobs", "0"]).parse_args()
        with self.assertRaises(errors.BuildError):
            make_driver().create_context(args)


class TestLayout(unittest.TestCase):
    """``jobcomp`` links to the python driver, ``jobcomp.sh`` is the bash."""

    def test_driver_exists_and_is_executable(self) -> None:
        self.assertTrue(DRIVER.is_file())
        self.assertTrue(os.access(DRIVER, os.X_OK))
        self.assertTrue(DRIVER.read_text().startswith("#!"))

    def test_jobcomp_is_a_symlink_to_the_driver(self) -> None:
        self.assertTrue(LINK.is_symlink())
        self.assertEqual(os.readlink(LINK), "jobcomp.py")

    def test_legacy_script_is_a_regular_file(self) -> None:
        self.assertFalse(LEGACY.is_symlink())
        self.assertTrue(LEGACY.is_file())

    def test_python_parts_live_in_their_own_directory(self) -> None:
        module = OCEAN / "python" / "jobcomp" / "jobcomp_driver.py"
        self.assertTrue(module.is_file())

    def test_each_component_has_its_own_module(self) -> None:
        python = OCEAN / "python" / "components"
        self.assertTrue((python / "__init__.py").is_file())
        self.assertTrue((python / "base.py").is_file())
        self.assertTrue((python / "pipelines.py").is_file())
        for name in ("sources", "agrif", "xios", "pisces", "mustang",
                     "obstruction", "stogen", "mpi", "oasis", "openmp",
                     "openacc"):
            self.assertTrue((python / f"{name}.py").is_file(), name)

    def test_the_pipelines_use_the_components(self) -> None:
        names = {
            type(component).__name__
            for group in (
                components.COPY_COMPONENTS,
                components.RUN_DIR_COMPONENTS,
                components.CONFIGURE_BEFORE,
                components.PREPARE_COMPONENTS,
                components.CONFIGURE_AFTER,
            )
            for component in group
        }
        self.assertEqual(
            names,
            {
                "ModelSources", "AgrifComponent", "XiosComponent",
                "PiscesComponent", "MustangComponent", "ObstructionComponent",
                "StogenComponent", "MpiComponent", "OasisComponent",
                "OpenMpComponent", "OpenAccComponent",
            },
        )

    def test_legacy_script_exists_and_is_executable(self) -> None:
        self.assertTrue(LEGACY.is_file())
        self.assertTrue(os.access(LEGACY, os.X_OK))
        self.assertTrue(LEGACY.read_text().startswith("#!/bin/bash"))

    def test_legacy_script_keeps_the_tunables_other_scripts_sed(self) -> None:
        # oce_compile.sh still rewrites an old-style jobcomp: those lines must
        # stay in jobcomp.sh exactly as they were in the original.
        text = LEGACY.read_text()
        tunables = text.split("END OF USER'S MODIFICATIONS")[0].splitlines()
        for prefix in ("SOURCE1=", "FC=gfortran", "MPIF90=", "PRISM_ROOT_DIR=",
                       "XIOS_ROOT_DIR="):
            self.assertEqual(
                sum(line.startswith(prefix) for line in tunables), 1, prefix
            )
        self.assertIn("NPROCS", text)

    def test_defaults_keep_the_historical_names(self) -> None:
        for name in ("SOURCE1", "FC", "MPIF90", "MPILIB", "MPIINC",
                     "PRISM_ROOT_DIR", "XIOS_ROOT_DIR"):
            self.assertTrue(hasattr(config.Config, name), name)

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_legacy_script_is_valid_bash(self) -> None:
        result = subprocess.run(
            ["bash", "-n", str(LEGACY)], capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


class TestShellHandOver(TempDirTest):
    """``--sh`` forwards the command line to the unmodified bash script."""

    def test_shell_flag_is_stripped_from_the_passthrough(self) -> None:
        argv = make_driver(
            ["--sh", "--fc", "ifort", "--jobs", "2"]
        ).legacy_argv()
        self.assertEqual(argv, ["--fc", "ifort", "--jobs", "2"])

    def test_no_shell_flag_means_no_hand_over(self) -> None:
        self.assertIsNone(make_driver(["--fc", "ifort"]).legacy_argv())

    def test_shell_flag_used_as_a_value_is_not_a_hand_over(self) -> None:
        self.assertIsNone(make_driver(["--fflags", "--sh"]).legacy_argv())

    def test_legacy_script_is_found_next_to_the_driver(self) -> None:
        found = config.Config(str(DRIVER)).legacy_script()
        self.assertEqual(found, LEGACY.resolve())

    def test_missing_legacy_script_is_an_error(self) -> None:
        cfg = config.Config(str(self.path("empty", "jobcomp.py")))
        with self.assertRaises(errors.BuildError):
            cfg.legacy_script()

    @unittest.skipUnless(shutil.which("bash"), "bash not available")
    def test_shell_flag_runs_the_bash_script(self) -> None:
        result = subprocess.run(
            [sys.executable, str(DRIVER), "--sh", "--help"],
            capture_output=True, text=True, cwd=str(self.dir),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("JOBCOMP script to compile CROCO", result.stdout)
        self.assertIn(str(LEGACY), result.stdout)


class TestConfFile(TempDirTest):
    """``jobcomp.conf`` is the generation-agnostic replacement of the seds."""

    def config_for(self, body: str) -> config.Config:
        """A Config whose script directory holds ``body`` as jobcomp.conf."""
        script = write(self.path("jobcomp.py"), "")
        write(self.path("jobcomp.conf"), body)
        return config.Config(str(script))

    def test_missing_conf_file_is_empty(self) -> None:
        self.assertEqual(config.Config().load_conf(self.path("nope.conf")), {})

    def test_parses_comments_quotes_and_blank_lines(self) -> None:
        conf = write(
            self.path("jobcomp.conf"),
            "# comment\n\nSOURCE1=/tmp/croco/OCEAN\nFC=\"ifort\"\n"
            "XIOS_ROOT_DIR='$HOME/xios'\n",
        )
        values = config.Config().load_conf(conf)
        self.assertEqual(values["SOURCE1"], "/tmp/croco/OCEAN")
        self.assertEqual(values["FC"], "ifort")
        self.assertEqual(values["XIOS_ROOT_DIR"], "$HOME/xios")

    def test_values_are_expanded_when_used(self) -> None:
        os.environ["JC_TEST_OCE"] = "/tmp/jc"
        self.addCleanup(os.environ.pop, "JC_TEST_OCE", None)
        cfg = self.config_for("SOURCE1=${JC_TEST_OCE}/OCEAN\n")
        self.assertEqual(cfg.pick(None, "CROCO_SRC"), "/tmp/jc/OCEAN")

    def test_precedence_cli_then_environment_then_conf(self) -> None:
        cfg = self.config_for("SOURCE1=/from/conf\n")
        self.assertEqual(cfg.pick("/cli", "CROCO_SRC"), "/cli")
        self.assertEqual(cfg.pick(None, "CROCO_SRC"), "/from/conf")
        os.environ["CROCO_SRC"] = "/env"
        self.addCleanup(os.environ.pop, "CROCO_SRC", None)
        self.assertEqual(cfg.pick(None, "CROCO_SRC"), "/env")
        # a key absent from the file falls back to the hard-coded default
        self.assertEqual(cfg.pick(None, "CROCO_MPIF90"), config.Config.MPIF90)

    def test_the_names_are_the_historical_ones(self) -> None:
        # create_config.bash writes these names: keep them stable.
        self.assertEqual(config.Config.CONF_KEYS["CROCO_SRC"], "SOURCE1")
        self.assertEqual(config.Config.CONF_KEYS["CROCO_CFT1"], "FC")
        self.assertEqual(config.Config.CONF_KEYS["CROCO_PRISM_ROOT_DIR"],
                         "PRISM_ROOT_DIR")

    def test_the_default_conf_is_parsed_once(self) -> None:
        cfg = self.config_for("SOURCE1=/from/conf\nFC=ifort\n")
        with mock.patch.object(
            cfg, "_parse_conf", wraps=cfg._parse_conf
        ) as parse:
            self.assertEqual(cfg.pick(None, "CROCO_SRC"), "/from/conf")
            self.assertEqual(cfg.pick(None, "CROCO_CFT1"), "ifort")
        self.assertEqual(parse.call_count, 1)


class TestPrintFflags(unittest.TestCase):
    """``--print-fflags`` replaces the old ``sed -e s/-O3/-O2/`` on jobcomp."""

    def run_main(self, argv) -> str:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = jobcomp_driver.main(argv, script_path=str(DRIVER))
            self.assertEqual(code, 0)
        return out.getvalue().strip()

    def test_prints_the_gfortran_default(self) -> None:
        flags = self.run_main(["--print-fflags"])
        self.assertIn("-mcmodel=medium", flags)
        self.assertIn("-fdefault-real-8", flags)

    def test_follows_the_requested_compiler(self) -> None:
        flags = self.run_main(["--print-fflags", "--fc", "ifort"])
        self.assertIn("-fno-alias", flags)

    def test_flags_are_joined_with_single_blanks(self) -> None:
        tools = make_toolchain()
        tools.add_fflags("")  # an empty fragment must not leave a blank
        tools.add_fflags("-I/opt/inc")
        self.assertNotIn("  ", tools.fflags)
        self.assertTrue(tools.fflags.endswith("-I/opt/inc"))

    def test_user_flags_win(self) -> None:
        self.assertEqual(
            self.run_main(["--print-fflags", "--fflags", "-O1 -g"]), "-O1 -g"
        )


# ---------------------------------------------------------------------------
# Makedefs substitution
# ---------------------------------------------------------------------------


class TestMakedefs(TempDirTest):
    def make_context(self, dry_run: bool = False) -> context.BuildContext:
        scratch = self.path("Compile")
        scratch.mkdir(exist_ok=True)
        source = self.path("src")
        source.mkdir(exist_ok=True)
        return context.BuildContext(
            paths=toolchain.Paths(
                run_dir=self.dir,
                source=source,
                root=self.dir,
                scratch=scratch,
                agrif=self.path("AGRIF"),
            ),
            tools=make_toolchain(),
            detector=features.FeatureDetector(
                make_shell(), "", scratch, source, self.dir
            ),
            shell=make_shell(dry_run),
            dry_run=dry_run,
        )

    def test_tokens_are_substituted(self) -> None:
        ctx = self.make_context()
        write(
            ctx.paths.scratch / "Makedefs.generic",
            "FFLAGS=$(FFLAGS1)\nLD=$(LDFLAGS1)\nCPP=$(CPP1)\n"
            "FC=$(CFT1)\nCPPFLAGS=$(CPPFLAGS1)\n",
        )
        makedefs.Makedefs(ctx.tools).write(
            ctx.paths.scratch, dry_run=ctx.dry_run
        )
        text = (ctx.paths.scratch / "Makedefs").read_text()
        self.assertNotIn("$(", text)
        self.assertIn("FC=gfortran", text)
        self.assertIn("cpp -traditional -DLinux", text)

    def test_dry_run_writes_nothing(self) -> None:
        ctx = self.make_context(dry_run=True)
        write(ctx.paths.scratch / "Makedefs.generic", "FC=$(CFT1)\n")
        makedefs.Makedefs(ctx.tools).write(
            ctx.paths.scratch, dry_run=ctx.dry_run
        )
        self.assertFalse((ctx.paths.scratch / "Makedefs").exists())


# ---------------------------------------------------------------------------
# end to end (dry run)
# ---------------------------------------------------------------------------

FAKE_MAKE = """#!/bin/sh
for arg in "$@"; do
    if [ "$arg" = "-f" ]; then cat >/dev/null; echo cc; exit 0; fi
done
printf '%s\\n' "$*" >> "$FAKE_MAKE_PLAN"
exit 0
"""


class TestDryRun(TempDirTest):
    def setUp(self) -> None:
        super().setUp()
        self.bin = self.path("bin")
        make_executable(self.bin / "cpp", FAKE_CPP)
        self.make = make_executable(self.bin / "fakemake", FAKE_MAKE)
        os.environ["PATH"] = f"{self.bin}{os.pathsep}{os.environ['PATH']}"
        os.environ["FAKE_CPP_LOG"] = str(self.path("cpp.log"))
        self.plan = self.path("plan.txt")
        os.environ["FAKE_MAKE_PLAN"] = str(self.plan)
        self.addCleanup(os.environ.pop, "FAKE_CPP_LOG", None)
        self.addCleanup(os.environ.pop, "FAKE_MAKE_PLAN", None)

        self.source = self.path("repo", "OCEAN")
        write(self.source / "testkeys.F", "")
        write(self.source / "cppdefs.h", "")
        write(self.source / "Makefile", "all:\n")
        write(self.source / "Makedefs.generic", "FC=$(CFT1)\n")
        write(self.source / "croco.F", "program croco\nend\n")
        write(self.source.parent / "AGRIF" / "Makefile", "all:\n")
        write(self.source.parent / "PISCES" / "p2z.F90", "! bio\n")

    def run_main(self, *extra: str) -> int:
        out = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            with contextlib.redirect_stdout(out):
                code = jobcomp_driver.main(
                    [
                        "--src", str(self.source),
                        "--make", str(self.make),
                        "--dry-run",
                        *extra,
                    ],
                    script_path=str(DRIVER),
                )
        finally:
            os.chdir(cwd)
        self.output = out.getvalue()
        return code

    def test_dry_run_does_not_touch_the_disk(self) -> None:
        self.assertEqual(self.run_main(), 0)
        self.assertFalse(self.path("Compile").exists())
        self.assertIn("CROCO is OK", self.output)
        self.assertIn("would run: ", self.output)

    def test_dry_run_plan(self) -> None:
        self.assertEqual(self.run_main(), 0)
        self.assertIn("would run: ", self.output)
        self.assertIn("depend", self.output)
        self.assertIn("-j 1", self.output)
        self.assertFalse(self.plan.exists())  # nothing is executed

    def test_agrif_plan(self) -> None:
        write(self.source / "testkeys.F", "AGRIFISDEFINED\nMPIISDEFINED\n")
        self.assertEqual(self.run_main(), 0)
        self.assertIn(" => AGRIF activated", self.output)
        self.assertIn(" => MPI activated", self.output)
        self.assertIn("-C AGRIF", self.output)
        self.assertIn("CPPFLAGS=-DAGRIF_MPI", self.output)
        self.assertIn("FC=mpif90", self.output)
        # The AGRIF library is built before the main compile, on --jobs 1
        self.assertIn("-C AGRIF", self.output)

    def test_openmp_and_xios_flags(self) -> None:
        write(self.source / "testkeys.F", "OPENMP\nXIOSISDEFINED\n")
        self.assertEqual(self.run_main("--xios-dir", "/opt/xios"), 0)
        self.assertIn(" => XIOS activated", self.output)
        self.assertIn("would symlink", self.output)

    def test_makedefs_plan_mentions_flags(self) -> None:
        self.assertEqual(self.run_main("-v"), 0)
        self.assertIn("would write", self.output)

    def test_parameters_are_resolved_without_writing(self) -> None:
        # The API another build system (scons, cmake, ...) would use.
        write(self.source / "testkeys.F", "MPIISDEFINED\n")
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            params = make_driver(
                [
                    "--src", str(self.source),
                    "--make", str(self.make),
                    "--netcdf-lib", "-L/nc/lib -lnetcdf",
                    "--netcdf-inc", "-I/nc/include",
                ]
            ).parameters()
        finally:
            os.chdir(cwd)
        self.assertTrue(params.features["mpi"])
        self.assertFalse(params.features["agrif"])
        self.assertEqual(params.fc, "mpif90")  # MPI switches the compiler
        self.assertEqual(params.makedefs["CFT1"], params.fc)
        self.assertEqual(params.makedefs["FFLAGS1"], params.fflags)
        self.assertIn("Linux", params.defines)
        self.assertIn("-L/nc/lib", params.ldflags)
        self.assertIn("netcdf", params.libraries)
        self.assertIn("/nc/include", params.include_dirs)
        self.assertEqual(params.executables, ("croco", "partit", "ncjoin"))
        self.assertEqual(params.scratch, Path("Compile"))
        # nothing was copied, written or run
        self.assertFalse(self.path("Compile").exists())
        self.assertFalse(self.plan.exists())


class TestFlagClassification(unittest.TestCase):
    """``parameters`` turns the flag strings into lists for build systems."""

    def test_splits_defines_includes_libraries(self) -> None:
        defines, includes, libdirs, libs = parameters._classify(
            "cpp -traditional -DLinux -DIfort",
            "-I/usr/include -ICROCOFILES/AGRIF_INC",
            "-O2 -fopenmp",
            "-LAGRIF -lagrif -L/opt/nc/lib -lnetcdf /opt/xios/lib/libxios.a",
        )
        self.assertEqual(defines, ("Linux", "Ifort"))
        self.assertEqual(includes, ("/usr/include", "CROCOFILES/AGRIF_INC"))
        self.assertEqual(libdirs, ("AGRIF", "/opt/nc/lib"))
        self.assertEqual(
            libs, ("agrif", "netcdf", "/opt/xios/lib/libxios.a")
        )

    def test_accepts_separated_values(self) -> None:
        defines, includes, libdirs, libs = parameters._classify(
            "-D KEY=1 -I /inc -L /lib -l netcdf"
        )
        self.assertEqual(defines, ("KEY=1",))
        self.assertEqual(includes, ("/inc",))
        self.assertEqual(libdirs, ("/lib",))
        self.assertEqual(libs, ("netcdf",))


FAKE_MAKE_BUILD = """#!/bin/sh
for arg in "$@"; do
    if [ "$arg" = "-f" ]; then cat >/dev/null; echo cc; exit 0; fi
done
case "$1" in
    depend) exit 0 ;;
esac
for exe in croco partit ncjoin; do
    printf '#!/bin/sh\n' > "$exe"
    chmod +x "$exe"
done
printf 'croco ascii banner\n' > croco_ascii.txt
exit 0
"""


class TestBuild(TempDirTest):
    """A real (non dry-run) build with a fake compiler toolchain."""

    def setUp(self) -> None:
        super().setUp()
        self.bin = self.path("bin")
        make_executable(self.bin / "cpp", FAKE_CPP)
        self.make = make_executable(self.bin / "fakemake", FAKE_MAKE_BUILD)
        os.environ["PATH"] = f"{self.bin}{os.pathsep}{os.environ['PATH']}"
        os.environ["FAKE_CPP_LOG"] = str(self.path("cpp.log"))
        self.addCleanup(os.environ.pop, "FAKE_CPP_LOG", None)

        self.source = self.path("repo", "OCEAN")
        write(self.source / "testkeys.F", "AGRIFISDEFINED\n")
        write(self.source / "cppdefs.h", "#define AGRIF\n")
        write(self.source / "croco.F", "program croco\nend\n")
        write(self.source / "amr.in", "amr stuff\n")
        write(self.source / "private_scratch.h", "header\n")
        write(self.source / "Makefile", "all:\n")
        write(self.source / "Makedefs.generic", "FC=$(CFT1)\nLD=$(LDFLAGS1)\n")
        write(
            self.source / "Makedefs.generic.AGRIF",
            "FC=$(CFT1)\nLD=$(LDFLAGS1)\nAGRIF=1\n",
        )
        write(self.source.parent / "AGRIF" / "Makefile", "all:\n")
        write(self.source.parent / "AGRIF" / "conv", "conv source\n")
        write(self.source.parent / "PISCES" / "p2z.F90", "! bio\n")
        write(self.source.parent / "STOGEN" / "namelist_sto", "&sto\n/\n")

    def test_build_produces_executables(self) -> None:
        out = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            with contextlib.redirect_stdout(out):
                code = jobcomp_driver.main(
                    [
                        "--src", str(self.source),
                        "--make", str(self.make),
                        "--jobs", "2",
                    ],
                    script_path=str(DRIVER),
                )
        finally:
            os.chdir(cwd)
        output = out.getvalue()
        self.assertEqual(code, 0, output)
        self.assertIn("CROCO is OK", output)
        self.assertIn("croco ascii banner", output)
        for exe in ("croco", "partit", "ncjoin"):
            self.assertTrue((self.dir / exe).exists(), exe)
        self.assertFalse(self.path("Compile", "croco").exists())

        makedefs = self.path("Compile", "Makedefs").read_text()
        self.assertIn("FC=gfortran", makedefs)
        self.assertNotIn("$(", makedefs)

        crocofiles = self.path("Compile", "CROCOFILES")
        self.assertTrue((crocofiles / "amr.scrum").exists())
        self.assertTrue((crocofiles / "conv").exists())
        self.assertTrue((crocofiles / "private_scratch.h").exists())
        self.assertFalse((crocofiles / "private_scratch_AMR.h").exists())
        self.assertIn("amr stuff", (crocofiles / "amr.scrum").read_text())

    def test_run_directory_files_win_over_the_sources(self) -> None:
        write(self.dir / "cppdefs.h", "#define LOCAL 1\n")
        out = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            with contextlib.redirect_stdout(out):
                code = jobcomp_driver.main(
                    ["--src", str(self.source), "--make", str(self.make)],
                    script_path=str(DRIVER),
                )
        finally:
            os.chdir(cwd)
        self.assertEqual(code, 0, out.getvalue())
        scratch = self.path("Compile", "cppdefs.h").read_text()
        self.assertIn("#define LOCAL 1", scratch)
        self.assertNotIn("#define AGRIF", scratch)

    def test_keep_scratch_does_not_wipe(self) -> None:
        scratch = self.path("Compile")
        scratch.mkdir()
        keep = write(scratch / "keep.me", "keep\n")
        out = io.StringIO()
        cwd = os.getcwd()
        os.chdir(self.dir)
        try:
            with contextlib.redirect_stdout(out):
                code = jobcomp_driver.main(
                    [
                        "--src", str(self.source),
                        "--make", str(self.make),
                        "--keep-scratch",
                    ],
                    script_path=str(DRIVER),
                )
        finally:
            os.chdir(cwd)
        self.assertEqual(code, 0, out.getvalue())
        self.assertTrue(keep.exists())


class TestMainExitCodes(TempDirTest):
    def test_unknown_option(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as caught:
                jobcomp_driver.main(["--bogus"])
        self.assertEqual(caught.exception.code, 1)

    def test_help_returns_zero(self) -> None:
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(jobcomp_driver.main(["-h"]), 0)
        self.assertIn("jobcomp [OPTIONS]", out.getvalue())

    def test_missing_sources_returns_one(self) -> None:
        args = make_driver(["--src", str(self.path("nope"))]).parse_args()
        with self.assertRaises(errors.BuildError):
            make_driver().create_context(args)


if __name__ == "__main__":
    unittest.main()
