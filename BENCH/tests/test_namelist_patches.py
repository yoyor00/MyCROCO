#!/bin/env python3
##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""
Unit tests for the namelist patching used by --debug and --restart.

Nothing is compiled here: the tests only check that BENCH patches the
`croco.nml` namelists of the test cases without touching anything else.
They cover the regression that made `--debug` abort with
"IndexError: list index out of range" (the old code was still parsing the
former `croco.in` card format) for every test case using
`&croco_use_calendar`.
"""

##########################################################
# python
import os
import shutil
from datetime import timedelta
from math import ceil
from pathlib import Path

# never open a window from the plot helpers imported by libs.croco
os.environ.setdefault("MPLBACKEND", "Agg")

# namelist
import f90nml

# tests
import pytest

# internal
from libs.croco import Croco
from libs.helpers import parse_datetime


##########################################################
REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_CASES = REPO_ROOT / "TEST_CASES"

# all the test cases of the repository defining &croco_use_calendar, i.e. the
# ones that used to crash in --debug
CALENDAR_NML_FILES = [
    "DUNE/croco_DUNE.nml",
    "ESTUARY/croco_ESTUARY.nml",
    "REALIST_EXAMPLES/croco_COASTAL.nml",
    "SED_TOY/croco_SED_TOY_FLOC_0D.nml",
    "SEAGRASS/croco_SEAGRASS.nml",
    "TIDAL_FLAT/croco_TIDAL_FLAT.nml",
]

# a case without &croco_use_calendar
NO_CALENDAR_NML_FILE = "SANDBAR/croco_SANDBAR.nml"


##########################################################
def make_croco(tmp_path: Path, nml_relpath: str, restarted: bool = False) -> Croco:
    """Build a Croco instance patching a copy of the given test-case nml."""
    shutil.copyfile(TEST_CASES / nml_relpath, tmp_path / "croco.nml")
    croco = Croco.__new__(Croco)
    croco.dirname = str(tmp_path)
    croco.croco_inputfile = "croco.nml"
    croco.restarted = restarted
    return croco


def expected_end_date(nml, ntimes: int) -> str:
    dt = float(nml["croco_time_stepping"]["dt"])
    duration = ceil(max(dt * ntimes, 1.0))
    start = parse_datetime(nml["croco_use_calendar"]["start_date"])
    return (start + timedelta(seconds=duration)).strftime("%Y-%m-%d %H:%M:%S")


def expected_dt_his(nml, ntimes: int) -> float:
    dt = float(nml["croco_time_stepping"]["dt"])
    duration = ceil(max(dt * ntimes, 1.0))
    return max(dt / 3600.0, duration / (ntimes * 3600))


##########################################################
@pytest.mark.parametrize("nml_relpath", CALENDAR_NML_FILES)
def test_debug_patches_calendar_case(tmp_path, nml_relpath):
    """--debug must shorten the run of every calendar-based test case."""
    croco = make_croco(tmp_path, nml_relpath)
    nml_file = tmp_path / "croco.nml"
    original_text = nml_file.read_text()
    before = f90nml.read(nml_file)

    croco.apply_debug_patches()

    after = f90nml.read(nml_file)
    assert after["croco_time_stepping"]["ntimes"] == 6
    assert after["croco_history"]["nwrt"] == 1
    assert after["croco_use_calendar"]["end_date"] == expected_end_date(before, 6)
    assert after["croco_use_calendar"]["dt_his"] == pytest.approx(
        expected_dt_his(before, 6)
    )

    # untouched values must stay untouched (f90nml.patch only rewrites the
    # patched values, a read+write would have turned 0.0d0 into 0.0)
    assert after["croco_time_stepping"]["dt"] == before["croco_time_stepping"]["dt"]
    if "0.0d0" in original_text:
        assert "0.0d0" in nml_file.read_text()

    # no temporary file left behind
    assert not list(tmp_path.glob("*.nml.tmp"))


def test_debug_patches_without_calendar(tmp_path):
    """--debug must work (and not invent sections) for a plain case."""
    croco = make_croco(tmp_path, NO_CALENDAR_NML_FILE)

    croco.apply_debug_patches()

    after = f90nml.read(tmp_path / "croco.nml")
    assert after["croco_time_stepping"]["ntimes"] == 6
    assert after["croco_history"]["nwrt"] == 1
    assert "croco_use_calendar" not in after


def test_change_nml_skips_unknown_keys(tmp_path):
    """An absent section/key must be reported as skipped, not written."""
    croco = make_croco(tmp_path, NO_CALENDAR_NML_FILE)
    nml_file = tmp_path / "croco.nml"
    before = f90nml.read(nml_file)

    assert croco.change_nml("croco.nml", "croco_diags_ek", "nwrtdiags_ek", 1) is False
    assert croco.change_nml("croco.nml", "croco_time_stepping", "ntimes", 6) is True

    after = f90nml.read(nml_file)
    assert "croco_diags_ek" not in after
    assert after["croco_time_stepping"]["ntimes"] == 6
    assert set(after.keys()) == set(before.keys())


@pytest.mark.parametrize(
    "nml_relpath", [NO_CALENDAR_NML_FILE, "TIDAL_FLAT/croco_TIDAL_FLAT.nml"]
)
def test_restart_patches(tmp_path, nml_relpath):
    """--restart must prepare both the forwarding run and the restarted one."""
    croco = make_croco(tmp_path, nml_relpath, restarted=True)
    nml_file = tmp_path / "croco.nml"
    before = f90nml.read(nml_file)

    croco.apply_restart_patches()

    nml = f90nml.read(nml_file)
    assert nml["croco_history"]["ldefhis"] is False
    assert nml["croco_time_stepping"]["ntimes"] == 3
    assert nml["croco_restart"]["nrst"] == 3
    assert nml["croco_restart"]["nrpfrst"] == 0
    assert nml["croco_restart"]["rstname"] == "croco_restart.nc"
    if "croco_use_calendar" in before:
        assert nml["croco_use_calendar"]["end_date"] == expected_end_date(before, 3)
        assert nml["croco_use_calendar"]["dt_rst"] == pytest.approx(
            ceil(max(float(before["croco_time_stepping"]["dt"]) * 3, 1.0)) / 3600.0
        )

    nml_rst = f90nml.read(tmp_path / "croco.nml_rst")
    assert nml_rst["croco_time_stepping"]["ntimes"] == 3
    assert nml_rst["croco_initial"]["nrrec"] == 2
    assert nml_rst["croco_initial"]["ininame"] == "croco_restart.nc"
