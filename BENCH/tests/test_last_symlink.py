#!/bin/env python3
##########################################################
#  CROCO build system, under CeCILL-C
#  From Sébastien Valat (INRIA & LJK) - 2023
#  CROCO website : http://www.croco-ocean.org
##########################################################

##########################################################
"""Tests for the `results/last` symlink handling of the benchmarking script."""

##########################################################
# python
import os

# internal
from libs.helpers import replace_symlink


##########################################################
def test_replace_symlink_creates_a_link(tmp_path):
    (tmp_path / "run-1").mkdir()
    link = tmp_path / "last"

    replace_symlink("run-1", str(link))

    assert link.is_symlink()
    assert os.readlink(link) == "run-1"


##########################################################
def test_replace_symlink_replaces_an_older_link(tmp_path):
    (tmp_path / "run-1").mkdir()
    (tmp_path / "run-2").mkdir()
    link = tmp_path / "last"
    os.symlink("run-1", link)

    replace_symlink("run-2", str(link))

    assert os.readlink(link) == "run-2"


##########################################################
def test_replace_symlink_replaces_a_dangling_link(tmp_path):
    """A removed run directory leaves a dangling link: os.path.exists is False."""
    link = tmp_path / "last"
    os.symlink("run-that-was-deleted", link)
    assert not os.path.exists(link)
    assert os.path.lexists(link)
    (tmp_path / "run-3").mkdir()

    # must not raise FileExistsError
    replace_symlink("run-3", str(link))

    assert os.readlink(link) == "run-3"
    assert (tmp_path / "run-3").is_dir()


##########################################################
def test_replace_symlink_replaces_a_directory(tmp_path):
    link = tmp_path / "last"
    link.mkdir()
    (link / "leftover.txt").write_text("x")

    replace_symlink("run-4", str(link))

    assert link.is_symlink()
    assert os.readlink(link) == "run-4"
