"""Run ``change_loops.py`` over fixed-form sources."""

from __future__ import annotations

from .compile_pass_base import CommandPass


class CrocoChangeLoops(CommandPass):
    """Apply the loop-reordering source pass."""

    stage = "cl"
    action = ['python3 ./change_loops.py "$SOURCE" "$TARGET"']
