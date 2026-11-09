"""Preprocess ``.F``/``.F90`` sources into compilable ``.f``/``.f90``."""

from __future__ import annotations

from SCons.Node.FS import File

from .compile_pass_base import CommandPass


class CrocoPreprocessor(CommandPass):
    """Run the CROCO preprocessor on each source."""

    stage = "pre"
    action = (
        "$CPP -P $CPPFLAGS $_CPPDEFFLAGS "
        "$_CPPINCFLAGS $SOURCE > $TARGET"
    )

    def target(self, source: File) -> File:
        target = super().target(source)
        if source.name.endswith(".F90"):
            return self.helpers.file_change_ext(target, "f90")
        if source.name.endswith(".F"):
            return self.helpers.file_change_ext(target, "f")
        raise ValueError(
            f"Unknown file ending for '{source}', expected .F or .F90"
        )
