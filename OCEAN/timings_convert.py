#!/usr/bin/env python3
"""Convert CROCO timing annotations into plain Fortran calls.

The input file is read once and the converted text is written once.  The
script is intentionally stateless: it does not need to know about other
source files and it never records region ids.  Annotated starts and named
stops are emitted with id ``-1``; the runtime module resolves the name
against its region table or open-region stack.

Only the ``!$croco timings ...`` form is recognized.  The optional trailing
text after a region name is treated as a comment.

Usage::

    python3 timings_convert.py INPUT OUTPUT
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# ``!$croco timings <keyword> ...`` with flexible blanks after ``!$`` and
# between the words.  A word boundary after the keyword prevents ``useful``
# from being parsed as ``use``.
ANNOTATION = re.compile(
    r"^(?P<indent>[ \t]*)!\$[ \t]*croco[ \t]+timings[ \t]+"
    r"(?P<keyword>use|include|init|start|stop|finalize)\b"
    r"(?P<rest>.*)$",
    re.IGNORECASE,
)

# Anything that looks like an annotation but uses a prefix we do not support.
WRONG_PREFIX = re.compile(
    r"^[ \t]*!(?:[ \t]*pragma[ \t]+croco\b|croco\b)",
    re.IGNORECASE,
)

# A supported prefix followed by an unknown keyword.
TIMING_PREFIX = re.compile(
    r"^[ \t]*!\$[ \t]*croco[ \t]+timings\b",
    re.IGNORECASE,
)

INTEGER = re.compile(r"^[+-]?\d+$")


class ConversionError(Exception):
    """An annotation could not be converted."""


def _parse_name(text: str, *, where: str) -> str:
    """Return a region name from the annotation payload.

    Names may be bare (no blanks) or quoted with single or double quotes.
    Trailing ``! ...`` text is a comment and is ignored.
    """
    text = text.split("!", 1)[0].strip()
    if not text:
        raise ConversionError(f"{where}: missing region name")

    if text[0] in ("'", '"'):
        quote = text[0]
        out: list[str] = []
        i = 1
        while i < len(text):
            char = text[i]
            if char == quote:
                if i + 1 < len(text) and text[i + 1] == quote:
                    out.append(quote)
                    i += 2
                    continue
                remainder = text[i + 1 :].strip()
                if remainder:
                    raise ConversionError(
                        f"{where}: unexpected text after quoted name: {remainder!r}"
                    )
                return "".join(out)
            out.append(char)
            i += 1
        raise ConversionError(f"{where}: unterminated quoted region name")

    name, _, extra = text.partition(" ")
    extra = extra.strip()
    if extra:
        raise ConversionError(
            f"{where}: unexpected text after region name: {extra!r} "
            "(quote names that contain blanks)"
        )
    return name


def _fortran_string(name: str) -> str:
    return "'" + name.replace("'", "''") + "'"


def _code_indent(indent: str) -> str:
    """Return indentation valid for both fixed- and free-form sources."""
    width = len(indent.expandtabs(8))
    if width < 6:
        # Fixed-form statements start in column 7; free-form ignores this.
        return " " * 6
    return " " * width


def convert_line(line: str, *, where: str) -> str:
    """Convert one source line, or return it unchanged."""
    match = ANNOTATION.match(line)
    if not match:
        if TIMING_PREFIX.match(line):
            raise ConversionError(
                f"{where}: unknown timings keyword; expected use, include, "
                "init, start, stop or finalize"
            )
        if WRONG_PREFIX.match(line):
            raise ConversionError(
                f"{where}: unsupported annotation prefix; use '!$croco timings ...'"
            )
        return line

    indent = _code_indent(match.group("indent"))
    keyword = match.group("keyword").lower()
    rest = match.group("rest")

    if keyword in ("use", "include"):
        return f"{indent}use timings\n"
    if keyword == "init":
        return f"{indent}call timings_init()\n"
    if keyword == "finalize":
        return f"{indent}call timings_finalize()\n"

    if keyword == "start":
        name = _parse_name(rest, where=where)
        return f"{indent}call timings_start({_fortran_string(name)}, -1)\n"

    # keyword == "stop"
    payload = rest.split("!", 1)[0].strip()
    if not payload:
        raise ConversionError(f"{where}: missing region id or name after 'stop'")
    if INTEGER.match(payload):
        return f"{indent}call timings_stop('', {int(payload)})\n"

    name = _parse_name(rest, where=where)
    return f"{indent}call timings_stop({_fortran_string(name)}, -1)\n"


def convert_text(text: str, *, source: str = "<input>") -> str:
    lines = text.splitlines(keepends=True)
    converted: list[str] = []
    for number, line in enumerate(lines, start=1):
        body = line.rstrip("\r\n")
        newline = line[len(body) :]
        converted.append(convert_line(body, where=f"{source}:{number}") + newline)
    return "".join(converted)


def convert_file(input_path: Path, output_path: Path) -> None:
    try:
        text = input_path.read_text(encoding="utf-8", errors="surrogateescape")
    except OSError as exc:
        raise ConversionError(f"cannot read {input_path}: {exc}") from exc

    converted = convert_text(text, source=str(input_path))

    try:
        output_path.write_text(
            converted, encoding="utf-8", errors="surrogateescape"
        )
    except OSError as exc:
        raise ConversionError(f"cannot write {output_path}: {exc}") from exc


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert '!$croco timings ...' annotations to Fortran calls."
    )
    parser.add_argument("input", type=Path, help="source file to read")
    parser.add_argument("output", type=Path, help="converted file to write")
    args = parser.parse_args(argv)

    try:
        convert_file(args.input, args.output)
    except ConversionError as exc:
        print(f"timings_convert.py: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
