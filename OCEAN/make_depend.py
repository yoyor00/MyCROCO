#!/usr/bin/env python3
# ======================================================================
# CROCO is derived from the ROMS-AGRIF branch of ROMS.
# ROMS-AGRIF was developed by IRD and Inria. CROCO also inherits
# from the UCLA branch (Shchepetkin et al.) and the Rutgers
# University branch (Arango et al.), both under MIT/X style license.
# Copyright (C) 2005-2026 CROCO Development Team
# License: CeCILL-2.1 - see LICENSE.txt
#
# CROCO website : https://www.croco-ocean.org
# ======================================================================
#
# make_depend.py: analyzer of dependencies to build table of
#               dependencies "Make.depend". Pure-Python, stdlib-only
#               replacement for cross_matrix.F.
"""
Scans Fortran/CPP source (.F, .F90, and any headers they transitively
#include/INCLUDE) and writes Make.depend, the same way cross_matrix.F
did, but as a plain script with no compile step, no fixed-size
limits, and a real module-name registry instead of guessing that a
module's name matches its defining file's name.

Usage: ./make_depend.py [-o OUTPUT] [--no-fast-filter] file.F90 file.F ...
(run from the directory make depend is invoked from, e.g. OCEAN/;
files are catalogued by basename, so this assumes the flat layout
jobcomp assembles in Compile/, which is the only supported one).

Key design decisions:
- USE'd modules resolve to a direct dependency on the defining
  file's .o; make's own graph handles further transitivity through
  that file's own Make.depend entry. #include'd headers are not
  build targets, so their full transitive closure is computed and
  listed explicitly instead (see resolve_dependencies).
- A file reached via #include keeps its own literal name (CROCO
  commonly splices in raw .F source this way, e.g. biology.F
  #include-ing bio_BioEBUS.F -- not a separate build target); only a
  file reached via USE is rewritten to its .o form.
- A missing #include'd file is warned about and excluded, matching
  cross_matrix.F's own EXCLUDE semantics: it is dropped from any
  emitted dependency line entirely, not just logged, since a still-
  emitted-but-unsatisfiable prerequisite makes `make` fail outright.
- The `fast_filter` pre-check (on by default, `--no-fast-filter` to
  disable) skips tokenizing a logical line that cannot possibly be a
  dependency-relevant statement; it is provably never stricter than
  the full parser (see TestFilterNeverStricterThanParser) and kept
  removable as a permanent escape hatch, not just trusted.

- A module that is USEd but defined by none of the scanned files is
  reported as a NOTE (see find_unresolved_uses). Usually that is a
  genuinely external module -- xios, netcdf, mpi, AGRIF -- so it is
  not a warning and changes nothing in the output; it is reported so
  that the other cause, renaming the file that defines a module,
  does not silently drop a dependency edge.

Known limitations:
- A USE target that is not a plain identifier (e.g. a macro-
  generated module name) is warned about and not recorded as a
  dependency -- not currently observed anywhere in this tree, but
  worth checking if this warning ever fires for real.
- Only the first statement on a logical line is parsed, so a
  semicolon-separated line like "X; USE modb" would miss the second
  statement -- not currently observed anywhere in this tree.
- Intra-target module build ordering (CMake dyndep-style: building
  A's .mod before compiling B, both in one target, within a single
  parallel `make -j` pass) is not implemented.
"""
from __future__ import annotations
import sys
from collections import deque
from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from typing import Iterator, List, Optional


@dataclass
class LogicalLine:
    text: str
    first_lineno: int


def _is_fixed_form(path: Path) -> bool:
    return path.suffix.lower() == ".f"


def _strip_fixed_form_comment(line: str) -> Optional[str]:
    if not line or line[0] in "cC*!":
        return None
    return line


def _strip_bang_comment(line: str) -> str:
    # quote-aware truncation at an unquoted "!" -- used for both
    # free-form lines (its only comment marker) and, in addition to
    # the column-1 whole-line check, fixed-form lines, since real
    # .F files in this tree also use trailing "! comment" (e.g.
    # main.F, step.F), not just the classic column-1 convention.
    out: List[str] = []
    quote: Optional[str] = None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            out.append(ch)
            continue
        if ch == "!":
            break
        out.append(ch)
    return "".join(out)


def _strip_block_comment(line: str, in_block: bool):
    # C-style /* ... */ comments, common in this tree's CPP-heavy
    # headers, quote-aware and can span multiple lines (`in_block`
    # is threaded across calls by the caller). Comment content is
    # replaced with spaces, not removed, so column positions used by
    # the fixed-form column-6-continuation check stay meaningful.
    out: List[str] = []
    quote: Optional[str] = None
    i, n = 0, len(line)
    while i < n:
        ch = line[i]
        if in_block:
            if ch == "*" and i + 1 < n and line[i + 1] == "/":
                in_block = False
                out.append("  ")
                i += 2
                continue
            out.append(" ")
            i += 1
            continue
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
            i += 1
            continue
        if ch in "'\"":
            quote = ch
            out.append(ch)
            i += 1
            continue
        if ch == "/" and i + 1 < n and line[i + 1] == "*":
            in_block = True
            out.append("  ")
            i += 2
            continue
        out.append(ch)
        i += 1
    return "".join(out), in_block


def _split_fixed_form(line: str) -> tuple[bool, str]:
    # Split a fixed-form line into (is_continuation, statement text).
    #
    # Besides the standard column-6 convention, this handles the
    # GNU/DEC "tab form" extension, which cross_matrix.F also accepted
    # (it tested explicitly for char(9) in the label field): a TAB
    # anywhere in columns 1-6 terminates the label field, so the
    # statement starts immediately after it -- or, if a nonzero digit
    # follows the TAB, that digit is the continuation marker and the
    # statement starts after the digit. Without this, a tab-indented
    # statement is sliced at the wrong offset and its first five
    # characters are silently eaten, which loses the dependency for a
    # tab-indented USE/INCLUDE.
    #
    # Both boundaries here were checked against gfortran itself:
    #   "\t1 5"       after "x = 1 +"  -> continuation, statement " 5"
    #   "     \t2"    after "x = 1 +"  -> a new statement "2"
    #   "\tuse mymod"                  -> a statement "use mymod"
    # so a TAB in column 6 ends the label field but is too late to be
    # followed by a continuation marker: only a TAB in columns 1-5,
    # the label field proper, can be.
    tab = line.find("\t")
    if 0 <= tab < 6:
        rest = line[tab + 1:]
        if tab < 5 and rest[:1].isdigit() and rest[0] != "0":
            return True, rest[1:]
        return False, rest
    is_cont = len(line) > 5 and line[5] not in " 0"
    return is_cont, line[6:] if len(line) > 6 else ""


def iter_logical_lines(path: Path) -> Iterator[LogicalLine]:
    fixed = _is_fixed_form(path)
    raw_lines = path.read_text(encoding="latin-1").splitlines()

    buf: List[str] = []
    start_lineno = 0
    in_block_comment = False
    for lineno, raw_line in enumerate(raw_lines, start=1):
        line, in_block_comment = _strip_block_comment(raw_line, in_block_comment)
        if fixed:
            kept = _strip_fixed_form_comment(line)
            if kept is None:
                continue
            kept = _strip_bang_comment(kept)
            if kept[:1] == "#":
                # a CPP directive: column 1 is the '#' itself, not a
                # label field, so it does not follow the fixed-form
                # column-6-continuation convention at all
                is_cont = False
                payload = kept
            else:
                is_cont, payload = _split_fixed_form(kept)
        else:
            kept = _strip_bang_comment(line)
            if not kept.strip():
                continue
            is_cont = bool(buf) and buf[-1].rstrip().endswith("&")
            if is_cont:
                buf[-1] = buf[-1].rstrip()[:-1]
            payload = kept.lstrip().lstrip("&") if is_cont else kept

        if is_cont and buf:
            buf.append(payload)
        else:
            if buf:
                yield LogicalLine("".join(buf), start_lineno)
            buf = [payload]
            start_lineno = lineno

    if buf:
        yield LogicalLine("".join(buf), start_lineno)


class TokKind(Enum):
    IDENT = auto()
    STRING = auto()
    COMMA = auto()
    COLON = auto()
    COLONCOLON = auto()
    HASH = auto()
    OTHER = auto()
    EOL = auto()


@dataclass
class Token:
    kind: TokKind
    value: str
    col: int


def tokenize(text: str) -> List[Token]:
    tokens: List[Token] = []
    i, n = 0, len(text)
    while i < n:
        ch = text[i]
        if ch in " \t":
            i += 1
            continue
        if ch in "'\"":
            quote, start = ch, i
            i += 1
            buf: List[str] = []
            while i < n:
                if text[i] == quote:
                    if i + 1 < n and text[i + 1] == quote:
                        buf.append(quote)
                        i += 2
                        continue
                    i += 1
                    break
                buf.append(text[i])
                i += 1
            tokens.append(Token(TokKind.STRING, "".join(buf), start))
            continue
        if ch.isalpha() or ch == "_":
            start = i
            while i < n and (text[i].isalnum() or text[i] == "_"):
                i += 1
            tokens.append(Token(TokKind.IDENT, text[start:i], start))
            continue
        if ch == ":" and i + 1 < n and text[i + 1] == ":":
            tokens.append(Token(TokKind.COLONCOLON, "::", i))
            i += 2
            continue
        if ch == ":":
            tokens.append(Token(TokKind.COLON, ":", i))
            i += 1
            continue
        if ch == ",":
            tokens.append(Token(TokKind.COMMA, ",", i))
            i += 1
            continue
        if ch == "#":
            tokens.append(Token(TokKind.HASH, "#", i))
            i += 1
            continue
        tokens.append(Token(TokKind.OTHER, ch, i))
        i += 1
    tokens.append(Token(TokKind.EOL, "", n))
    return tokens


_RELEVANT_KEYWORDS = {"use", "include", "module"}


def looks_relevant(text: str) -> bool:
    s = text.lstrip()
    if not s:
        return False
    if s[0] == "#":
        return True
    i = 0
    n = len(s)
    while i < n and (s[i].isalnum() or s[i] == "_"):
        i += 1
    return s[:i].lower() in _RELEVANT_KEYWORDS


@dataclass
class UseStmt:
    module: str


@dataclass
class IncludeStmt:
    filename: str


@dataclass
class ModuleDefStmt:
    name: str


@dataclass
class UnresolvedUse:
    raw: str


class TokenStream:
    def __init__(self, tokens: List[Token]):
        self._tokens = tokens
        self._pos = 0

    def peek(self) -> Token:
        return self._tokens[self._pos]

    def advance(self) -> Token:
        tok = self._tokens[self._pos]
        self._pos += 1
        return tok

    def at_keyword(self, word: str) -> bool:
        tok = self.peek()
        return tok.kind is TokKind.IDENT and tok.value.lower() == word

    def eat_keyword(self, word: str) -> bool:
        if self.at_keyword(word):
            self.advance()
            return True
        return False


_MODULE_NON_DEF = {"procedure", "subroutine", "function", "pure", "elemental"}


def parse_statement(tokens: List[Token]):
    ts = TokenStream(tokens)
    if ts.peek().kind is TokKind.HASH:
        ts.advance()
        if not ts.eat_keyword("include"):
            return None
        return _parse_include(ts)
    if ts.eat_keyword("include"):
        return _parse_include(ts)
    if ts.eat_keyword("module"):
        return _parse_module(ts)
    if ts.eat_keyword("use"):
        return _parse_use(ts)
    return None


def _parse_include(ts: TokenStream):
    tok = ts.peek()
    if tok.kind is TokKind.STRING:
        ts.advance()
        return IncludeStmt(tok.value)
    return None


def _parse_module(ts: TokenStream):
    tok = ts.peek()
    if tok.kind is TokKind.IDENT and tok.value.lower() in _MODULE_NON_DEF:
        return None
    if tok.kind is TokKind.IDENT:
        ts.advance()
        if ts.peek().kind is TokKind.EOL:
            return ModuleDefStmt(tok.value.lower())
    return None


def _parse_use(ts: TokenStream):
    if ts.peek().kind is TokKind.COMMA:
        ts.advance()
        if ts.at_keyword("intrinsic") or ts.at_keyword("non_intrinsic"):
            ts.advance()
    if ts.peek().kind is TokKind.COLONCOLON:
        ts.advance()
    tok = ts.peek()
    if tok.kind is TokKind.IDENT:
        ts.advance()
        if ts.peek().kind in (TokKind.EOL, TokKind.COMMA):
            return UseStmt(tok.value.lower())
    return UnresolvedUse(tok.value)


@dataclass
class FileInfo:
    path: Path
    includes: set[str]
    uses: set[str]
    defines: set[str]
    is_given: bool
    found: bool = True


@dataclass
class ScanResult:
    includes: List[tuple[str, int]]
    uses: List[tuple[str, int]]
    defines: List[tuple[str, int]]
    unresolved_uses: List[tuple[str, int]]


def scan_file(path: Path, fast_filter: bool = True) -> ScanResult:
    result = ScanResult([], [], [], [])
    for ll in iter_logical_lines(path):
        if fast_filter and not looks_relevant(ll.text):
            continue
        stmt = parse_statement(tokenize(ll.text))
        if isinstance(stmt, IncludeStmt):
            result.includes.append((stmt.filename, ll.first_lineno))
        elif isinstance(stmt, ModuleDefStmt):
            result.defines.append((stmt.name, ll.first_lineno))
        elif isinstance(stmt, UseStmt):
            result.uses.append((stmt.module, ll.first_lineno))
        elif isinstance(stmt, UnresolvedUse):
            result.unresolved_uses.append((stmt.raw, ll.first_lineno))
    return result


def build_catalog(
    given_paths: List[Path], fast_filter: bool = True
) -> tuple[dict[str, FileInfo], List[str]]:
    catalog: dict[str, FileInfo] = {}
    warnings: List[str] = []
    given_names = {p.name for p in given_paths}
    queue: deque[Path] = deque(given_paths)

    while queue:
        path = queue.popleft()
        if path.name in catalog:
            continue
        if not path.exists():
            warnings.append(f"WARNING: file not found, excluded: {path.name}")
            catalog[path.name] = FileInfo(
                path, set(), set(), set(), path.name in given_names, found=False
            )
            continue

        result = scan_file(path, fast_filter=fast_filter)
        for name, lineno in result.unresolved_uses:
            warnings.append(
                f"WARNING: {path.name}:{lineno}: USE target {name!r} is not "
                "a plain identifier (macro-generated module name?); "
                "dependency NOT recorded -- check manually"
            )
        catalog[path.name] = FileInfo(
            path=path,
            includes={n for n, _ in result.includes},
            uses={n for n, _ in result.uses},
            defines={n for n, _ in result.defines},
            is_given=path.name in given_names,
        )
        for inc_name, _ in result.includes:
            if inc_name not in catalog:
                queue.append(path.parent / inc_name)

    return catalog, warnings


def build_module_registry(
    catalog: dict[str, FileInfo]
) -> tuple[dict[str, str], List[str]]:
    registry: dict[str, str] = {}
    warnings: List[str] = []
    for name, info in catalog.items():
        for mod in info.defines:
            if mod in registry and registry[mod] != name:
                warnings.append(
                    f"WARNING: module {mod!r} defined in both "
                    f"{registry[mod]!r} and {name!r}; keeping {registry[mod]!r}"
                )
                continue
            registry[mod] = name
    return registry, warnings


def _object_name(filename: str) -> str:
    # matches cross_matrix.F's own scope: CROCO sources are always
    # uppercase .F/.F90 (the preprocessed .f/.f90 are intermediate
    # build products, never scanned); anything else is returned
    # unchanged rather than guessed at
    if filename.endswith(".F90"):
        return filename[:-4] + ".o"
    if filename.endswith(".F"):
        return filename[:-2] + ".o"
    return filename


def _header_closure(catalog: dict[str, FileInfo], start: str) -> set[str]:
    seen: set[str] = set()
    stack = [start]
    while stack:
        cur = stack.pop()
        info = catalog.get(cur)
        if info is None:
            continue
        for inc in info.includes:
            # exclude `start` itself: a cycle back to the file
            # being resolved must not make it its own prerequisite
            if inc != start and inc not in seen:
                seen.add(inc)
                stack.append(inc)
    return seen


def _used_modules(
    catalog: dict[str, FileInfo], name: str, headers: set[str]
) -> set[str]:
    # a USE inside an #include'd header is a USE by the including file,
    # since the header is spliced in before the compiler sees it
    mods = set(catalog[name].uses)
    for h in headers:
        mods |= catalog[h].uses
    return mods


def find_unresolved_uses(
    catalog: dict[str, FileInfo], registry: dict[str, str]
) -> List[str]:
    """
    Report modules that are USEd but defined by none of the scanned
    files, as NOTE lines.

    Most are legitimately external -- xios, netcdf, mpi, the AGRIF
    library -- so this is a NOTE, not a WARNING, and it does not
    change the emitted dependencies. It is reported at all because the
    alternative is silence, and silence is exactly how renaming the
    file that defines a module drops its dependency edge without a
    trace. cross_matrix.F did at least print something here, as a side
    effect of guessing a filename and then failing to open it.
    """
    referrers: dict[str, set[str]] = {}
    for name, info in catalog.items():
        if not info.is_given:
            continue
        headers = _header_closure(catalog, name)
        for mod in _used_modules(catalog, name, headers):
            if mod not in registry:
                referrers.setdefault(mod, set()).add(name)

    notes: List[str] = []
    for mod, files in sorted(referrers.items()):
        shown = sorted(files)
        suffix = ", ..." if len(shown) > 3 else ""
        notes.append(
            f"NOTE: module {mod!r} is USEd by {len(shown)} scanned file(s) "
            f"({', '.join(shown[:3])}{suffix}) but is defined by none of "
            "them; no dependency recorded (external module?)"
        )
    return notes


def resolve_dependencies(
    catalog: dict[str, FileInfo], registry: dict[str, str]
) -> dict[str, set[str]]:
    deps: dict[str, set[str]] = {}

    for name, info in catalog.items():
        if not info.is_given:
            continue
        headers = _header_closure(catalog, name)
        used_modules = _used_modules(catalog, name, headers)
        # convert module dependencies to their .o form here, where the
        # distinction between "reached via #include" (headers, kept
        # as literal source filenames -- e.g. a CROCO idiom like
        # biology.F #include-ing bio_BioEBUS.F as raw spliced-in
        # source, not a separately compiled+linked unit) and "reached
        # via USE" (a real .o build target) is still known; the
        # writer must not have to guess this back from a filename
        # suffix alone, since both can end in .F/.F90
        module_files = {
            _object_name(registry[m]) for m in used_modules
            if m in registry and registry[m] != name
        }
        # a header that was never found on disk (see build_catalog's
        # EXCLUDE handling) must not be emitted as a literal
        # prerequisite: make has no rule to satisfy it and would
        # fail outright ("No rule to make target ..."), whereas
        # simply not tracking it only means a change to that file
        # (which was never found in the first place) cannot trigger
        # a rebuild -- matching cross_matrix.F's own EXCLUDE
        # philosophy of dropping such files from the dependency list
        found_headers = {h for h in headers if catalog[h].found}
        deps[name] = found_headers | module_files

    return deps


def write_make_depend(out_path: Path, given_paths: List[Path],
                       deps: dict[str, set[str]], warnings: List[str]) -> None:
    with out_path.open("w") as fh:
        fh.write(
            "# Make.depend: list of dependencies generated by\n"
            "# make_depend.py.\n"
            "# WARNING: THIS IS A MACHINE GENERATED FILE: DO NOT EDIT !!!!!\n"
            "# To create or update this file use:\n"
            "#        ./make_depend.py *.F90 *.F\n"
            "# or      make depend\n"
        )
        for w in warnings:
            fh.write(f"# {w}\n")
        for path in given_paths:
            name = path.name
            target = _object_name(name)
            prereqs = sorted(deps.get(name, set()))
            fh.write(f"#\n{target}: {name} " + " ".join(prereqs) + "\n")


_USAGE = """\
usage: make_depend.py [-h] [-o OUTPUT] [--no-fast-filter] file.F90 file.F ...

Scan Fortran/CPP source files and write a make dependency file,
replacing cross_matrix.F (see the module docstring for design notes).

  -h, --help         show this message and exit
  -o, --output PATH  write the dependency file to PATH
                     (default: Make.depend in the current directory)
  --no-fast-filter   disable the statement pre-filter and tokenize
                     every logical line unconditionally (see the
                     module docstring); the two must always agree,
                     this is an escape hatch, not a normal option
"""


def _take_output_option(argv: List[str]) -> tuple[Optional[str], List[str]]:
    # minimal hand-rolled option parsing, to keep the CLI stdlib-only
    # and dependency-free in the same spirit as the rest of the tool
    out: Optional[str] = None
    rest: List[str] = []
    i = 0
    while i < len(argv):
        arg = argv[i]
        if arg in ("-o", "--output"):
            if i + 1 >= len(argv):
                raise ValueError(f"{arg} requires an argument")
            out = argv[i + 1]
            i += 2
            continue
        if arg.startswith("--output="):
            out = arg[len("--output="):]
            i += 1
            continue
        rest.append(arg)
        i += 1
    return out, rest


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if "-h" in argv or "--help" in argv:
        print(_USAGE, end="")
        return 0
    fast_filter = "--no-fast-filter" not in argv
    argv = [a for a in argv if a != "--no-fast-filter"]
    try:
        output, argv = _take_output_option(argv)
    except ValueError as exc:
        print(f"ERROR IN make_depend: {exc}", file=sys.stderr)
        return 1
    if not argv:
        print("ERROR IN make_depend: NO FILE NAMES ARE GIVEN", file=sys.stderr)
        return 1
    # a Makefile pattern rule can pass a stray *.mod glob alongside the
    # real sources (mirrors cross_matrix.F's own handling of this);
    # .mod files are compiler-generated module files, never a source
    # to scan
    given_paths = [Path(a) for a in argv if not a.endswith(".mod")]
    catalog, cat_warnings = build_catalog(given_paths, fast_filter=fast_filter)
    registry, reg_warnings = build_module_registry(catalog)
    warnings = cat_warnings + reg_warnings
    notes = find_unresolved_uses(catalog, registry)
    deps = resolve_dependencies(catalog, registry)
    write_make_depend(
        Path(output or "Make.depend"), given_paths, deps, warnings + notes
    )
    for line in warnings + notes:
        print(line, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
