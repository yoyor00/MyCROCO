#!/usr/bin/env python3
"""
Differential check of make_depend.py against the legacy cross_matrix.F,
run against the real OCEAN/*.F90 *.F tree.

Both tools run in a scratch directory of symlinks, so neither one can
clobber a real Make.depend in the source tree. The check is that
make_depend.py's dependency set for every shared target is either a
superset of cross_matrix.F's, or that any "missing" .o-suffixed
dependency is still reachable through make_depend.py's own .o
dependency chain (the documented design: module dependencies rely on
make's own transitivity, only #include chains need an explicit
closure -- see the make_depend.py module docstring).

Why the baseline is patched
---------------------------
cross_matrix.F is compiled here with max_name_size raised from 62 to
300, matching the declared length of `character*300 fname`. As
shipped, it blanks only the first 62 columns of a 300-character
filename buffer and then passes the whole buffer to open(), so under
gfortran every file discovered through #include fails to open and is
dropped, leaving almost every target with no header dependencies at
all. Comparing against that baseline would be vacuous -- almost any
output is a superset of nearly nothing. Raising the parameter repairs
the buffer without touching any of the logic being compared, and
gives a baseline that carries the header dependencies the original
algorithm was meant to produce.

Measured on 2026-10-13 over the 254 sources in OCEAN/: unpatched, the
baseline emits 13 header prerequisites in total and leaves 249 of its
254 targets without a single one; patched, it emits 2867, against
2813 from make_depend.py. Re-measure rather than trust these figures
if the tree has moved on.

Earlier revisions of this docstring quoted 108 of 108 discovered
files and 342 of 362 targets. Those numbers are real, but they were
measured on the flat Compile/ tree jobcomp assembles -- which also
carries PISCES/, MUSTANG/, STOGEN/ and XIOS/ -- not on the OCEAN/
-only scope this check actually runs in. Keep the two apart when
quoting either.

Note that the repaired baseline still silently drops one prerequisite
every time an output line reaches max_string_size, and still resolves
a USE by guessing that the module name is the filename. Those are
real cross_matrix.F defects, not make_depend.py regressions, which is
why this check is one-directional: new >= old, never old >= new.

Scope: this runs in OCEAN/ alone, not in the flat tree jobcomp
assembles in Compile/, so modules living in PISCES/, MUSTANG/ and
STOGEN/ are simply absent here and unresolvable for both tools. That
is fine for a regression gate, and it keeps the check from depending
on jobcomp's copy list.

Needs gfortran. Run it manually, or as a CI step:

    python3 OCEAN/tests/check_differential.py

Exits 0 if there are no unexplained differences, 1 otherwise.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
OCEAN = REPO_ROOT / "OCEAN"

# see "Why the baseline is patched" above
_BASELINE_PATCH = ("max_name_size=62", "max_name_size=300")


def parse_make_depend(path: Path) -> dict[str, set[str]]:
    deps: dict[str, set[str]] = defaultdict(set)
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        m = re.match(r"^(\S+\.o)\s*:\s*(.*)$", line)
        if not m:
            continue
        deps[m.group(1)] |= set(m.group(2).split())
    return deps


def transitive_closure(start: str, deps: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    stack = [start]
    while stack:
        cur = stack.pop()
        for d in deps.get(cur, set()):
            if d not in seen:
                seen.add(d)
                stack.append(d)
    return seen


def patched_baseline_source(dest: Path) -> None:
    text = (OCEAN / "cross_matrix.F").read_text()
    old, new = _BASELINE_PATCH
    if text.count(old) != 1:
        raise SystemExit(
            f"cross_matrix.F no longer contains exactly one {old!r}; "
            "the baseline patch in check_differential.py needs revisiting"
        )
    dest.write_text(text.replace(old, new))


def main() -> int:
    if shutil.which("gfortran") is None:
        print("gfortran not found; cannot build the cross_matrix.F baseline",
              file=sys.stderr)
        return 1

    sources = sorted(OCEAN.glob("*.F90")) + sorted(OCEAN.glob("*.F"))
    if not sources:
        print(f"no *.F90/*.F found under {OCEAN}", file=sys.stderr)
        return 1

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        # a scratch tree of symlinks: cross_matrix.F always writes
        # Make.depend into its own working directory, so running it
        # in OCEAN/ would clobber a real one
        for pattern in ("*.F90", "*.F", "*.h"):
            for src in OCEAN.glob(pattern):
                (work / src.name).symlink_to(src)

        # the baseline is built from a separate, differently-named
        # copy: cross_matrix.F itself stays in the analysed source
        # list, exactly as it is when `make depend` globs *.F
        baseline_src = work / "cross_matrix_baseline.F"
        patched_baseline_source(baseline_src)
        baseline = work / "cross_matrix_baseline"
        subprocess.run(
            ["gfortran", "-O2", "-o", str(baseline), str(baseline_src)],
            check=True,
        )

        names = [p.name for p in sources]
        old_out = work / "Make.depend"
        subprocess.run([str(baseline), *names], cwd=work, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        old_deps = parse_make_depend(old_out)
        old_out.unlink()

        new_out = work / "new.depend"
        subprocess.run(
            [sys.executable, str(OCEAN / "make_depend.py"),
             "-o", str(new_out), *names],
            cwd=work, check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        new_deps = parse_make_depend(new_out)

    only_old = set(old_deps) - set(new_deps)
    only_new = set(new_deps) - set(old_deps)
    unexplained: list[tuple[str, str]] = []
    for target in sorted(set(old_deps) & set(new_deps)):
        missing = old_deps[target] - new_deps[target]
        if not missing:
            continue
        closure = transitive_closure(target, new_deps)
        for dep in missing:
            if dep not in closure:
                unexplained.append((target, dep))

    old_edges = sum(len(v) for v in old_deps.values())
    new_edges = sum(len(v) for v in new_deps.values())
    print(f"targets: old={len(old_deps)} new={len(new_deps)}")
    print(f"prerequisites: old={old_edges} new={new_edges}")
    print(f"targets only in old: {len(only_old)}")
    print(f"targets only in new: {len(only_new)}")
    print(f"unexplained missing dependencies: {len(unexplained)}")
    for target, dep in unexplained:
        print(f"  {target}: missing {dep!r} (not reachable via new tool's own .o graph)")

    ok = not only_old and not only_new and not unexplained
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
