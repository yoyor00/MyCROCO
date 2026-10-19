#!/usr/bin/env python3
"""
gen_cppdefs_info.py -- Generate the cppdefs_info Fortran 90 module from the
active CPP keys, resolved by the preprocessor at build time.

Replaces cppcheck.F + check_switches1.F.

Pipeline:
  1. Run `cpp -dM` on cppdefs.h with the same CPP/CPPFLAGS as the build.
     Yields the full set of macros active after preprocessing.
  2. Keep value-less macros that CROCO itself knows: identifiers found in
     cppdefs.h / cppdefs_dev.h (#ifdef, #ifndef, defined(), #define).
     System macros (__GNUC__, ...) are filtered out.
  3. Collect build information:
       - code version: `git describe --tags --always --dirty` in --srcdir,
         'unknown' if there is no git or no repository;
       - local files: source files overridden from the run directory
         by jobcomp (list given by --local-files), with a short md5;
         'unknown' if the list is not available, 'none' if empty.
  4. Emit cppdefs_info.F90, a free-form module with:
       - subroutine print_cppdefs  : fills Coptions (NetCDF global
         attribute, see put_global_atts) and prints the active keys
         on stdout
       - code_version, local_files : written as NetCDF global attributes
"""

import argparse
import hashlib
import os
import re
import shlex
import subprocess
import sys
from datetime import datetime

RE_IFDEF = re.compile(r"#\s*if(?:n?def)\s+([A-Za-z_]\w*)")
RE_DEFINED = re.compile(r"defined\s*\(?\s*([A-Za-z_]\w*)")
RE_DEFINE = re.compile(r"#\s*define\s+([A-Za-z_]\w*)")


def read_whitelist(headers):
    """Extract all CPP identifiers mentioned in the given header files."""
    whitelist = set()
    for path in headers:
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                text = f.read()
        except FileNotFoundError:
            print(
                f"[gen_cppdefs_info] warning: {path} not found, skipped",
                file=sys.stderr,
            )
            continue
        for rx in (RE_IFDEF, RE_DEFINED, RE_DEFINE):
            whitelist.update(rx.findall(text))
    whitelist = {k for k in whitelist if not (k.startswith("__") and k.endswith("__"))}
    return whitelist


def active_flags(cpp, cppflags, includes, cppdefs_h):
    """Run `cpp -dM` and return the set of value-less defined macro names."""
    cmd = shlex.split(cpp) + ["-dM", "-E"]
    if cppflags:
        cmd.extend(cppflags.split())
    cmd.extend(includes)
    cmd.append(cppdefs_h)
    try:
        out = subprocess.check_output(cmd, text=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print(f"[gen_cppdefs_info] cpp -dM failed:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)

    flags = set()
    for line in out.splitlines():
        # Only value-less "#define FOO" lines are CROCO toggles.
        m = re.match(r"#\s*define\s+([A-Za-z_]\w*)\s*$", line)
        if m:
            flags.add(m.group(1))
    return flags


def get_code_version(srcdir):
    """`git describe --tags --always --dirty` in srcdir, or 'unknown'."""
    if not srcdir:
        return "unknown"
    try:
        out = subprocess.check_output(
            ["git", "-C", srcdir, "describe", "--tags", "--always", "--dirty"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return out or "unknown"


def _md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def _same_as_source(name, digest, srcdir):
    """True if a file with the same name and content exists in the sources."""
    if not srcdir:
        return False
    root = os.path.dirname(os.path.abspath(srcdir))  # CROCO root
    for dirpath, _, files in os.walk(root):
        if name in files:
            try:
                if _md5(os.path.join(dirpath, name)) == digest:
                    return True
            except OSError:
                pass
    return False


def get_local_files(listfile, srcdir=""):
    """
    'name(md5) ...' for the files copied from the run directory by jobcomp
    that differ from the sources (or do not exist there).
    cppdefs.h and param.h always differ (they define the configuration)
    and are not listed.
    """
    if not listfile or not os.path.isfile(listfile):
        return "unknown"
    with open(listfile) as f:
        names = [n for n in f.read().split() if n not in ("cppdefs.h", "param.h")]
    items = []
    for name in sorted(set(names)):
        try:
            digest = _md5(name)
        except OSError:
            items.append(f"{name}(?)")
            continue
        if not _same_as_source(name, digest, srcdir):
            items.append(f"{name}({digest[:8]})")
    return " ".join(items) if items else "none"


# ---------------------------------------------------------------------------
# Fortran 90 module template
# ---------------------------------------------------------------------------

MODULE_TEMPLATE = """\
#include "cppdefs.h"
!=======================================================================
!  MODULE cppdefs_info  --  auto-generated, DO NOT EDIT
!
!  Regenerated at build time by gen_cppdefs_info.py from the active CPP
!  keys resolved on cppdefs.h (+ cppdefs_dev.h) by the preprocessor.
!
!  Build stamp   : {stamp}
!  Code version  : {code_version}
!  Active keys   : {n}
!=======================================================================

MODULE cppdefs_info
   implicit none
   private

   public :: print_cppdefs
   public :: code_version, local_files

   integer, parameter :: n_active_keys = {n}
   integer, parameter :: key_len       = {key_len}

   ! git describe of the sources, 'unknown' if not available
   character(len={cv_len}), parameter :: code_version = '{code_version}'
   ! Files overridden from the run directory (filled by print_cppdefs)
   character(len={lf_len}) :: local_files = ' '

contains

   !---------------------------------------------------------------------
   !  print_cppdefs
   !  Fill Coptions (written by put_global_atts as a NetCDF global
   !  attribute) and print the active keys on the master rank.
   !---------------------------------------------------------------------
   subroutine print_cppdefs(ierr)
      use param, ONLY: stdout
#if defined MPI
      use scalars, ONLY: mynode
#endif
      implicit none
      integer, intent(inout) :: ierr
      integer :: i, is, ie
      character(len=key_len) :: active_keys(n_active_keys)
#include "strings.h"

      ! One assignment per key: a single long array constructor makes
      ! the AGRIF conv tool crash (statement too long).
{keys_list}
{local_files_fill}

      Coptions = ' '
      is = 1
      do i = 1, n_active_keys
         ie = is + len_trim(active_keys(i)) - 1
         if (ie >= max_opt_size) then
            MPI_master_only write (stdout, '(/1x,A,A/14x,A)') &
               'CPPDEFS_INFO -- ERROR: insufficient size of string Coptions ', &
               'in file "strings.h".', 'Increase the size and recompile.'
            ierr = ierr + 1
            return
         end if
         Coptions(is:ie) = trim(active_keys(i))
         is = ie + 2
      end do

      MPI_master_only write (stdout, '(/1x,A/)') 'Activated C-preprocessing Options:'
      do i = 1, n_active_keys
         MPI_master_only write (stdout, '(10x,A)') trim(active_keys(i))
      end do
      MPI_master_only write (stdout, '(/1x,A,1x,A)') 'Code version:', trim(code_version)
      MPI_master_only write (stdout, '(1x,A,1x,A/)') 'Local files :', trim(local_files)
   end subroutine print_cppdefs

END MODULE cppdefs_info
"""


def format_string_fill(var, text, width=50):
    """Assign text to var by short substrings (keeps statements short)."""
    return "\n".join(
        f"      {var}({i + 1}:{i + len(text[i:i + width])}) = " f"'{text[i:i + width]}'"
        for i in range(0, len(text), width)
    )


def format_keys_list(keys):
    return "\n".join(
        f"      active_keys({i}) = '{k}'" for i, k in enumerate(keys, start=1)
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--cppdefs", default="cppdefs.h")
    ap.add_argument("--dev", default="cppdefs_dev.h")
    ap.add_argument("--cpp", default="cpp")
    ap.add_argument("--cppflags", default="")
    ap.add_argument("-I", dest="includes", action="append", default=["-I."])
    ap.add_argument("-o", dest="out", default="cppdefs_info.F90")
    ap.add_argument(
        "--srcdir", default="", help="CROCO source directory (for git describe)"
    )
    ap.add_argument(
        "--local-files",
        default="",
        help="file listing the sources overridden by jobcomp",
    )
    args = ap.parse_args()

    whitelist = read_whitelist([args.cppdefs, args.dev])
    flags = active_flags(args.cpp, args.cppflags, args.includes, args.cppdefs)
    keys = sorted(flags & whitelist)

    if not keys:
        print("[gen_cppdefs_info] ERROR: no active CPP keys found.", file=sys.stderr)
        sys.exit(1)

    key_len = max(len(k) for k in keys)
    keys_list = format_keys_list(keys)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    version = get_code_version(args.srcdir)
    locals_ = get_local_files(args.local_files, args.srcdir)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(
            MODULE_TEMPLATE.format(
                stamp=stamp,
                n=len(keys),
                key_len=key_len,
                keys_list=keys_list,
                code_version=version,
                cv_len=len(version),
                lf_len=len(locals_),
                local_files_fill=format_string_fill("local_files", locals_),
            )
        )

    print(
        f"[gen_cppdefs_info] {len(keys)} active CPP keys, "
        f"version={version} -> {args.out}"
    )


if __name__ == "__main__":
    main()
