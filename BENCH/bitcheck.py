#!/usr/bin/env python3
"""Independent bit-level comparison of CROCO netCDF outputs.

Compares a reference file against one or more candidate files, variable by
variable, in float64 after casting the stored values, and reports:
  - number of differing elements
  - max abs/rel diff
  - first record (time) index where a difference shows up
  - the sha256 of the raw variable buffer (true bitwise check)
"""
import sys
import hashlib
import numpy as np
from netCDF4 import Dataset

SKIP = {"spherical", "hc", "time_step"}


def var_bytes(path, name):
    with Dataset(path, "r") as ds:
        v = ds.variables[name][...]
    a = np.ascontiguousarray(v)
    return hashlib.sha256(a.tobytes()).hexdigest(), a


def compare(ref, cand):
    print(f"\n### ref = {ref}\n### cur = {cand}")
    with Dataset(ref) as dr, Dataset(cand) as dc:
        names = [n for n in dr.variables if n not in SKIP]
        all_ok = True
        for n in names:
            if n not in dc.variables:
                print(f"  {n}: MISSING in candidate")
                all_ok = False
                continue
            a = dr.variables[n][...].astype(np.float64)
            b = dc.variables[n][...].astype(np.float64)
            if a.shape != b.shape:
                print(f"  {n}: SHAPE DIFF {a.shape} vs {b.shape}")
                all_ok = False
                continue
            if a.size == 0:
                continue
            ha = hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
            hb = hashlib.sha256(np.ascontiguousarray(b).tobytes()).hexdigest()
            if ha == hb:
                print(f"  {n:>12}: BITWISE IDENTICAL ({a.size} values)")
                continue
            all_ok = False
            d = np.abs(a - b)
            nz = d > 0
            ndiff = int(nz.sum())
            maxabs = float(d.max())
            scale = np.maximum(np.abs(a), np.abs(b))
            rel = np.where(scale > 0, d / np.maximum(scale, 1e-300), 0.0)
            maxrel = float(rel.max())
            # first time index (if any dim can be time)
            first_t = None
            if a.ndim >= 2:
                per_t = nz.reshape(a.shape[0], -1).any(axis=1)
                idx = np.nonzero(per_t)[0]
                if idx.size:
                    first_t = int(idx[0])
            print(
                f"  {n:>12}: DIFFERS  n={ndiff}/{a.size} "
                f"({100.0*ndiff/a.size:.3f}%)  max_abs={maxabs:.3e}  "
                f"max_rel={maxrel:.3e}  first_rec={first_t}"
            )
        return all_ok


if __name__ == "__main__":
    ref = sys.argv[1]
    ok = True
    for cand in sys.argv[2:]:
        ok = compare(ref, cand) and ok
    print("\nRESULT:", "ALL BITWISE IDENTICAL" if ok else "DIFFERENCES FOUND")
