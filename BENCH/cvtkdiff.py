#!/usr/bin/env python3
"""Compare CROCO CVTK_DEBUG_WRITE dumps between a sequential run (1 rank) and
a MPI run (N ranks) to find WHERE (which checkpoint, which point) the two
start to differ -- including in the ghost cells vs the interior.

CROCO writes, for each `check_tab2d/3d` call, four Fortran unformatted records:
    i1,i2,j1,j2          (global index range covered by the check)
    lbound(tab)[+offset] (global lower bound of the dumped array)
    ubound(tab)[+offset] (global upper bound of the dumped array)
    tab                  (the whole array, ghost cells included)
"""
import sys
import struct
import numpy as np

LM = 60  # LLm of the BASIN case
MM = 50  # MMm


def records(path):
    """Yield the raw payload of each Fortran sequential-unformatted record."""
    with open(path, "rb") as f:
        data = f.read()
    pos = 0
    n = len(data)
    while pos < n:
        (ln,) = struct.unpack_from("<i", data, pos)
        pos += 4
        payload = data[pos:pos + ln]
        pos += ln
        (ln2,) = struct.unpack_from("<i", data, pos)
        pos += 4
        assert ln == ln2, (ln, ln2)
        yield payload


def checkpoints(path):
    """Parse the dump file into a list of (i1,i2,j1,j2,lb,ub,array)."""
    out = []
    recs = list(records(path))
    assert len(recs) % 4 == 0, len(recs)
    for k in range(0, len(recs), 4):
        i1, i2, j1, j2 = struct.unpack("<4i", recs[k])
        lb = struct.unpack("<%di" % (len(recs[k + 1]) // 4), recs[k + 1])
        ub = struct.unpack("<%di" % (len(recs[k + 2]) // 4), recs[k + 2])
        dims = [ub[d] - lb[d] + 1 for d in range(len(lb))]
        a = np.frombuffer(recs[k + 3], dtype="<f8").reshape(dims, order="F")
        out.append((i1, i2, j1, j2, lb, ub, a))
    return out


def place(canvas, ox, oy, a, lb):
    """Copy a dumped array into `canvas` (whose origin is global (ox,oy))."""
    ix = lb[0] - ox
    iy = lb[1] - oy
    canvas[ix:ix + a.shape[0], iy:iy + a.shape[1]] = a


def main(seq_path, par_paths):
    seq = checkpoints(seq_path)
    par = [checkpoints(p) for p in par_paths]
    n = len(seq)
    for p in par:
        assert len(p) == n, (len(p), n)
    print(f"{n} checkpoints in each file")

    for k in range(n):
        i1, i2, j1, j2, lb, ub, ref = seq[k]
        # global window covering all dumps of this checkpoint
        all_blocks = [(lb, ref)]
        for r in range(len(par)):
            all_blocks.append((par[r][k][4], par[r][k][6]))
        ox = min(b[0][0] for b in all_blocks)
        oy = min(b[0][1] for b in all_blocks)
        ex = max(b[0][0] + b[1].shape[0] for b in all_blocks)
        ey = max(b[0][1] + b[1].shape[1] for b in all_blocks)
        shape = [ex - ox, ey - oy] + list(ref.shape[2:])
        A = np.full(shape, np.nan)
        B = np.full(shape, np.nan)
        place(A, ox, oy, ref, lb)
        for r in range(len(par)):
            plb = par[r][k][4]
            place(B, ox, oy, par[r][k][6], plb)
        both = ~np.isnan(A) & ~np.isnan(B)
        d = np.where(both, np.abs(np.nan_to_num(A) - np.nan_to_num(B)), 0.0)

        ii = np.arange(ox, ex)
        jj = np.arange(oy, ey)
        interior = np.zeros(d.shape[:2], dtype=bool)
        interior[np.ix_((ii >= 1) & (ii <= LM), (jj >= 1) & (jj <= MM))] = True
        if d.ndim == 3:
            interior = np.repeat(interior[:, :, None], d.shape[2], axis=2)
        dint = np.where(interior, d, 0.0)

        ndiff = int((d > 0).sum())
        ndiff_int = int((dint > 0).sum())
        if ndiff == 0:
            print(f"[{k:4d}] identical (window i=[{ox},{ex-1}] j=[{oy},{ey-1}])")
            continue
        print(
            f"[{k:4d}] DIFFERS n={ndiff} (interior n={ndiff_int}) "
            f"max_abs={d.max():.3e} max_abs_interior={dint.max():.3e}"
        )
        if ndiff_int == 0:
            print("       -> ghost cells only")
            continue
        idx = np.unravel_index(int(dint.ravel().argmax()), dint.shape)
        extra = tuple(int(v) for v in idx[2:]) if dint.ndim == 3 else ()
        print(
            f"       -> largest interior diff at global i={ii[idx[0]]} "
            f"j={jj[idx[1]]} k={extra}: seq={A[idx]:.17e} mpi={B[idx]:.17e} "
            f"diff={d[idx]:.3e}"
        )
        return k
    print("no interior difference found in any checkpoint")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
