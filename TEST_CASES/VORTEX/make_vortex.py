#!/usr/bin/env python3
"""
Build a CROCO VORTEX configuration.
Create a grid file and initial/climatology files for the vortex experiment
(for the parent and the child grid).
The vortex is defined such that there is no motion below a defined depth H0.
"""

import argparse
import os
import sys
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import netCDF4 as nc4
import numpy as np
from scipy.interpolate import RectBivariateSpline

# croco_utils.py lives one level up (TEST_CASES/)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from croco_utils import scoordinate, tridim, u2rho_2d, v2rho_2d, zlevs  # noqa: I001


# ═══════════════════════════════════════════════════════════════════════════════
# USER-DEFINED PARAMETERS  — edit here
# ═══════════════════════════════════════════════════════════════════════════════

title = "VORTEX"
parent_grd = "vortex_grd.nc"
child_grd = "vortex_grd.nc.1"
parent_ini = "vortex_ini.nc"
child_ini = "vortex_ini.nc.1"
parent_clm = "vortex_clm.nc"
child_clm = "vortex_clm.nc.1"

dx = 30e3  # Horizontal resolution (m)
xmax = 900e3  # Domain half-length (m)
H0 = 5000.0  # Depth (m)
H = 2500.0  # Level of no-motion (m)
theta = 38.5  # Latitude (degrees, beta-plane)
R = 6367442.76  # Earth radius (m)
Pa = 1013e2  # Atmospheric pressure (Pa)
rho0 = 1024.4  # Mean ocean density (kg/m³)
umax = 1.0  # Max velocity (m/s)
lam = np.sqrt(2) * 60e3  # Vortex radius (m)
g = 9.81  # Gravity (m/s²)
N2_sq = (0.003) ** 2  # Brunt-Väisälä frequency squared (s⁻²)

# Vertical grid (vtransform=2 = NEW_S_COORD)
N = 10
theta_s = 1.0
theta_b = 0.0
hc = H0  # no hc=min(hmin,hc) constraint for vtransform=2
vtransform = 2

# Nesting
refinecoeff = 3
jmin, jmax = 21, 40
imin, imax = 21, 40

# ═══════════════════════════════════════════════════════════════════════════════


# ── Helpers ──────────────────────────────────────────────────────────────────


def write_grd(fname, title, X, Y, h, f, dx, dy, *, grd_pos=None, refine_coef=None):
    """Write a CROCO grid NetCDF file."""
    Mp, Lp = h.shape
    M, L = Mp - 1, Lp - 1

    with nc4.Dataset(fname, "w") as ds:
        ds.title = title
        ds.date = datetime.now(timezone.utc).strftime("%d-%b-%Y")
        ds.type = "CROCO grid file"

        ds.createDimension("xi_rho", Lp)
        ds.createDimension("eta_rho", Mp)
        ds.createDimension("xi_psi", L)
        ds.createDimension("eta_psi", M)
        ds.createDimension("one", 1)
        if grd_pos is not None:
            ds.createDimension("four", 4)

        ds.createVariable("xl", "f8", ("one",))[:] = dx * (L - 1)
        ds.createVariable("el", "f8", ("one",))[:] = dy * (M - 1)
        sph = ds.createVariable("spherical", "S1", ("one",))
        sph[0] = b"F"
        ds.createVariable("h", "f8", ("eta_rho", "xi_rho"))[:] = h
        ds.createVariable("f", "f8", ("eta_rho", "xi_rho"))[:] = f
        ds.createVariable("pm", "f8", ("eta_rho", "xi_rho"))[:] = 1.0 / dx
        ds.createVariable("pn", "f8", ("eta_rho", "xi_rho"))[:] = 1.0 / dy
        ds.createVariable("x_rho", "f8", ("eta_rho", "xi_rho"))[:] = X
        ds.createVariable("y_rho", "f8", ("eta_rho", "xi_rho"))[:] = Y
        ds.createVariable("mask_rho", "f8", ("eta_rho", "xi_rho"))[:] = np.ones(
            (Mp, Lp)
        )

        if grd_pos is not None:
            # grd_pos: [imin, imax, jmin, jmax] position of this child in parent grid (AGRIF)
            ds.createVariable("grd_pos", "i4", ("four",))[:] = grd_pos
            ds.createVariable("refine_coef", "i4", ("one",))[:] = refine_coef


def write_ini(fname, grd_fname, title, theta_s, theta_b, hc, N, vtransform):
    """Write a CROCO initial-conditions NetCDF file (fields initialised to zero)."""
    with nc4.Dataset(grd_fname, "r") as grd:
        Mp, Lp = grd.variables["h"].shape

    M, L = Mp - 1, Lp - 1
    sc_r, Cs_r, _sc_w, _Cs_w = scoordinate(theta_s, theta_b, N, hc, vtransform)
    Vstretching = 4 if vtransform == 2 else 1

    with nc4.Dataset(fname, "w") as ds:
        ds.title = title
        ds.date = datetime.now(timezone.utc).strftime("%d-%b-%Y")
        ds.clim_file = fname
        ds.grd_file = grd_fname
        ds.type = "INITIAL file"
        ds.history = "CROCO"

        ds.createDimension("one", 1)
        ds.createDimension("s_rho", N)
        ds.createDimension("time", None)  # unlimited
        ds.createDimension("eta_rho", Mp)
        ds.createDimension("xi_rho", Lp)
        ds.createDimension("eta_u", Mp)
        ds.createDimension("xi_u", L)
        ds.createDimension("eta_v", M)
        ds.createDimension("xi_v", Lp)

        sph = ds.createVariable("spherical", "S1", ("one",))
        sph[0] = b"F"
        ds.createVariable("Vtransform", "i4", ("one",))[0] = vtransform
        ds.createVariable("Vstretching", "i4", ("one",))[0] = Vstretching
        ds.createVariable("theta_s", "f8", ("one",))[0] = theta_s
        ds.createVariable("theta_b", "f8", ("one",))[0] = theta_b
        ds.createVariable("Tcline", "f8", ("one",))[0] = hc
        ds.createVariable("hc", "f8", ("one",))[0] = hc
        ds.createVariable("sc_r", "f8", ("s_rho",))[:] = sc_r
        ds.createVariable("Cs_r", "f8", ("s_rho",))[:] = Cs_r
        ot = ds.createVariable("ocean_time", "f8", ("time",))
        ot[0] = 0.0
        ot.units    = "second"
        ot.calendar = "360.0 days in every year"
        st = ds.createVariable("scrum_time", "f8", ("time",))
        st[0] = 0.0
        st.units    = "second"
        st.calendar = "360.0 days in every year"

        ds.createVariable("u", "f8", ("time", "s_rho", "eta_u", "xi_u"))
        ds.createVariable("v", "f8", ("time", "s_rho", "eta_v", "xi_v"))
        ds.createVariable("ubar", "f8", ("time", "eta_u", "xi_u"))
        ds.createVariable("vbar", "f8", ("time", "eta_v", "xi_v"))
        ds.createVariable("zeta", "f8", ("time", "eta_rho", "xi_rho"))
        ds.createVariable("temp", "f8", ("time", "s_rho", "eta_rho", "xi_rho"))
        ds.createVariable("salt", "f8", ("time", "s_rho", "eta_rho", "xi_rho"))


def write_clm(
    fname,
    grd_fname,
    title,
    theta_s,
    theta_b,
    hc,
    N,
    clim_times,
    cycle_length,
    vtransform,
):
    """Write a CROCO climatology NetCDF file (fields initialised to zero)."""
    with nc4.Dataset(grd_fname, "r") as grd:
        Mp, Lp = grd.variables["h"].shape

    M, L = Mp - 1, Lp - 1
    nt = len(clim_times)
    sc_r, Cs_r, sc_w, Cs_w = scoordinate(theta_s, theta_b, N, hc, vtransform)
    Vstretching = 4 if vtransform == 2 else 1
    times = np.asarray(clim_times, dtype=float)

    with nc4.Dataset(fname, "w") as ds:
        ds.title = title
        ds.date = datetime.now(timezone.utc).strftime("%d-%b-%Y")
        ds.clim_file = fname
        ds.grd_file = grd_fname
        ds.type = "CLIMATOLOGY file"
        ds.history = "CROCO"

        ds.createDimension("one", 1)
        ds.createDimension("s_rho", N)
        ds.createDimension("s_w", N + 1)
        ds.createDimension("eta_rho", Mp)
        ds.createDimension("xi_rho", Lp)
        ds.createDimension("eta_u", Mp)
        ds.createDimension("xi_u", L)
        ds.createDimension("eta_v", M)
        ds.createDimension("xi_v", Lp)
        for tname in (
            "tclm_time",
            "sclm_time",
            "uclm_time",
            "vclm_time",
            "ssh_time",
            "zeta_time",
            "temp_time",
            "salt_time",
            "v2d_time",
            "v3d_time",
        ):
            ds.createDimension(tname, nt)

        sph = ds.createVariable("spherical", "S1", ("one",))
        sph[0] = b"F"
        ds.createVariable("Vtransform", "i4", ("one",))[0] = vtransform
        ds.createVariable("Vstretching", "i4", ("one",))[0] = Vstretching
        ds.createVariable("theta_s", "f8", ("one",))[0] = theta_s
        ds.createVariable("theta_b", "f8", ("one",))[0] = theta_b
        ds.createVariable("Tcline", "f8", ("one",))[0] = hc
        ds.createVariable("hc", "f8", ("one",))[0] = hc
        ds.createVariable("sc_r", "f8", ("s_rho",))[:] = sc_r
        ds.createVariable("sc_w", "f8", ("s_w",))[:] = sc_w
        ds.createVariable("Cs_r", "f8", ("s_rho",))[:] = Cs_r
        ds.createVariable("Cs_w", "f8", ("s_w",))[:] = Cs_w

        tvar_longnames = {
            "tclm_time": "time for temperature climatology",
            "sclm_time": "time for salinity climatology",
            "uclm_time": "time for u-momentum climatology",
            "vclm_time": "time for v-momentum climatology",
            "ssh_time": "time for sea surface height",
            "zeta_time": "time for sea surface height",
            "temp_time": "time for temperature climatology",
            "salt_time": "time for salinity climatology",
            "v2d_time": "time for 2D velocity climatology",
            "v3d_time": "time for 3D velocity climatology",
        }
        for tname in (
            "tclm_time",
            "sclm_time",
            "uclm_time",
            "vclm_time",
            "ssh_time",
            "zeta_time",
            "temp_time",
            "salt_time",
            "v2d_time",
            "v3d_time",
        ):
            tv = ds.createVariable(tname, "f8", (tname,))
            tv[:] = times
            tv.long_name = tvar_longnames[tname]
            tv.units = "day"
            tv.calendar = "360.0 days in every year"
            tv.cycle_length = float(cycle_length)

        ds.createVariable("temp", "f8", ("tclm_time", "s_rho", "eta_rho", "xi_rho"))
        ds.createVariable("salt", "f8", ("sclm_time", "s_rho", "eta_rho", "xi_rho"))
        ds.createVariable("u", "f8", ("uclm_time", "s_rho", "eta_u", "xi_u"))
        ds.createVariable("v", "f8", ("vclm_time", "s_rho", "eta_v", "xi_v"))
        ds.createVariable("ubar", "f8", ("uclm_time", "eta_u", "xi_u"))
        ds.createVariable("vbar", "f8", ("vclm_time", "eta_v", "xi_v"))
        ds.createVariable("SSH", "f8", ("ssh_time", "eta_rho", "xi_rho"))
        ds.createVariable("zeta", "f8", ("zeta_time", "eta_rho", "xi_rho"))


def fill_ini(fname, u, v, ubar, vbar, zeta, t):
    """Write vortex fields into an already-created initial file."""
    with nc4.Dataset(fname, "r+") as ds:
        ds["u"][0] = u
        ds["v"][0] = v
        ds["ubar"][0] = ubar
        ds["vbar"][0] = vbar
        ds["zeta"][0] = zeta
        ds["temp"][0] = t
        ds["salt"][0] = np.zeros_like(t)


def fill_clm(fname, u, v, ubar, vbar, zeta, t):
    """Write vortex fields into an already-created climatology file (both time steps)."""
    with nc4.Dataset(fname, "r+") as ds:
        for i in range(ds.dimensions["tclm_time"].size):
            ds["u"][i] = u
            ds["v"][i] = v
            ds["ubar"][i] = ubar
            ds["vbar"][i] = vbar
            ds["zeta"][i] = zeta
            ds["SSH"][i] = zeta
            ds["temp"][i] = t


def write_agrif(imin, imax, jmin, jmax, rcoeff):
    """Write AGRIF_FixedGrids.in."""
    fname = "AGRIF_FixedGrids.in"
    with open(fname, "w") as f:
        f.write("    1\n")
        f.write(
            f"    {imin}    {imax}    {jmin}    {jmax}"
            f"    {rcoeff}    {rcoeff}    {rcoeff}    {rcoeff}\n"
        )
        f.write("    0\n")
        f.write("# number of children per parent\n")
        f.write("# imin imax jmin jmax spacerefx spacerefy timerefx timerefy\n")
        f.write("# [all coordinates are relative to each parent grid!]\n")
        f.write("~\n")
    print(f"  Written {fname}")


def interp_child(field, eta_p, xi_p, jrchild, irchild):
    """Cubic interpolation of a parent-grid field onto child rho points."""
    spl = RectBivariateSpline(eta_p, xi_p, field, kx=3, ky=3)
    return spl(jrchild, irchild)


# ── Vortex physics ────────────────────────────────────────────────────────────


def barocvortex(
    X,
    Y,
    h0,
    theta_s,
    theta_b,
    hc,
    N,
    vtransform,
    rho0,
    Pa,
    f0,
    umax,
    lam,
    g,
    N2_sq,
    H,
    geostrophic=False,
):
    """
    Compute zeta, u, v, ubar, vbar, temp for the baroclinic vortex.

    Parameters
    ----------
    geostrophic : bool
        False (default) — gradient-wind correction (current barocvortex.m).
        True  — pure geostrophy evaluated directly at u/v points.

    Returns
    -------
    zeta : (Mp, Lp)
    u    : (N, Mp, L)
    v    : (N, M, Lp)
    ubar : (Mp, L)
    vbar : (M, Lp)
    t    : (N, Mp, Lp)
    xr   : (N, Mp, Lp)   3D x-coordinates
    zr   : (N, Mp, Lp)   3D z-coordinates
    """
    P0 = rho0 * f0 * umax * lam * np.sqrt(np.e / 2.0)
    r2 = X**2 + Y**2
    P1 = Pa + P0 * np.exp(-r2 / lam**2)

    a = -P0 * (1.0 - np.exp(-H)) / (g * (H - 1.0 + np.exp(-H)))
    rho1 = rho0 + a * np.exp(-r2 / lam**2)
    zeta = (P1 - Pa) / (g * rho1)

    zw = zlevs(h0, zeta, theta_s, theta_b, hc, N, "w", vtransform)  # (N+1, Mp, Lp)
    zr = zlevs(h0, zeta, theta_s, theta_b, hc, N, "r", vtransform)  # (N, Mp, Lp)

    xr = tridim(X, N)  # (N, Mp, Lp)
    yr = tridim(Y, N)

    # Density and temperature
    # Clip exp argument to avoid overflow; values where zr < -H are masked out anyway
    exp_arg = np.clip(-zr - H, None, 500.0)
    rho = rho0 * (1.0 - N2_sq * zr / g)
    rhodyn = (
        -P0
        * (1.0 - np.exp(exp_arg))
        * np.exp(-(xr**2 + yr**2) / lam**2)
        / (g * (H - 1.0 + np.exp(-H)))
    )
    rho[zr > -H] += rhodyn[zr > -H]

    R0, TCOEF = 30.0, 0.28
    t = (-rho + 1000.0 + R0) / TCOEF

    if geostrophic:
        # Pure geostrophy evaluated directly at staggered u/v points.
        # Matches the old commented-out formula in barocvortex.m and the
        # reference .nc files shipped in the repo.
        a_coef = 2.0 * P0 / (f0 * rho0 * lam**2)

        zu = 0.5 * (zr[:, :, :-1] + zr[:, :, 1:])  # (N, Mp, L)
        xu = 0.5 * (xr[:, :, :-1] + xr[:, :, 1:])
        yu = 0.5 * (yr[:, :, :-1] + yr[:, :, 1:])
        exp_u = np.clip(-zu - H, None, 500.0)
        Fu = (H - 1.0 + zu + np.exp(exp_u)) / (H - 1.0 + np.exp(-H))
        Fu[zu < -H] = 0.0
        u = a_coef * Fu * yu * np.exp(-(xu**2 + yu**2) / lam**2)

        zv = 0.5 * (zr[:, :-1, :] + zr[:, 1:, :])  # (N, M, Lp)
        xv = 0.5 * (xr[:, :-1, :] + xr[:, 1:, :])
        yv = 0.5 * (yr[:, :-1, :] + yr[:, 1:, :])
        exp_v = np.clip(-zv - H, None, 500.0)
        Fv = (H - 1.0 + zv + np.exp(exp_v)) / (H - 1.0 + np.exp(-H))
        Fv[zv < -H] = 0.0
        v = -a_coef * Fv * xv * np.exp(-(xv**2 + yv**2) / lam**2)

    else:
        # Gradient-wind correction (current barocvortex.m).
        F = (H - 1.0 + zr + np.exp(exp_arg)) / (H - 1.0 + np.exp(-H))
        F[zr < -H] = 0.0
        r = np.sqrt(xr**2 + yr**2)
        r = np.where(r == 0.0, 1e-10, r)  # avoid 0/0

        Vg = -(2.0 * r * P0 * F / (rho0 * f0 * lam**2)) * np.exp(-(r**2) / lam**2)
        a_c = 1.0 + 4.0 * Vg / (f0 * r)
        n_neg = int(np.sum(a_c < 0))
        if n_neg:
            print(
                f"  {n_neg} points with no gradient-wind solution (fall back to geostrophy)"
            )
        a_c = np.where(a_c < 0, 1.0, a_c)
        Vgr = 2.0 * Vg / (1.0 + np.sqrt(a_c))

        ur_3d = -Vgr * yr / r  # (N, Mp, Lp)
        vr_3d = Vgr * xr / r

        u = 0.5 * (ur_3d[:, :, :-1] + ur_3d[:, :, 1:])  # (N, Mp, L)
        v = 0.5 * (vr_3d[:, :-1, :] + vr_3d[:, 1:, :])  # (N, M, Lp)

    # Barotropic speeds
    dz = zw[1:] - zw[:-1]  # (N, Mp, Lp)
    dzu = 0.5 * (dz[:, :, :-1] + dz[:, :, 1:])  # (N, Mp, L)
    dzv = 0.5 * (dz[:, :-1, :] + dz[:, 1:, :])  # (N, M, Lp)
    ubar = np.sum(dzu * u, axis=0) / np.sum(dzu, axis=0)  # (Mp, L)
    vbar = np.sum(dzv * v, axis=0) / np.sum(dzv, axis=0)  # (M, Lp)

    return zeta, u, v, ubar, vbar, t, xr, zr


# ── Main ──────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Build CROCO VORTEX grid, initial and climatology files.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("--no-show", action="store_true", help="Suppress plots")
    parser.add_argument(
        "--no-plots", action="store_true", help="Skip plotting entirely"
    )
    parser.add_argument(
        "--geostrophic",
        action="store_true",
        help="Use pure geostrophic formula for u/v instead of default "
        "gradient-wind correction (current barocvortex.m).",
    )
    args = parser.parse_args()

    dy = dx

    # Coriolis (beta-plane)
    deg2rad = np.pi / 180.0
    omega = 2.0 * np.pi / (24.0 * 3600.0)
    f0 = 2.0 * omega * np.sin(deg2rad * theta)
    beta = 2.0 * omega / R * np.cos(deg2rad * theta)

    # ── Parent grid ──────────────────────────────────────────────────────────
    print("Building parent grid …")
    x = np.arange(-xmax - dx / 2.0, xmax + dx / 2.0 + dx, dx)
    X, Y = np.meshgrid(x, x)
    h0 = H0 * np.ones_like(X)
    f = f0 + beta * Y

    write_grd(parent_grd, title, X, Y, h0, f, dx, dy)

    # Parent vortex fields
    print("Computing parent vortex …")
    zeta, u, v, ubar, vbar, t, xr, zr = barocvortex(
        X,
        Y,
        h0,
        theta_s,
        theta_b,
        hc,
        N,
        vtransform,
        rho0,
        Pa,
        f0,
        umax,
        lam,
        g,
        N2_sq,
        H,
        geostrophic=args.geostrophic,
    )

    # Parent ini
    write_ini(parent_ini, parent_grd, title, theta_s, theta_b, hc, N, vtransform)
    fill_ini(parent_ini, u, v, ubar, vbar, zeta, t)

    # Parent clm
    write_clm(
        parent_clm,
        parent_grd,
        title,
        theta_s,
        theta_b,
        hc,
        N,
        [25.0, 75.0],
        100,
        vtransform,
    )
    fill_clm(parent_clm, u, v, ubar, vbar, zeta, t)

    # ── Child grid ───────────────────────────────────────────────────────────
    print("Building child grid …")
    with nc4.Dataset(parent_grd, "r") as grd_p:
        xr_p = grd_p["x_rho"][:]
        yr_p = grd_p["y_rho"][:]
        f_p = grd_p["f"][:]
        pm_p = grd_p["pm"][:]
        pn_p = grd_p["pn"][:]

    Mp_p, Lp_p = xr_p.shape

    # 1-indexed coordinate arrays matching MATLAB interp2 source grids
    xi_p = np.arange(1, Lp_p + 1, dtype=float)
    eta_p = np.arange(1, Mp_p + 1, dtype=float)

    # Child index arrays (rho and psi points)
    # rho points: imin+0.5*(1-1/r) : 1/r : imax+0.5*(1+1/r)  → (imax-imin)*r + 2 points
    step = 1.0 / refinecoeff
    n_rho_child = (imax - imin) * refinecoeff + 2  # = 59 for defaults
    irchild = np.linspace(imin + 0.5 - 0.5 * step, imax + 0.5 + 0.5 * step, n_rho_child)
    jrchild = np.linspace(jmin + 0.5 - 0.5 * step, jmax + 0.5 + 0.5 * step, n_rho_child)

    # Cubic interpolation onto child rho grid
    xrchild = interp_child(xr_p, eta_p, xi_p, jrchild, irchild)
    yrchild = interp_child(yr_p, eta_p, xi_p, jrchild, irchild)
    fchild = interp_child(f_p, eta_p, xi_p, jrchild, irchild)
    pmchild = interp_child(pm_p, eta_p, xi_p, jrchild, irchild)
    pnchild = interp_child(pn_p, eta_p, xi_p, jrchild, irchild)

    Mp_c, Lp_c = xrchild.shape
    M_c, L_c = Mp_c - 1, Lp_c - 1

    # Override xl/el and pm/pn for child
    dx_c = dx / refinecoeff
    dy_c = dy / refinecoeff

    file_title = (
        f"Grid embedded in {parent_grd} - positions in parent: "
        f"{imin}-{imax}-{jmin}-{jmax}; refine={refinecoeff}"
    )
    print(f"  {file_title}")
    print(f"  Child size: L={L_c}, M={M_c}")

    with nc4.Dataset(child_grd, "w") as ds:
        ds.title = file_title
        ds.date = datetime.now(timezone.utc).strftime("%d-%b-%Y")
        ds.type = "CROCO grid file"

        ds.createDimension("xi_rho", Lp_c)
        ds.createDimension("eta_rho", Mp_c)
        ds.createDimension("xi_psi", L_c)
        ds.createDimension("eta_psi", M_c)
        ds.createDimension("one", 1)
        ds.createDimension("four", 4)

        ds.createVariable("xl", "f8", ("one",))[:] = dx_c * (L_c - 1)
        ds.createVariable("el", "f8", ("one",))[:] = dy_c * (M_c - 1)
        sph = ds.createVariable("spherical", "S1", ("one",))
        sph[0] = b"F"
        ds.createVariable("h", "f8", ("eta_rho", "xi_rho"))[:] = H0
        ds.createVariable("f", "f8", ("eta_rho", "xi_rho"))[:] = fchild
        ds.createVariable("pm", "f8", ("eta_rho", "xi_rho"))[:] = refinecoeff * pmchild
        ds.createVariable("pn", "f8", ("eta_rho", "xi_rho"))[:] = refinecoeff * pnchild
        ds.createVariable("x_rho", "f8", ("eta_rho", "xi_rho"))[:] = xrchild
        ds.createVariable("y_rho", "f8", ("eta_rho", "xi_rho"))[:] = yrchild
        ds.createVariable("mask_rho", "f8", ("eta_rho", "xi_rho"))[:] = np.ones(
            (Mp_c, Lp_c)
        )
        ds.createVariable("grd_pos", "i4", ("four",))[:] = [imin, imax, jmin, jmax]
        ds.createVariable("refine_coef", "i4", ("one",))[0] = refinecoeff

    # Child vortex fields
    print("Computing child vortex …")
    with nc4.Dataset(child_grd, "r") as grd_c:
        h0_c = grd_c["h"][:]
        X_c = grd_c["x_rho"][:]
        Y_c = grd_c["y_rho"][:]

    zeta_c, u_c, v_c, ubar_c, vbar_c, t_c, _xr_c, _zr_c = barocvortex(
        X_c,
        Y_c,
        h0_c,
        theta_s,
        theta_b,
        hc,
        N,
        vtransform,
        rho0,
        Pa,
        f0,
        umax,
        lam,
        g,
        N2_sq,
        H,
        geostrophic=args.geostrophic,
    )

    write_ini(child_ini, child_grd, title, theta_s, theta_b, hc, N, vtransform)
    fill_ini(child_ini, u_c, v_c, ubar_c, vbar_c, zeta_c, t_c)

    write_clm(
        child_clm,
        child_grd,
        title,
        theta_s,
        theta_b,
        hc,
        N,
        [25.0, 75.0],
        100,
        vtransform,
    )
    fill_clm(child_clm, u_c, v_c, ubar_c, vbar_c, zeta_c, t_c)

    # ── AGRIF_FixedGrids.in ──────────────────────────────────────────────────
    write_agrif(imin, imax, jmin, jmax, refinecoeff)

    # ── Plots ────────────────────────────────────────────────────────────────
    if args.no_plots:
        print("Done.")
        return

    # Surface velocity (parent grid)
    ur_surf = u2rho_2d(u[-1, :, :])
    vr_surf = v2rho_2d(v[-1, :, :])
    spd = np.sqrt(ur_surf**2 + vr_surf**2)

    _fig1, ax1 = plt.subplots()
    pc = ax1.pcolormesh(X * 1e-3, Y * 1e-3, spd, shading="auto")
    ax1.quiver(X * 1e-3, Y * 1e-3, ur_surf, vr_surf, color="k", scale=20)
    plt.colorbar(pc, ax=ax1, label="Speed (m/s)")
    ax1.set_aspect("equal")
    ax1.set_xlabel("X (km)")
    ax1.set_ylabel("Y (km)")
    ax1.set_title("Surface velocity – parent grid")
    ax1.get_figure().savefig("vortex_spd.png", dpi=150, bbox_inches="tight")

    # Temperature vertical section at mid-y (parent grid)
    jmid = (X.shape[0] - 2) // 2
    _fig2, ax2 = plt.subplots()
    pc2 = ax2.pcolormesh(
        xr[:, jmid, :] * 1e-3,
        zr[:, jmid, :],
        t[:, jmid, :],
        shading="auto",
        cmap="RdYlBu_r",
    )
    plt.colorbar(pc2, ax=ax2, label="Temp (°C)")
    ax2.set_xlabel("X (km)")
    ax2.set_ylabel("Depth (m)")
    ax2.set_title("Temperature section – parent grid")
    ax2.get_figure().savefig("vortex_temp.png", dpi=150, bbox_inches="tight")

    # SSH (parent grid)
    _fig3, ax3 = plt.subplots()
    pc3 = ax3.pcolormesh(X * 1e-3, Y * 1e-3, zeta, shading="auto")
    plt.colorbar(pc3, ax=ax3, label="SSH (m)")
    ax3.set_aspect("equal")
    ax3.set_xlabel("X (km)")
    ax3.set_ylabel("Y (km)")
    ax3.set_title("SSH – parent grid")
    ax3.get_figure().savefig("vortex_ssh.png", dpi=150, bbox_inches="tight")

    if not args.no_show:
        plt.show()
    else:
        plt.close("all")

    print("Done.")


if __name__ == "__main__":
    main()
