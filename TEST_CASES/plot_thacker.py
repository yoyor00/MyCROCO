#!/usr/bin/env python3
#  CROCO — Coastal and Regional Ocean COmmunity model
#  Copyright (C) 2005-2026 CROCO Development Team
#
#  SPDX-License-Identifier: CECILL-2.1

"""
Plot results from the THACKER test case.

Auto-detects 2DV (THACKER_2DV, slab in x) vs 3D (paraboloid bowl)
and produces the appropriate plots:
  - 2DV: vertical section with U-error + analytical/numerical zeta
  - 3D:  plan views of zeta at several instants + analytical comparison
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# ── CLI ──────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Plot results from the THACKER test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="thacker_his.nc",
    help="Path to the NetCDF file",
)
parser.add_argument(
    "--tindex", type=int, default=None,
    help="Time index to plot (default: last time step)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Save plots as PDF files",
)
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files",
)
parser.add_argument(
    "--no-show", action="store_true", help="Suppress plot display",
)
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files",
)
args = parser.parse_args()
os.makedirs(args.output_dir, exist_ok=True)

# ── Physical parameters (must match ana_initial.F) ───────

eta = 0.1       # nondimensional amplitude
D0 = 10.0       # max depth at rest (m)
Lt = 80.0e3     # length scale (m)
g = 9.81


def zeta_analytical_2dv(xr, t, omega):
    """Analytical zeta for THACKER_2DV (1D, f=0)."""
    return 2 * eta * D0 * (xr * np.cos(omega * t) / Lt - 0.5 * eta / Lt)


def zeta_analytical_3d(xr, yr, t, omega):
    """Analytical zeta for THACKER 3D (paraboloid, f!=0)."""
    return 2 * eta * D0 * (
        xr * np.cos(omega * t) / Lt
        + yr * np.sin(omega * t) / Lt
        - 0.5 * eta / Lt
    )


def u_analytical(t, omega):
    """Analytical u-velocity (uniform, both modes)."""
    return -eta * omega * Lt * np.sin(omega * t)


def v_analytical(t, omega):
    """Analytical v-velocity (uniform, 3D only; zero for 2DV)."""
    return -eta * omega * Lt * np.cos(omega * t)


# ── Read data ────────────────────────────────────────────

try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

scrum_time = nc.variables["scrum_time"][:]
nt = len(scrum_time)
if args.tindex is None:
    tindex = nt - 1
else:
    tindex = min(args.tindex, nt - 1)

# Detect 2DV vs 3D
Mm = nc.dimensions["eta_rho"].size
is_2dv = (Mm <= 3)
print(f"Detected mode: {'2DV' if is_2dv else '3D'} (eta_rho = {Mm})")

# Grid
xr = nc.variables["x_rho"][:]    # (M, L) in meters
yr = nc.variables["y_rho"][:]
hr = nc.variables["h"][:]
f0 = float(nc.variables["f"][0, 0])
Dcrit = float(nc.variables["Dcrit"][:])

# Vertical grid parameters
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
hc = float(nc.hc)
N = len(nc.dimensions["s_rho"])

# Omega
omega = 0.5 * f0 + np.sqrt(0.25 * f0**2 + 2 * g * D0 / Lt**2)
T_period = 2 * np.pi / omega

print(f"  f = {f0:.2e}, omega = {omega:.6e}, T = {T_period/3600:.2f} h")


# ── Helper: mask dry cells ───────────────────────────────

def mask_dry(zeta_field, h_field, dcrit):
    """Set dry cells to NaN."""
    out = zeta_field.copy()
    D = h_field + out
    out[D < dcrit + 0.01] = np.nan
    return out


# ==================================================================
#  2DV MODE
# ==================================================================

if is_2dv:
    j0 = Mm // 2
    xr_1d = xr[j0, :] / 1000.0   # km
    hr_1d = hr[j0, :]

    time_val = scrum_time[tindex]

    # Read fields
    zeta_num = np.squeeze(nc.variables["zeta"][tindex, j0, :])
    zeta_num = mask_dry(zeta_num, hr_1d, Dcrit)

    # Analytical
    zeta_ana = zeta_analytical_2dv(xr[j0, :], time_val, omega)
    zeta_ana[zeta_ana < -hr_1d] = np.nan

    # Vertical section: U error
    u_num = np.squeeze(nc.variables["u"][tindex, :, j0, :])
    zeta_for_z = np.squeeze(nc.variables["zeta"][tindex, j0, :])
    zr = cr.zlevs(hr_1d, zeta_for_z, theta_s, theta_b, hc, N, "r", 2)

    u_ana_val = u_analytical(time_val, omega)
    u_ana_2d = np.full_like(u_num, u_ana_val)

    # u is on u-grid (L-1 points), pad to L for display
    u_num = np.concatenate([u_num, u_num[:, -1:]], axis=1)
    u_ana_2d = np.concatenate([u_ana_2d, u_ana_2d[:, -1:]], axis=1)

    # Error in %
    with np.errstate(divide="ignore", invalid="ignore"):
        uerr = np.full_like(u_num, np.nan)
        mask = np.abs(u_ana_2d) > 1e-3
        uerr[mask] = 100 * (u_num[mask] - u_ana_2d[mask]) / u_ana_2d[mask]
    # Mask dry + boundaries
    D_1d = hr_1d + zeta_for_z
    D_2d = np.tile(D_1d, (N, 1))
    uerr[D_2d < Dcrit] = np.nan
    uerr[:, :3] = np.nan
    uerr[:, -3:] = np.nan

    xr_2d = np.tile(xr_1d, (N, 1))

    # ── Read data for second plot (sea level at 60h, 62h, 64h) ──

    # Original instants: every 30min output → index = hours * 2
    t_compare = []
    for hour in [60, 62, 64]:
        ti = hour * 2   # assuming NWRT gives 1 record every 30 min
        if ti < nt:
            t_compare.append(ti)

    if not t_compare:
        # Fallback: use last 3 records if run is too short
        t_compare = list(range(max(0, nt - 3), nt))

    zeta_compare = []
    for ti in t_compare:
        zm = np.squeeze(nc.variables["zeta"][ti, j0, :])
        time_ti = scrum_time[ti]
        za = zeta_analytical_2dv(xr[j0, :], time_ti, omega)
        za[za < -hr_1d] = np.nan
        zm = mask_dry(zm, hr_1d, Dcrit)
        zeta_compare.append((ti, time_ti, zm, za))

    nc.close()

    # ── Plot 1: Vertical section + U error ───────────────

    fig1 = plt.figure(figsize=(12, 7))

    levels = np.linspace(-100, 100, 21)
    cf = plt.contourf(xr_2d, zr, uerr, levels=levels,
                      cmap="RdBu_r", extend="both")
    cbar = plt.colorbar(cf, label="U Error [%]")
    cbar.set_ticks(np.linspace(-100, 100, 11))

    plt.plot(xr_1d, -hr_1d, "k-", linewidth=3)
    plt.plot(xr_1d, zeta_ana, "g-", linewidth=3, label="Analytical")
    plt.plot(xr_1d, zeta_num, "r-", linewidth=3, label="Numerical")

    plt.xlim(xr_1d[0], xr_1d[-1])
    plt.ylim(-D0 - 1, D0 * eta * 3)
    plt.xlabel("X (km)")
    plt.ylabel("Z (m)")
    plt.legend()
    plt.grid(True, alpha=0.3)

    thour = time_val / 3600
    plt.title(f"THACKER 2DV — η (m) and U error (%) at t = {thour:.1f} h",
              fontsize=13, fontweight="bold")

    if args.makepng:
        p = os.path.join(args.output_dir, "thacker_uerr.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PNG saved: {p}")
    if args.makepdf:
        p = os.path.join(args.output_dir, "thacker_uerr.pdf")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PDF saved: {p}")

    # ── Plot 2: Sea level at multiple times ──────────────

    fig2 = plt.figure(figsize=(12, 6))
    plt.plot(xr_1d, -hr_1d, "k-", linewidth=4)

    for idx, (ti, time_ti, zm, za) in enumerate(zeta_compare):
        tlabel = f"{time_ti/3600:.0f} h"
        plt.plot(xr_1d, za, "g-", linewidth=2,
                 label="Analytical" if idx == 0 else "")
        plt.plot(xr_1d, zm, "r-", linewidth=2,
                 label="Numerical" if idx == 0 else "")
        # Label on curve
        valid_pts = np.where(~np.isnan(zm))[0]
        if len(valid_pts) > 10:
            ilbl = valid_pts[int(0.8 * len(valid_pts))]
            plt.text(xr_1d[ilbl], zm[ilbl] + 0.15, tlabel, fontsize=10)

    plt.xlim(20, 100)
    plt.ylim(-3, 3)
    plt.xlabel("X (km)")
    plt.ylabel("Z (m)")
    plt.title("THACKER 2DV — Sea level at various times",
              fontsize=13, fontweight="bold")
    plt.legend()
    plt.grid(True, alpha=0.3)

    if args.makepng:
        p = os.path.join(args.output_dir, "thacker_zcomp.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PNG saved: {p}")
    if args.makepdf:
        p = os.path.join(args.output_dir, "thacker_zcomp.pdf")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PDF saved: {p}")

    if not args.no_show:
        plt.show()
    else:
        plt.close("all")

# ==================================================================
#  3D MODE
# ==================================================================

else:
    xr_km = xr / 1000.0
    yr_km = yr / 1000.0

    # Choose 4 instants: 0, T/4, T/2, tindex
    t_indices = [0, nt // 4, nt // 2, tindex]
    # Remove duplicates, keep order
    seen = set()
    t_indices = [t for t in t_indices if not (t in seen or seen.add(t))]
    npanels = len(t_indices)

    nc_zeta = nc.variables["zeta"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    axes = axes.ravel()

    for idx in range(4):
        ax = axes[idx]

        if idx < npanels:
            ti = t_indices[idx]
            time_val = scrum_time[ti]

            zeta_num = np.squeeze(nc_zeta[ti, :, :])
            zeta_num = mask_dry(zeta_num, hr, Dcrit)

            zeta_ana = zeta_analytical_3d(xr, yr, time_val, omega)
            zeta_ana[zeta_ana < -hr] = np.nan

            # Numerical zeta (color fill)
            vmax = eta * D0 * 2.5 / Lt * Lt  # ~ 2*eta*D0
            vmax = max(np.nanmax(np.abs(zeta_num)), 0.5)
            cf = ax.pcolormesh(xr_km, yr_km, zeta_num,
                               cmap="RdBu_r", shading="auto",
                               vmin=-vmax, vmax=vmax)
            fig.colorbar(cf, ax=ax, label="m", shrink=0.9)

            # Analytical zeta (contour lines)
            clevels = np.linspace(-vmax, vmax, 11)
            clevels = clevels[clevels != 0]
            ax.contour(xr_km, yr_km, zeta_ana,
                       levels=clevels, colors="k", linewidths=0.8,
                       linestyles="dashed")

            # Dry area boundary (h=0 contour = rim of the bowl)
            ax.contour(xr_km, yr_km, hr, levels=[Dcrit],
                       colors="grey", linewidths=1.5)

            thour = time_val / 3600
            ax.set_title(f"t = {thour:.1f} h")
            ax.set_aspect("equal")

            if idx >= 2:
                ax.set_xlabel("X (km)")
            else:
                ax.tick_params(labelbottom=False)
            if idx % 2 == 0:
                ax.set_ylabel("Y (km)")
        else:
            ax.set_visible(False)

    # Read final zeta for error computation before closing
    zeta_final = np.squeeze(nc_zeta[tindex, :, :])

    # ── Read data for vertical section at y=0 (j = Mm/2) ──

    j0 = Mm // 2
    xr_1d = xr[j0, :] / 1000.0   # km
    hr_1d = hr[j0, :]
    time_sec = scrum_time[tindex]

    zeta_sec_num = np.squeeze(nc.variables["zeta"][tindex, j0, :])
    zeta_sec_for_z = zeta_sec_num.copy()
    zeta_sec_num = mask_dry(zeta_sec_num, hr_1d, Dcrit)

    # Analytical zeta along section (y=0 → sin term vanishes)
    zeta_sec_ana = zeta_analytical_2dv(xr[j0, :], time_sec, omega)
    zeta_sec_ana[zeta_sec_ana < -hr_1d] = np.nan

    # Vertical grid along section
    zr_sec = cr.zlevs(hr_1d, zeta_sec_for_z, theta_s, theta_b, hc, N, "r", 2)

    # U velocity along section
    u_sec = np.squeeze(nc.variables["u"][tindex, :, j0, :])
    u_ana_val = u_analytical(time_sec, omega)
    u_ana_sec = np.full_like(u_sec, u_ana_val)

    # Pad u from u-grid (L-1) to rho-grid (L) for display
    u_sec = np.concatenate([u_sec, u_sec[:, -1:]], axis=1)
    u_ana_sec = np.concatenate([u_ana_sec, u_ana_sec[:, -1:]], axis=1)

    # Error in %
    with np.errstate(divide="ignore", invalid="ignore"):
        uerr_sec = np.full_like(u_sec, np.nan)
        mask_u = np.abs(u_ana_sec) > 1e-3
        uerr_sec[mask_u] = 100 * (u_sec[mask_u] - u_ana_sec[mask_u]) / u_ana_sec[mask_u]
    D_1d = hr_1d + zeta_sec_for_z
    D_2d = np.tile(D_1d, (N, 1))
    uerr_sec[D_2d < Dcrit] = np.nan
    uerr_sec[:, :3] = np.nan
    uerr_sec[:, -3:] = np.nan

    xr_sec_2d = np.tile(xr_1d, (N, 1))

    nc.close()

    # ── Figure 1: Plan views ─────────────────────────────

    # Compute error at final time for suptitle
    zeta_ana_final = zeta_analytical_3d(xr, yr, scrum_time[tindex], omega)
    valid = hr > Dcrit
    if np.sum(valid) > 0:
        err_l2 = np.sqrt(np.nanmean((zeta_final[valid] - zeta_ana_final[valid])**2))
        err_linf = np.nanmax(np.abs(zeta_final[valid] - zeta_ana_final[valid]))
    else:
        err_l2 = err_linf = 0.0

    fig.suptitle(
        f"THACKER 3D — Sea surface elevation\n"
        f"Color: numerical, dashed: analytical    "
        f"L2 = {err_l2:.2e} m, L∞ = {err_linf:.2e} m",
        fontsize=12, fontweight="bold",
    )
    plt.subplots_adjust(hspace=0.25, top=0.90)

    if args.makepng:
        p = os.path.join(args.output_dir, "thacker_plan.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PNG saved: {p}")
    if args.makepdf:
        p = os.path.join(args.output_dir, "thacker_plan.pdf")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PDF saved: {p}")

    # ── Figure 2: Vertical section at y=0 ────────────────

    fig2 = plt.figure(figsize=(12, 7))

    levels = np.linspace(-100, 100, 21)
    cf = plt.contourf(xr_sec_2d, zr_sec, uerr_sec, levels=levels,
                      cmap="RdBu_r", extend="both")
    cbar = plt.colorbar(cf, label="U Error [%]")
    cbar.set_ticks(np.linspace(-100, 100, 11))

    plt.plot(xr_1d, -hr_1d, "k-", linewidth=3)
    plt.plot(xr_1d, zeta_sec_ana, "g-", linewidth=3, label="Analytical")
    plt.plot(xr_1d, zeta_sec_num, "r-", linewidth=3, label="Numerical")

    plt.xlim(xr_1d[0], xr_1d[-1])
    plt.ylim(-D0 - 1, D0 * eta * 3)
    plt.xlabel("X (km)")
    plt.ylabel("Z (m)")
    plt.legend()
    plt.grid(True, alpha=0.3)

    thour = time_sec / 3600
    plt.title(f"THACKER 3D — Section at y=0: η (m) and U error (%) "
              f"at t = {thour:.1f} h",
              fontsize=13, fontweight="bold")

    if args.makepng:
        p = os.path.join(args.output_dir, "thacker_section.png")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PNG saved: {p}")
    if args.makepdf:
        p = os.path.join(args.output_dir, "thacker_section.pdf")
        plt.savefig(p, dpi=300, bbox_inches="tight")
        print(f"PDF saved: {p}")

    if not args.no_show:
        plt.show()
    else:
        plt.close("all")
