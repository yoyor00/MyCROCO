#!/usr/bin/env python3

"""
Plot results from the SINGLE_COLUMN test cases.

Produces:
- A top panel with forcing time series (surface stress, heat flux,
  bottom stress) when these variables are present in the file.
- Hovmöller diagrams (time vs depth) of the non-trivial 3D fields
  (temperature, velocity, turbulent viscosity/diffusivity).

The script auto-detects which fields contain meaningful data and
only plots those, so it works for all seven SCM variants.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# ── 3D fields to inspect (Hovmöller), in display order ──
#
# (netcdf_name, label, unit, colormap, symmetric, grid)

CANDIDATE_3D = [
    ("temp", "Temperature", "°C", "RdYlBu_r", False, "r"),
    ("u", "U-velocity", "m/s", "RdBu_r", True, "r"),
    ("v", "V-velocity", "m/s", "RdBu_r", True, "r"),
    ("AKv", "Vertical viscosity (Akv)", "m²/s", "viridis", False, "w"),
    ("AKt", "Vertical diffusivity (Akt)", "m²/s", "viridis", False, "w"),
]

# ── 2D forcing fields to inspect (time series), in display order ──
#
# (netcdf_name, label, unit, color, group)
#   group: 'surface' or 'bottom' — plotted on separate sub-axes

CANDIDATE_FORCING = [
    ("sustr", "Surface u-stress (τx)", "N/m²", "tab:blue", "surface"),
    ("svstr", "Surface v-stress (τy)", "N/m²", "tab:cyan", "surface"),
    ("shflux", "Net heat flux", "W/m²", "tab:red", "surface"),
    ("swrad", "Shortwave radiation", "W/m²", "tab:orange", "surface"),
    ("bostr", "Bottom stress", "N/m²", "tab:brown", "bottom"),
    ("bustr", "Bottom u-stress", "N/m²", "tab:brown", "bottom"),
    ("bvstr", "Bottom v-stress", "N/m²", "tab:olive", "bottom"),
]

# ── CLI ──────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Plot results from the SINGLE_COLUMN test case.",
    epilog="Example usage:\n"
    + "  python plot_single_column.py --file scm_his.nc\n"
    + "  python plot_single_column.py --makepng --no-show --output-dir ./plots",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="scm_his.nc",
    help="Path to the NetCDF history file (default: scm_his.nc)",
)
parser.add_argument(
    "--makepdf",
    action="store_true",
    help="Generate a PDF of the plots",
)
parser.add_argument(
    "--makepng",
    action="store_true",
    help="Generate a PNG of the plots",
)
parser.add_argument(
    "--no-show",
    action="store_true",
    help="Do not display the plots on the screen",
)
parser.add_argument(
    "--output-dir",
    type=str,
    default=".",
    help="Directory to save the output files (default: current directory)",
)
args = parser.parse_args()


# ── Helpers ──────────────────────────────────────────────


def is_nontrivial(arr, rtol=1e-10, require_variation=True):
    """Return True if the array has meaningful content."""
    arr = np.asarray(arr, dtype=float)
    if np.all(np.isnan(arr)):
        return False
    vmin, vmax = np.nanmin(arr), np.nanmax(arr)
    span = vmax - vmin
    scale = max(abs(vmax), abs(vmin), 1.0)
    has_variation = span / scale > rtol
    has_signal = abs(vmax) > rtol
    if require_variation:
        return has_variation  # for 3D fields: must vary
    else:
        return has_variation or has_signal  # for forcing: non-zero suffit


# ── Read data ────────────────────────────────────────────

try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# Time in hours
time = nc.variables["scrum_time"][:] / 3600.0

# Grid parameters
h = nc.variables["h"][:]
theta_s = float(nc.theta_s)
theta_b = float(nc.theta_b)
N = len(nc.dimensions["s_rho"])

vtrans_var = nc.variables.get("Vtransform", None)
vtransform = int(np.squeeze(vtrans_var[:])) if vtrans_var is not None else 2

hc_var = nc.variables.get("hc", None)
if hc_var is not None:
    hc = float(hc_var[:])
else:
    hc = float(nc.hc) if hasattr(nc, "hc") else float(nc.variables["Tcline"][:])

# Central column indices
j0 = h.shape[0] // 2
i0 = h.shape[1] // 2

# Compute depths at rho and w points for central column
zeta0 = np.zeros_like(h)
zr = cr.zlevs(h, zeta0, theta_s, theta_b, hc, N, "r", vtransform)
zw = cr.zlevs(h, zeta0, theta_s, theta_b, hc, N, "w", vtransform)
z_col_r = zr[:, j0, i0]
z_col_w = zw[:, j0, i0]

# ── Scan forcing fields (2D time series) ─────────────────

forcing_curves = []
seen_forcing = set()  # avoid duplicates (sustr vs Ustr)

for vname, label, unit, color, group in CANDIDATE_FORCING:
    if vname not in nc.variables:
        continue

    var = nc.variables[vname]
    # Forcing fields are (time, eta, xi) — extract central column
    if len(var.shape) == 3:
        raw = np.squeeze(var[:, j0, i0])
    elif len(var.shape) == 2:
        j = min(j0, var.shape[-1] - 1)
        raw = np.squeeze(var[:, j])
    elif len(var.shape) == 1:
        raw = var[:]
    else:
        continue

    if not is_nontrivial(raw, require_variation=True):
        continue

    # Avoid plotting both sustr and Ustr (same physical quantity)
    dedup_key = (group, label.split("(")[0].strip())
    if dedup_key in seen_forcing:
        continue
    seen_forcing.add(dedup_key)

    forcing_curves.append(
        {
            "vname": vname,
            "values": raw,
            "label": label,
            "unit": unit,
            "color": color,
            "group": group,
        }
    )
    print(f"  Forcing: '{vname}' ({label}).")

# Group forcing by surface/bottom
forcing_groups = {}
for curve in forcing_curves:
    g = curve["group"]
    if g not in forcing_groups:
        forcing_groups[g] = []
    forcing_groups[g].append(curve)

# ── Scan 3D fields (Hovmöller) ──────────────────────────

hovmoller_panels = []
for vname, label, unit, cmap, sym, grid in CANDIDATE_3D:
    if vname not in nc.variables:
        continue

    var = nc.variables[vname]
    if len(var.shape) < 3:
        continue

    j = min(j0, var.shape[-2] - 1)
    i = min(i0, var.shape[-1] - 1)
    raw = np.squeeze(var[:, :, j, i])

    if not is_nontrivial(raw, require_variation=True):
        print(f"  Skipping '{vname}' (uniform/zero).")
        continue

    z_col = z_col_w if grid == "w" else z_col_r

    hovmoller_panels.append(
        {
            "vname": vname,
            "values": raw,
            "z": z_col,
            "label": label,
            "unit": unit,
            "cmap": cmap,
            "symmetric": sym,
        }
    )
    print(f"  Plotting '{vname}'.")

# Title from NetCDF attribute
title = nc.title.strip() if hasattr(nc, "title") else "SINGLE COLUMN"
nc.close()

# ── Layout (2-column GridSpec: plots | colorbars) ────────
#
# All plot axes occupy the left column with the same width.
# Hovmöller colorbars go in the right column.
# Forcing panels leave the right column empty.

n_forcing = len(forcing_groups)  # 0, 1, or 2 panels
n_hovmoller = len(hovmoller_panels)
n_total = n_forcing + n_hovmoller

if n_total == 0:
    print("No fields to plot.")
    exit(0)

# Height ratios: forcing panels are shorter (1), Hovmöller panels are taller (3)
heights = [1] * n_forcing + [3] * n_hovmoller

fig = plt.figure(figsize=(10, sum(heights) * 1.0 + 1))
gs = fig.add_gridspec(
    n_total, 2, height_ratios=heights, width_ratios=[1, 0.03], hspace=0.35, wspace=0.05
)

panel_idx = 0
plot_axes = []  # all left-column axes, in order

# ── Forcing panels (time series) ─────────────────────────

for group_name in ["surface", "bottom"]:
    if group_name not in forcing_groups:
        continue

    ax = fig.add_subplot(gs[panel_idx, 0])
    plot_axes.append(ax)
    curves = forcing_groups[group_name]

    for curve in curves:
        ax.plot(
            time,
            curve["values"],
            color=curve["color"],
            label=f"{curve['label']} [{curve['unit']}]",
        )

    ax.set_ylabel(curves[0]["unit"])
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)

    if panel_idx == 0:
        ax.set_title(f"{title}\n{group_name.capitalize()} forcing")
    else:
        ax.set_title(f"{group_name.capitalize()} forcing")

    panel_idx += 1

# ── Hovmöller panels ─────────────────────────────────────

for panel in hovmoller_panels:
    ax = fig.add_subplot(gs[panel_idx, 0])
    cax = fig.add_subplot(gs[panel_idx, 1])
    plot_axes.append(ax)

    T, Z = np.meshgrid(time, panel["z"], indexing="ij")

    kwargs = {"cmap": panel["cmap"], "shading": "auto"}
    if panel["symmetric"]:
        vmax = np.nanmax(np.abs(panel["values"]))
        if vmax > 0:
            kwargs["vmin"] = -vmax
            kwargs["vmax"] = vmax

    cf = ax.pcolormesh(T, Z, panel["values"], **kwargs)
    fig.colorbar(cf, cax=cax, label=panel["unit"])
    ax.set_ylabel("Depth (m)")

    if panel_idx == 0:
        ax.set_title(f"{title}\n{panel['label']}")
    else:
        ax.set_title(panel["label"])

    panel_idx += 1

# ── Align all x-axes ─────────────────────────────────────

for ax in plot_axes:
    ax.set_xlim(time[0], time[-1])

# Hide tick labels on all but the last panel
for ax in plot_axes[:-1]:
    ax.tick_params(labelbottom=False)

plot_axes[-1].set_xlabel("Time (hours)")

# ── Save / Show ──────────────────────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

if args.makepdf:
    pdf_path = os.path.join(output_dir, "single_column_plots.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "single_column_plots.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()
