#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot results from the ESTUARY test case.
Displays the maximum bed shear stress (tauskin) over the simulation duration,
masked where the water depth is zero.
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import xarray
import xarray.plot
import argparse

# ── CLI (standard BENCH-compatible interface) ────────────

parser = argparse.ArgumentParser(
    description="Plot results from the ESTUARY test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="estuary_his.nc",
    help="Path to the NetCDF history file",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate a PDF of the plots",
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate a PNG of the plots",
)
parser.add_argument(
    "--no-show", action="store_true",
    help="Do not display the plots on the screen",
)
parser.add_argument(
    "--output-dir", type=str, default=".",
    help="Directory to save the output files",
)
args = parser.parse_args()

# ── Read data ────────────────────────────────────────────

if not os.path.isfile(args.file):
    print(f"Error: File '{args.file}' not found.")
    exit(1)

print(f"Reading file: {args.file}")
data = xarray.open_dataset(args.file)

tauskin = data.TAUSKIN
bathy = data.h
zeta = data.zeta
hwater = bathy + zeta

dx = data.xl / (data.xi_rho.shape[0] - 2)
dy = data.el / (data.eta_rho.shape[0] - 2)

x_axis = (data.xi_rho - 1) * dx
y_axis = (data.eta_rho - 1) * dy

tauskinmask = tauskin.where(hwater > 0.0)
tauskinmask_max = tauskinmask.max(axis=0)

tauskinmask_max_xy = xarray.DataArray(
    tauskinmask_max[1:-1, 1:-1],
    coords=[y_axis[1:-1], x_axis[1:-1]],
    dims=["y", "x"],
)

# ── Plot ─────────────────────────────────────────────────

cmap_lev = matplotlib.colors.ListedColormap(
    [
        (0.0, 0.0, 0.3),
        (0.2, 0.2, 1.0),
        (0.3, 0.5, 1.0),
        (0.1, 1.0, 1.0),
        (0.2, 0.9, 0.1),
        (1.0, 1.0, 0.0),
        (1.0, 0.7, 0.0),
        (1.0, 0.0, 0.0),
        (0.75, 0.0, 0.0),
        (0.9, 0.0, 0.9),
    ],
    name="from_list",
    N=None,
)
lev = [0, 0.1, 0.2, 0.5, 1, 2, 3, 4, 5, 10]

fig, ax = plt.subplots(figsize=(10, 5))
p = tauskinmask_max_xy.plot.pcolormesh(
    cmap=cmap_lev, add_colorbar=False, levels=lev, extend="max"
)
fmt = lambda x, pos: "{:.4}".format(x)
cb = matplotlib.pyplot.colorbar(p, ax=ax, format=matplotlib.ticker.FuncFormatter(fmt))
cb.ax.tick_params(labelsize=10)
cb.set_ticks(lev)
ax.set_title("Maximum tauskin (N/m²) — dx=%.0fm dy=%.0fm" % (dx, dy))
ax.set_xlabel("X (m)")
ax.set_ylabel("Y (m)")
ax.grid()

plt.suptitle("ESTUARY test case — maximum bed shear stress",
             fontsize=13, fontweight="bold")
plt.subplots_adjust(top=0.88)

# ── Save / Show (do not modify) ──────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

if args.makepdf:
    pdf_path = os.path.join(output_dir, "estuary_plots.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "estuary_plots.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()
