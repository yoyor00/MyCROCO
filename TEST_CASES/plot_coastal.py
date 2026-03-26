#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot results from the COASTAL test case.
Displays the sediment thickness difference (end minus start) and the maximum
suspended sediment concentration for MUD1, SAND and GRAV fractions.
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import xarray
import xarray.plot
import cartopy
import cartopy.crs as ccrs
import argparse


def fmt(x, pos):
    return "{:.4}".format(x)


# ── CLI (standard BENCH-compatible interface) ────────────

parser = argparse.ArgumentParser(
    description="Plot results from the COASTAL test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="vilaine_his.nc",
    help="Path to the NetCDF history file",
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
    help="Directory to save the output files",
)
args = parser.parse_args()

# ── Read data ────────────────────────────────────────────

if not os.path.isfile(args.file):
    print(f"Error: File '{args.file}' not found.")
    exit(1)

print(f"Reading file: {args.file}")
data = xarray.open_dataset(args.file)

hsed = data.HSED
bathy = data.h
lon = data.lon_rho.isel(eta_rho=0)
lat = data.lat_rho.isel(xi_rho=0)

diffhsed = hsed.isel(time=-1) - hsed.isel(time=0)
diffhsedmask = diffhsed.where(bathy > -5.0)
d_hsedmask = xarray.DataArray(
    diffhsedmask[:, :], coords=[lat, lon], dims=["lat", "lon"]
)

# ── Helper functions ──────────────────────────────────────


def add_gridlines(ax):
    """Add formatted lat/lon gridlines to a cartopy axis."""
    gl = ax.gridlines(draw_labels=True, linewidth=0.5, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 9, "color": "gray"}
    gl.ylabel_style = {"size": 9, "color": "gray"}
    gl.xformatter = cartopy.mpl.gridliner.LONGITUDE_FORMATTER
    gl.yformatter = cartopy.mpl.gridliner.LATITUDE_FORMATTER
    gl.xlocator = matplotlib.ticker.AutoLocator()
    gl.ylocator = matplotlib.ticker.AutoLocator()


# ── Plot ─────────────────────────────────────────────────

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

cmap_lev = matplotlib.colors.ListedColormap(
    [
        (0.9, 0.9, 1.0),
        (0.1, 0.1, 1.0),
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
lev = [0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 5, 10]

# Panel 1: sediment thickness difference
fig = plt.figure(figsize=(8, 4))
ax = plt.axes(projection=ccrs.PlateCarree())
d_hsedmask.plot.pcolormesh(
    transform=ccrs.PlateCarree(),
    cmap="RdBu_r",
    add_colorbar=True,
    vmin=-0.001,
    vmax=0.001,
)
ax.set_title(
    "COASTAL test case — HSED difference (end − start)", fontsize=13, fontweight="bold"
)
add_gridlines(ax)

# ── Save / Show (do not modify) ──────────────────────────

if args.makepdf:
    pdf_path = os.path.join(output_dir, "coastal_plots_diffhsed.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(output_dir, "coastal_plots_diffhsed.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close(fig)

# Panels 2–4: maximum concentration per sediment fraction
listvar = ["MUD1", "SAND", "GRAV"]

for var in listvar:
    mud = data[var]
    cmaxmud = mud.max(axis=0).max(axis=0)
    cmaxmudmask = cmaxmud.where(bathy > -5.0)
    cmax_mud = xarray.DataArray(
        cmaxmudmask[:, :], coords=[lat, lon], dims=["lat", "lon"]
    )

    fig = plt.figure(figsize=(8, 4))
    ax = plt.axes(projection=ccrs.PlateCarree())
    p = cmax_mud.plot.pcolormesh(
        transform=ccrs.PlateCarree(),
        cmap=cmap_lev,
        add_colorbar=False,
        levels=lev,
        extend="max",
    )
    cb = matplotlib.pyplot.colorbar(
        p, ax=ax, format=matplotlib.ticker.FuncFormatter(fmt)
    )
    cb.ax.tick_params(labelsize=10)
    cb.set_ticks(lev)
    cb.set_label("kg/m³", size=10)
    ax.set_title(
        f"COASTAL test case — maximum concentration of {var}",
        fontsize=13,
        fontweight="bold",
    )
    add_gridlines(ax)

    # ── Save / Show (do not modify) ──────────────────────────

    if args.makepdf:
        pdf_path = os.path.join(output_dir, f"coastal_plots_cmax{var}.pdf")
        plt.savefig(pdf_path, bbox_inches="tight", transparent=True)
        print(f"PDF file '{pdf_path}' has been created.")

    if args.makepng:
        png_path = os.path.join(output_dir, f"coastal_plots_cmax{var}.png")
        plt.savefig(png_path, bbox_inches="tight", dpi=300)
        print(f"PNG file '{png_path}' has been created.")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)
