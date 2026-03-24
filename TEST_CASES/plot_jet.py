#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the JET test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="jet_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex", type=int, default=19, help="Time index to plot (default: 19)"
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument(
    "--makepng", action="store_true", help="Save individual plots as PNG files"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Constants from Matlab script
R0 = 30
TCOEF = 0.28

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

print("Reading grid and variables...")

# Read coordinates (convert to km)
xr = 1e-3 * nc.variables["x_rho"][:, :]
yr = 1e-3 * nc.variables["y_rho"][:, :]
h = nc.variables["h"][:, :]
pm = nc.variables["pm"][:, :]
pn = nc.variables["pn"][:, :]
f = nc.variables["f"][:, :]

# Dimensions
N = len(nc.dimensions["s_rho"])
tlen = len(nc.variables["scrum_time"][:])
tindex = min(tlen - 1, args.tindex)  # Ensure valid index (0-based)

# Time
time = round(nc.variables["scrum_time"][tindex] / (24 * 3600))
print(f"Day: {time}  index: {tindex}")

# Read variables at specified time
zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
u = np.squeeze(nc.variables["u"][tindex, :, :, :])
v = np.squeeze(nc.variables["v"][tindex, :, :, :])
rho = np.squeeze(nc.variables["temp"][tindex, :, :, :])

# Read initial conditions
u0 = np.squeeze(nc.variables["u"][0, :, :, :])
rho0 = np.squeeze(nc.variables["temp"][0, :, :, :])

print("Processing data...")

# Convert grids using croco_utils functions
xu, xv, xp = cr.rho2uvp(xr)
yu, yv, yp = cr.rho2uvp(yr)
fu, fv, fp = cr.rho2uvp(f)

# Convert rho to temperature
t = (-rho + R0) / TCOEF
t0 = (-rho0 + R0) / TCOEF
sst = np.squeeze(t[N - 1, :, :])  # Surface temperature (Python: N-1)

# Vertical grid
zr = cr.zlevs(h, zeta, 5, 0, 100, N, "r", 2)

# Surface velocities
us = np.squeeze(u[N - 1, :, :])  # Python: N-1
vs = np.squeeze(v[N - 1, :, :])

# Debug dimensions
print(f"pm shape: {pm.shape}")
print(f"pn shape: {pn.shape}")
print(f"us shape: {us.shape}")
print(f"vs shape: {vs.shape}")

# Vorticity
vort = cr.vorticity(us, vs, pm, pn)
vort = vort / fp

# Interpolate velocities to rho grid
ur = cr.u2rho_2d(np.squeeze(u[N - 1, :, :]))
vr = cr.v2rho_2d(np.squeeze(v[N - 1, :, :]))

# Create 3D coordinate grids for vertical sections
XR = cr.tridim(xr, N)
YR = cr.tridim(yr, N)
YU = 0.5 * (YR[:, :, :-1] + YR[:, :, 1:])
zu = 0.5 * (zr[:, :, :-1] + zr[:, :, 1:])

# Middle section for vertical plots
M, L = yr.shape
imid = round(L / 2)
yrs = np.squeeze(YR[:, :, imid])
zrs = np.squeeze(zr[:, :, imid])
yus = np.squeeze(YU[:, :, imid])
zus = np.squeeze(zu[:, :, imid])
u0s = np.squeeze(u0[:, :, imid])
t0s = np.squeeze(t0[:, :, imid])

print("Creating plots...")

# ====================================================================
# Initial fields plot
# ====================================================================
fig1 = plt.figure(figsize=(10, 8))

# Initial Temperature
plt.subplot(2, 1, 1)
temp_contour_levels = np.arange(8, 20, 0.5)  # Contours: 8 to 18 by steps of 2
temp_fill_levels = np.linspace(8, 20, 100)  # Smooth colorbar: many levels
cs1 = plt.contourf(yrs, zrs, t0s, levels=temp_fill_levels, cmap="jet")
plt.contour(yrs, zrs, t0s, levels=temp_contour_levels, colors="black", linewidths=0.5)
plt.clim(8, 20)
plt.ylabel("Z [m]", fontsize=15)
v1 = np.arange(8, 22, 2)
plt.colorbar(cs1, ticks=v1)
plt.title("Initial Temperature", fontsize=15)
plt.gca().tick_params(labelsize=15)

# Initial Zonal Velocity
plt.subplot(2, 1, 2)
vel_contour_levels = np.arange(-0.1, 0.4, 0.01)  # Contours: -0.1 to 0.3 by steps of 0.1
vel_fill_levels = np.linspace(-0.1, 0.3, 100)  # Smooth colorbar: many levels
cs2 = plt.contourf(yus, zus, u0s, levels=vel_fill_levels, cmap="jet")
plt.contour(yus, zus, u0s, levels=vel_contour_levels, colors="black", linewidths=0.5)
plt.clim(-0.1, 0.3)
plt.xlabel("Y [km]", fontsize=15)
plt.ylabel("Z [m]", fontsize=15)
v1 = np.arange(-0.1, 0.4, 0.1)
plt.colorbar(cs2, ticks=v1)
plt.title("Initial Zonal Velocity", fontsize=15)
plt.gca().tick_params(labelsize=15)

plt.tight_layout()

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "jet_init.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "jet_init.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG saved: {png_path}")

if not args.no_show:
    plt.show()

# ====================================================================
# Final dynamic fields plot
# ====================================================================
fig2 = plt.figure(figsize=(15, 6))

# SSH
plt.subplot(1, 3, 1)
ssh_contour_levels = np.arange(-1.2, 1.2, 0.05)  # Contours: -1 to 1 by steps of 0.2
ssh_fill_levels = np.linspace(-1, 1, 100)  # Smooth colorbar: many levels
cs3 = plt.contourf(xr, yr, zeta, levels=ssh_fill_levels, cmap="jet")
plt.contour(xr, yr, zeta, levels=ssh_contour_levels, colors="black", linewidths=0.5)
v1 = np.arange(-1, 1.2, 0.2)
plt.colorbar(cs3, ticks=v1)
plt.clim(-1, 1)
plt.xlim(0, 500)
plt.ylim(0, 2000)
plt.xlabel("X [km]", fontsize=12)
plt.ylabel("Y [km]", fontsize=12)
plt.title(f"SSH - day={time}", fontsize=12)
plt.gca().tick_params(labelsize=12)

# Vorticity/f
plt.subplot(1, 3, 2)
vort_fill_levels = np.linspace(-0.3, 0.3, 100)  # Smooth colorbar: many levels
cs4 = plt.contourf(xp, yp, vort, levels=vort_fill_levels, cmap="jet")
v1 = np.arange(-0.3, 0.4, 0.1)
plt.colorbar(cs4, ticks=v1)
plt.clim(-0.3, 0.3)
plt.xlim(0, 500)
plt.ylim(0, 2000)
plt.xlabel("X [km]", fontsize=12)
plt.title(f"ξ/f - day={time}", fontsize=12)
plt.gca().tick_params(labelsize=12)

# Speed with quiver
plt.subplot(1, 3, 3)
spd = np.sqrt(ur**2 + vr**2)
speed_fill_levels = np.linspace(0, 0.7, 100)  # Smooth colorbar: many levels
cs5 = plt.contourf(xr, yr, spd, levels=speed_fill_levels, cmap="jet")
# More arrows - reduce skip factor even more
skip = 2  # Reduced from 4 to 3 for even more arrows
plt.quiver(
    xr[::skip, ::skip],
    yr[::skip, ::skip],
    ur[::skip, ::skip],
    vr[::skip, ::skip],
    scale=10,
)
plt.clim(0, 0.7)
plt.xlim(0, 500)
plt.ylim(0, 2000)
plt.xlabel("X [km]", fontsize=12)
v1 = np.arange(0, 0.8, 0.1)
plt.colorbar(cs5, ticks=v1)
plt.title(f"Speed - day={time}", fontsize=12)
plt.gca().tick_params(labelsize=12)

plt.tight_layout()

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "jet.pdf")
    plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
    print(f"PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, "jet.png")
    plt.savefig(png_path, bbox_inches="tight", dpi=300)
    print(f"PNG saved: {png_path}")

if not args.no_show:
    plt.show()

# ====================================================================
# Time evolution of SST (if tindex > 17, Python: 0-based)
# ====================================================================
if tindex > 17:
    print("Creating SST time evolution plot...")

    fig3 = plt.figure(figsize=(20, 7))  # Even larger figure

    # Time indices to plot (Python: 0-based, so subtract 1)
    time_indices = [10, 14, 18]  # Matlab [11;15;19] -> Python [10;14;18]

    axes = []
    for i, tidx in enumerate(time_indices):
        if tidx < tlen:  # Check if time index exists
            # Read SST at this time
            rhos = np.squeeze(nc.variables["temp"][tidx, N - 1, :, :])
            sst_t = (-rhos + R0) / TCOEF

            # Create subplot with more spacing
            ax = plt.subplot(1, 3, i + 1)
            axes.append(ax)

            sst_contour_levels = np.arange(
                13, 21, 0.25
            )  # Contours: 13 to 20 by steps of 1
            sst_fill_levels = np.linspace(13, 20, 100)  # Smooth colorbar: many levels
            cs = plt.contourf(xr, yr, sst_t, levels=sst_fill_levels, cmap="jet")
            plt.contour(
                xr, yr, sst_t, levels=sst_contour_levels, colors="black", linewidths=0.5
            )
            plt.clim(13, 20)
            plt.xlabel("X [km]", fontsize=12)
            if i == 0:
                plt.ylabel("Y [km]", fontsize=12)

            time_day = round(nc.variables["scrum_time"][tidx] / (24 * 3600))
            plt.title(f"SST - day={time_day}", fontsize=12)

    # Add much larger colorbar with better positioning
    if axes:
        # Less compression - keep more space between plots
        for i, ax in enumerate(axes):
            pos = ax.get_position()
            # Much less compression: 0.95 instead of 0.92
            ax.set_position([0.95 * pos.x0, pos.y0, 0.9 * pos.width, pos.height])

        # Much larger colorbar
        cbar_ax = fig3.add_axes([0.92, 0.15, 0.015, 0.7])  # Narrow but very tall
        v1 = np.arange(13, 21, 1)
        cbar = plt.colorbar(cs, cax=cbar_ax, ticks=v1)
        cbar.ax.tick_params(labelsize=12)  # Larger tick labels

    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "jet_multivor.pdf")
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        print(f"PDF saved: {pdf_path}")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "jet_multivor.png")
        plt.savefig(png_path, bbox_inches="tight", dpi=300)
        print(f"PNG saved: {png_path}")

    if not args.no_show:
        plt.show()

# Close NetCDF file
nc.close()

print("Script completed successfully!")
