#!/usr/bin/env python3
"""
Moving Bathy Test Case

Internal gravity waves are produced over an oscillating ridge
and Brunt-Vaisala frequency anomaly is compared with lab experiment
(Fig. 1 in Auclair et al., 2014)

Reference:
----------
Auclair et al., 2014: Implementation of a time-dependent bathymetry
 in a free-surface ocean model: Application to internal wave
 generation, Ocean Modelling, 80, 1-9.

Further Information:
http://www.croco-ocean.org

This file is part of CROCOTOOLS

CROCOTOOLS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

CROCOTOOLS is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston,
MA  02111-1307  USA

Patrick Marchesiello - 2020

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Moving Bathy test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="movbat_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex",
    type=int,
    default=49,
    help="Time index to plot (default: 49 - after 10 periods)",
)
parser.add_argument(
    "--j", type=int, default=2, help="Y index for cross-section (default: 2)"
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument("--makepng", action="store_true", help="Save plots as PNG files")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
parser.add_argument(
    "--additional-plots",
    action="store_true",
    help="Generate additional diagnostic plots",
)
parser.add_argument(
    "--passive-tracer", action="store_true", help="Plot passive tracer if available"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

print("=== Moving Bathy Test Case ===")
print("Internal gravity waves over oscillating ridge")
if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
fname = args.file
tndx = args.tindex
j = args.j - 1  # Convert to Python 0-based indexing
g = 9.8  # gravity acceleration

print(f"Using file: {fname}")
print(f"Y-section index: {j + 1} (Matlab indexing)")

# Read data
print("Reading data from CROCO file...")
try:
    with Dataset(fname, "r") as nc:
        # Check time index
        time = nc.variables["scrum_time"][:]
        tndx = min(tndx, len(time) - 1)  # Convert to 0-based and check bounds
        print(f"tndx = {tndx}")

        # Grid and bathymetry (time-dependent)
        h = np.squeeze(nc.variables["hmorph"][tndx, j, :])
        x = np.squeeze(nc.variables["x_rho"][j, :])
        zeta = np.squeeze(nc.variables["zeta"][tndx, j, :])
        L = len(zeta)

        # Density fields
        r = np.squeeze(nc.variables["rho"][tndx, :, j, :])  # Current density
        r0 = np.squeeze(nc.variables["rho"][0, :, j, :])  # Initial density

        # Velocity fields
        u = np.squeeze(nc.variables["u"][tndx, :, j, :])
        w = np.squeeze(nc.variables["w"][tndx, :, j, :])

        # Vertical grid parameters
        theta_s = nc.theta_s
        theta_b = nc.theta_b
        rho0 = nc.rho0
        hc = nc.hc
        try:
            Vtrans = nc.variables["Vtransform"][:]
        except KeyError:
            Vtrans = 1  # Default value

        N = len(nc.dimensions["s_rho"])

        print(f"Grid dimensions: N={N}, L={L}")
        print(f"Time: {time[tndx]:.2f} s")

except FileNotFoundError:
    print(f"Error: File '{fname}' not found.")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print("Computing vertical coordinates...")
# Compute vertical coordinates
zr = np.squeeze(cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", Vtrans))
zw = np.squeeze(cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "w", Vtrans))

# Create u-grid coordinates
zu = np.zeros((N, L - 1))
zu[:, :] = 0.5 * (zr[:, :-1] + zr[:, 1:])

# Create horizontal grids
xr = np.tile(x, (N, 1))
xw = np.tile(x, (N + 1, 1))
xu = np.zeros((N, L - 1))
xu[:, :] = 0.5 * (xr[:, :-1] + xr[:, 1:])

print("Computing Brunt-Vaisala frequency anomaly...")
# Compute Brunt-Vaisala frequency anomaly
# bvf = -g/rho0 * d(rho-rho_ini)/dz
r_anom = r - r0
bvf = np.zeros_like(zw)

# Compute gradient (avoiding boundaries)
bvf[1:-1, :] = -g / rho0 * (r_anom[1:, :] - r_anom[:-1, :]) / (zr[1:, :] - zr[:-1, :])

# Fill boundaries
bvf[0, :] = bvf[1, :]
bvf[-1, :] = bvf[-2, :]

# Remove mean
bvf = bvf - np.mean(bvf)

print("Creating main plot - Brunt-Vaisala frequency...")
# Main plot - Brunt-Vaisala frequency
fig = plt.figure(figsize=(8.5, 5))  # 600x350 equivalent

# Create 80 levels like MATLAB
levels = np.linspace(-0.05, 0.05, 81)
cs = plt.contourf(xw, zw, bvf, levels=levels, cmap="gray")

# Plot bathymetry line
plt.plot(xr[0, :], -h, "k-", linewidth=3)

# Colorbar with exact MATLAB range
cbar = plt.colorbar(cs)
cbar.set_ticks(np.arange(-0.05, 0.06, 0.025))  # Ticks every 0.025
cbar.set_ticklabels([f"{x:.3f}" for x in np.arange(-0.05, 0.06, 0.025)])

plt.clim(-0.05, 0.05)
plt.axis([-0.47, 0.04, -0.4, 0])

# Set exact ticks like MATLAB
plt.xticks(np.arange(-0.4, 0.01, 0.1))
plt.yticks(np.arange(-0.4, 0.01, 0.1))

plt.xlabel("X [m]")
plt.ylabel("Z [m]")
plt.title("Moving Bathy - N² (dx=8mm)")
plt.tick_params(labelsize=15)

# Save main plot
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "MovBat.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "MovBat.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()

# Additional diagnostic plots
if args.additional_plots:
    print("Creating additional diagnostic plots...")

    # Density anomaly plot
    fig2 = plt.figure(figsize=(10, 6))
    levels_rho = np.linspace(np.min(r_anom), np.max(r_anom), 21)
    cs2 = plt.contourf(xr, zr, r_anom, levels=levels_rho, cmap="viridis")
    plt.plot(xr[0, :], -h, "k-", linewidth=2)
    plt.colorbar(cs2, label="Density anomaly")
    plt.axis([-0.5, 0.5, -0.4, 0])
    plt.xlabel("X [m]")
    plt.ylabel("Z [m]")
    plt.title("Moving Bathy - Density Anomaly")

    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "MovBat_rho.pdf")
        plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"PDF file '{pdf_path}' has been created.")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "MovBat_rho.png")
        plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
        print(f"PNG file '{png_path}' has been created.")

    if not args.no_show:
        plt.show()
    else:
        plt.close()

    # U velocity plot
    fig3 = plt.figure(figsize=(10, 6))
    levels_u = np.linspace(np.min(u), np.max(u), 21)
    cs3 = plt.contourf(xu, zu, u, levels=levels_u, cmap="RdBu_r")
    plt.plot(xr[0, :], -h, "k-", linewidth=2)
    plt.colorbar(cs3, label="U velocity [m/s]")
    plt.axis([-0.5, 0.5, -0.4, 0])
    plt.xlabel("X [m]")
    plt.ylabel("Z [m]")
    plt.title("Moving Bathy - U Velocity")

    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "MovBat_u.pdf")
        plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"PDF file '{pdf_path}' has been created.")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "MovBat_u.png")
        plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
        print(f"PNG file '{png_path}' has been created.")

    if not args.no_show:
        plt.show()
    else:
        plt.close()

    # W velocity plot
    fig4 = plt.figure(figsize=(10, 6))
    # Handle case where w might have different dimensions
    if w.shape[0] < zw.shape[0]:
        xw_plot = xr
        zw_plot = zr
        w_plot = w
    else:
        xw_plot = xw
        zw_plot = zw
        w_plot = w

    levels_w = np.linspace(np.min(w_plot), np.max(w_plot), 21)
    cs4 = plt.contourf(xw_plot, zw_plot, w_plot, levels=levels_w, cmap="RdBu_r")
    plt.plot(xr[0, :], -h, "k-", linewidth=2)
    plt.colorbar(cs4, label="W velocity [m/s]")
    plt.axis([-0.5, 0.5, -0.4, 0])
    plt.xlabel("X [m]")
    plt.ylabel("Z [m]")
    plt.title("Moving Bathy - W Velocity")

    if args.makepdf:
        pdf_path = os.path.join(args.output_dir, "MovBat_w.pdf")
        plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"PDF file '{pdf_path}' has been created.")

    if args.makepng:
        png_path = os.path.join(args.output_dir, "MovBat_w.png")
        plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
        print(f"PNG file '{png_path}' has been created.")

    if not args.no_show:
        plt.show()
    else:
        plt.close()

# Passive tracer plot (if requested and available)
if args.passive_tracer:
    print("Checking for passive tracer...")
    try:
        with Dataset(fname, "r") as nc:
            if "tpas" in nc.variables:
                print("Plotting passive tracer...")
                t_pas = np.squeeze(nc.variables["tpas"][tndx, :, j, :])

                fig5 = plt.figure(figsize=(10, 6))
                cs5 = plt.pcolormesh(xr, zr, t_pas, shading="flat", cmap="viridis")
                plt.colorbar(cs5, label="Passive tracer")
                plt.axis([-0.5, 0.5, -0.4, 0])
                plt.xlabel("X [m]")
                plt.ylabel("Z [m]")
                plt.title("Moving Bathy - Passive Tracer")

                if args.makepdf:
                    pdf_path = os.path.join(args.output_dir, "MovBat_tpas.pdf")
                    plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
                    print(f"PDF file '{pdf_path}' has been created.")

                if args.makepng:
                    png_path = os.path.join(args.output_dir, "MovBat_tpas.png")
                    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
                    print(f"PNG file '{png_path}' has been created.")

                if not args.no_show:
                    plt.show()
                else:
                    plt.close()
            else:
                print("Passive tracer 'tpas' not found in file.")
    except Exception as e:
        print(f"Error reading passive tracer: {e}")

print("Script completed successfully!")
