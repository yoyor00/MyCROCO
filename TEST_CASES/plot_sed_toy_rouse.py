#!/usr/bin/env python3
"""
Sed_toy Rouse Test Case

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

Patrick Marchesiello - 2012

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Sed_toy Rouse test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="sed_toy_rouse_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--usgs",
    action="store_true",
    help="Use USGS model (default). If not set, uses MUSTANG model",
)
parser.add_argument(
    "--mustang", action="store_true", help="Use MUSTANG model instead of USGS"
)
parser.add_argument(
    "--nt1", type=int, default=12, help="Time index to use (default: 12)"
)
parser.add_argument("--idy", type=int, default=4, help="Y grid index (default: 4)")
parser.add_argument("--idx", type=int, default=4, help="X grid index (default: 4)")
parser.add_argument(
    "--dt", type=float, default=1.0, help="Time step for display (default: 1.0 s)"
)
parser.add_argument(
    "--nz",
    type=int,
    default=100,
    help="Number of vertical levels for display (default: 100)",
)
parser.add_argument("--makepdf", action="store_true", help="Save plots as PDF files")
parser.add_argument("--makepng", action="store_true", help="Save plots as PNG files")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

print("=== Sed_toy Rouse Test Case ===")
print("Comparing sediment concentration profiles with Rouse analytical solution")

# Determine model type
if args.mustang:
    usgs = False
    model_name = "Mustang"
    tau_var = "TAUSKIN"
    print("Using MUSTANG model")
else:
    usgs = True  # Default
    model_name = "Usgs"
    tau_var = "bostr"
    print("Using USGS model")

if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
fname = args.file
dt = args.dt
nz = args.nz
it = 1  # Fixed at 1 like MATLAB
nt1 = args.nt1
idy = args.idy - 1  # Convert to Python 0-based indexing
idx = args.idx - 1  # Convert to Python 0-based indexing

print(f"Using file: {fname}")
print(f"Grid point: idy={idy + 1}, idx={idx + 1} (Matlab indexing)")
print(f"Time index: nt1={nt1}")
print(f"Model parameters: dt={dt}s, nz={nz}")

# Read data
print("Reading data from CROCO file...")
try:
    with Dataset(fname, "r") as nc:
        # Grid and bathymetry
        h = np.squeeze(nc.variables["h"][:, :])
        x = np.squeeze(nc.variables["x_rho"][:, :])
        y = np.squeeze(nc.variables["y_rho"][:, :])
        zeta = np.squeeze(nc.variables["zeta"][nt1, idy, idx])

        # Bottom stress
        if tau_var in nc.variables:
            tauskin = np.squeeze(nc.variables[tau_var][nt1, idy, idx])
            print(f"Using bottom stress variable: {tau_var}")
        else:
            print(f"Warning: Variable '{tau_var}' not found, trying alternative...")
            if "bostr" in nc.variables:
                tauskin = np.squeeze(nc.variables["bostr"][nt1, idy, idx])
                tau_var = "bostr"
            elif "TAUSKIN" in nc.variables:
                tauskin = np.squeeze(nc.variables["TAUSKIN"][nt1, idy, idx])
                tau_var = "TAUSKIN"
            else:
                raise KeyError("Neither 'bostr' nor 'TAUSKIN' found in file")
            print(f"Using bottom stress variable: {tau_var}")

        # Velocity
        u = np.squeeze(nc.variables["u"][nt1, :, idy, idx])
        N = u.shape[0]

        # Vertical grid parameters
        theta_s = nc.theta_s
        theta_b = nc.theta_b
        hc = nc.hc
        try:
            Vtrans = nc.variables["Vtransform"][:]
        except KeyError:
            Vtrans = 1  # Default value

        # Convert masked arrays if needed
        if hasattr(h, "mask"):
            h = np.asarray(h)
        if hasattr(tauskin, "mask"):
            tauskin = np.asarray(tauskin)
        if hasattr(u, "mask"):
            u = np.asarray(u)

        print(f"Grid dimensions: {h.shape}")
        print(f"Vertical levels: N={N}")
        print(f"Bottom stress: {tauskin:.6f}")

except FileNotFoundError:
    print(f"Error: File '{fname}' not found.")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print("Computing vertical coordinates...")
# Compute vertical coordinates
depth = h[idy, idx]
print(f"Depth at selected point: {depth:.3f} m")

# For zlevs, we need to provide arrays, even if they're 1x1
depth_array = np.array([[depth]])
zeta_array = np.array([[zeta]])

zr = cr.zlevs(depth_array, zeta_array, theta_s, theta_b, hc, N, "r", Vtrans)
zr = np.squeeze(zr)
zu = 0.5 * (zr[:-1] + zr[1:])

h_total = depth + zeta
z = -zr

print("Setting up sediment parameters...")
# Sediment parameters
ws = np.array([0.001, 0.01, 0.02, 0.04, 0.08, 0.1])  # settling velocity (m/s)

if usgs:
    sandstr = ["mud_01", "mud_02", "mud_03", "mud_04", "mud_05", "mud_06"]
else:
    sandstr = ["SED1", "SED2", "SED3", "SED4", "SED5", "SED6"]

# Physical constants
rho = 1030.0  # density
vk = 0.41  # von Karman constant
vk_hydro = 0.41
z0_val = 0.0001
nsand = 6
k_aref = 0  # Python 0-based indexing (k=1 in MATLAB becomes k=0)
a = z[0]  # reference level

print(f"Reference level a = {a:.4f} m")
print(f"Total depth h = {h_total:.4f} m")

# Compute Rouse number
print("Computing Rouse numbers...")
uet = np.sqrt(tauskin / rho)
rouse = ws / (vk * uet)

print("Rouse numbers for each sediment class:")
for i, (w, r) in enumerate(zip(ws, rouse)):
    print(f"  Class {i + 1}: ws={w:.3f} m/s, Rouse={r:.3f}")

# Read sediment concentrations and compute analytical solution
print("Reading sediment concentrations and computing analytical profiles...")
c_ana = np.zeros((nsand, nz))
c_sand = np.zeros((nsand, nz))

try:
    with Dataset(fname, "r") as nc:
        for isand in range(nsand):
            sd = sandstr[isand]
            if sd in nc.variables:
                sandn = np.squeeze(nc.variables[sd][nt1, :, idy, idx])

                # Convert masked array if needed
                if hasattr(sandn, "mask"):
                    sandn = np.asarray(sandn)

                c_sand[isand, :] = sandn[:, it - 1] if sandn.ndim > 1 else sandn

                # Compute analytical Rouse profile
                for k in range(nz):
                    if k < len(z):
                        z_k = z[k]
                        if z_k != a and h_total - z_k > 0 and z_k > 0:
                            # Rouse profile: C(z) = C_ref * ((h-z)/z * a/(h-a))^P
                            ratio = ((h_total - z_k) / z_k) * (a / (h_total - a))
                            if ratio > 0:
                                c_ana[isand, k] = c_sand[isand, k_aref] * (
                                    ratio ** rouse[isand]
                                )
                            else:
                                c_ana[isand, k] = c_sand[isand, k_aref]
                        else:
                            c_ana[isand, k] = c_sand[isand, k_aref]
                    else:
                        c_ana[isand, k] = c_sand[isand, k_aref]

                print(f"  Processed sediment class: {sd}")
            else:
                print(f"  Warning: Sediment variable '{sd}' not found in file")

except Exception as e:
    print(f"Error reading sediment data: {e}")
    exit(1)

print("Normalizing concentrations...")
# Normalize concentrations
rapC_ana = np.zeros((nsand, nz))
rapC_sand = np.zeros((nsand, nz))

for isand in range(nsand):
    # Normalize by first level (bottom)
    if c_sand[isand, 0] != 0:
        rapC_sand[isand, :] = c_sand[isand, :] / c_sand[isand, 0]
    else:
        rapC_sand[isand, :] = c_sand[isand, :]

    # Normalize analytical (flip to match MATLAB indexing)
    if c_ana[isand, -1] != 0:  # Use last element (surface)
        rapC_ana[isand, :] = np.flip(c_ana[isand, :] / c_ana[isand, -1])
    else:
        rapC_ana[isand, :] = np.flip(c_ana[isand, :])

print("Creating vertical profile plots...")
# Plot vertical profiles
fig = plt.figure(figsize=(21, 7))  # 1500x500 equivalent

title0 = f"Vertical profile of sand concentrations at equilibrium / Rouse profile (model:dt={dt}s,nz={nz})"
fig.suptitle(title0, fontsize=14)

for isand in range(nsand):
    ax = plt.subplot(1, nsand, isand + 1)

    # Plot model data (green)
    plt.plot(np.flip(rapC_sand[isand, :]), z, "g-", linewidth=2, label=model_name)

    # Plot analytical Rouse profile (dashed)
    plt.plot(np.flip(rapC_ana[isand, :]), z, "--", linewidth=2, label="Rouse")

    plt.grid(True)
    plt.xlabel("C/C[k=1]", fontsize=12)

    if isand == 0:
        plt.ylabel("Water level [m]", fontsize=12)

    plt.title(f"WS={ws[isand]:.3f} m/s", fontsize=12)
    plt.xlim([0, 1])

    if isand == nsand - 1:  # Add legend to last subplot
        plt.legend()

plt.tight_layout()

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "sed_toy.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "sed_toy.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
print(f"Compared {nsand} sediment classes with Rouse analytical profiles")

# Print summary statistics
print("\n=== Rouse Profile Analysis ===")
print("Sediment classes and their characteristics:")
for i in range(nsand):
    print(f"Class {i + 1}: ws={ws[i]:.3f} m/s, Rouse number={rouse[i]:.3f}")

print(f"\nShear velocity u*: {uet:.6f} m/s")
print(f"Bottom stress: {tauskin:.6f} N/m²")
