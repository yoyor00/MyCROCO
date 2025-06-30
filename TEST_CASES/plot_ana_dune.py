#!/usr/bin/env python3
"""
Analytical Dune Test Case

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
import argparse
import os
from scipy.interpolate import interp1d

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Analytical Dune test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="ana_dune_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--tindex",
    type=int,
    default=2,
    help="Time index to plot (default: 2 - half-hourly outputs)",
)
parser.add_argument(
    "--j", type=int, default=2, help="Y index for cross-section (default: 2)"
)
parser.add_argument(
    "--usgs",
    action="store_true",
    help="Use USGS model (default). If not set, uses MUSTANG model",
)
parser.add_argument(
    "--mustang", action="store_true", help="Use MUSTANG model instead of USGS"
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

print("=== Analytical Dune Test Case ===")
print("Comparing analytical and numerical bed evolution")

# Determine model type
if args.mustang:
    usgs = False
    model_name = "MUSTANG"
    var_hmorph = "Hm"
    print("Using MUSTANG model")
else:
    usgs = True  # Default
    model_name = "USGS"
    var_hmorph = "hmorph"
    print("Using USGS model")

if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
fname = args.file
j = args.j - 1  # Convert to Python 0-based indexing
tndx = args.tindex

print(f"Using file: {fname}")
print(f"Y-section index: {j + 1} (Matlab indexing)")

# Read Model data
print("Reading model data...")
try:
    with Dataset(fname, "r") as nc:
        t0 = nc.variables["scrum_time"][:]

        # Compute relative time
        nt = len(t0)
        t = np.zeros(nt)
        for i in range(1, nt):
            t[i] = t0[i] - t0[0]

        tndx = min(tndx, len(t) - 1)  # Ensure valid index
        print(f"tndx = {tndx} - Time = {t[tndx] / 3600:.1f} hr")

        x = np.squeeze(nc.variables["x_rho"][j, :])

        # Convert from masked array to regular array if needed
        if hasattr(x, "mask"):
            x = np.asarray(x)

        # Read morphological height
        if var_hmorph in nc.variables:
            hm = np.squeeze(nc.variables[var_hmorph][:, j, :])
        else:
            print(f"Warning: Variable '{var_hmorph}' not found, trying alternative...")
            # Try alternative names
            if "hmorph" in nc.variables:
                hm = np.squeeze(nc.variables["hmorph"][:, j, :])
                var_hmorph = "hmorph"
            elif "Hm" in nc.variables:
                hm = np.squeeze(nc.variables["Hm"][:, j, :])
                var_hmorph = "Hm"
            else:
                raise KeyError(f"Neither 'hmorph' nor 'Hm' found in file")

        # Convert from masked array to regular array if needed
        if hasattr(hm, "mask"):
            hm = np.asarray(hm)

        print(f"Using variable: {var_hmorph}")
        print(f"Data shape: t={len(t)}, x={len(x)}")

except FileNotFoundError:
    print(f"Error: File '{fname}' not found.")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print("Computing analytical solution...")
# Analytical solution
# From Marrieu PhD thesis 2007 and
# Long, Wen, James T. Kirby, and Zhiyu Shao, A numerical scheme for
# morphological bed level calculations.
# Coastal Engineering, 2007, doi:10.1016/j.coastaleng.2007.09.009.

# Model parameters
alpha = 0.001
gamma = 0.01
beta = 3
h0 = 6
xc = 150
Q = 10
poros = 0.4

# Initial bathymetry
h = h0 - 2 * np.exp(-gamma * (x - xc) ** 2)

# Initial phase speed of bedform
# dq/dh*1/(1-poros), with bedload q=alpha*u^beta
C = alpha * beta * Q**beta / (-h) ** (beta + 1) / (1 - poros)

# Set array of final bathymetry
ha = np.full((len(t), len(x)), np.nan)
ha[0, :] = h

# Mass budget under bedform
sed_dune = np.sum(h0 - h)
print(f"Initial sediment mass under dune: {sed_dune:.2f}")

# Compute analytical evolution
print("Computing bed evolution...")
for k in range(1, len(t)):
    if k % 10 == 0:  # Progress indicator
        print(f"  Processing time step {k}/{len(t) - 1}")

    # Find location of bedform: hai and hbi
    # hai and hbi are the 2 possible solutions
    # At any time, a location can originate from a position
    # located both upstream and downstream of the bedform

    hai = np.full(len(x), np.nan)
    hbi = np.full(len(x), np.nan)
    d = np.zeros(len(x))

    # Position of first half of bedform at t
    mid_point = len(x) // 2
    for i in range(mid_point):
        d[i] = x[i] + C[i] * t[k]
        idx = np.where(x > d[i])[0]
        if len(idx) > 0:
            hai[idx[0]] = h[i]

    # Position of second half of bedform at t
    for i in range(mid_point, len(x)):
        d[i] = x[i] + C[i] * t[k]
        idx = np.where(x > d[i])[0]
        if len(idx) > 0:
            hbi[idx[0]] = h[i]

    # Interpolate bathymetry
    # For hai
    valid_hai = ~np.isnan(hai)
    if np.sum(valid_hai) > 1:
        # Convert masked arrays to regular arrays
        x_interp = np.asarray(x[valid_hai])
        hai_interp = np.asarray(hai[valid_hai])

        f_hai = interp1d(
            x_interp, hai_interp, kind="linear", bounds_error=False, fill_value=np.nan
        )
        hai = f_hai(np.asarray(x))

    # For hbi
    valid_hbi = ~np.isnan(hbi)
    if np.sum(valid_hbi) > 1:
        # Convert masked arrays to regular arrays
        x_interp = np.asarray(x[valid_hbi])
        hbi_interp = np.asarray(hbi[valid_hbi])

        f_hbi = interp1d(
            x_interp, hbi_interp, kind="linear", bounds_error=False, fill_value=np.nan
        )
        hbi = f_hbi(np.asarray(x))

    # Locate shock: where hai and hbi solutions should merge
    # hai and hbi are the solutions before and after the shock

    # Initialize with the first point of the second half of bedform
    first_hbi = np.where(~np.isnan(hbi))[0]
    if len(first_hbi) > 0:
        id_shock = max(0, first_hbi[0] - 1)
    else:
        id_shock = mid_point

    # Look for the solution (hai or hbi) that maintains mass budget
    mass_sed = 0
    mass_sedp1 = 0
    max_iter = len(x) - 1

    while mass_sed <= sed_dune and mass_sed <= mass_sedp1 and id_shock < max_iter:
        id_shock += 1

        # Calculate mass for current shock position
        hai_part = hai[:id_shock]
        hbi_part = hbi[id_shock:]

        mass_sed = np.nansum(h0 - hai_part) + np.nansum(h0 - hbi_part)

        # Calculate mass for next shock position
        if id_shock < len(x) - 1:
            hai_part_p1 = hai[: id_shock + 1]
            hbi_part_p1 = hbi[id_shock + 1 :]
            mass_sedp1 = np.nansum(h0 - hai_part_p1) + np.nansum(h0 - hbi_part_p1)
        else:
            mass_sedp1 = np.nansum(h0 - hai)

    # Final solution combining hai and hbi
    ha[k, :] = np.concatenate([hai[:id_shock], hbi[id_shock:]])

print("Creating plot...")
# Plot analytical and numerical bed evolution
fig = plt.figure(figsize=(13, 8))  # 920x550 equivalent

# Plot every 3rd time step up to half the total time
n_steps = len(t) // 2
step_size = max(1, n_steps // 10)  # Ensure we don't plot too many lines

colors_num = []
colors_ana = []

for n in range(0, n_steps, step_size):
    # Numerical solution (black)
    line_num = plt.plot(x, -hm[n, :], "k-", linewidth=2)
    if n == 0:
        colors_num = line_num

    # Analytical solution (red)
    line_ana = plt.plot(x, -ha[n, :], "r-", linewidth=2)
    if n == 0:
        colors_ana = line_ana

# Add legend
plt.text(120, -3.6, "Numerical", color="k", fontsize=20)
plt.text(120, -3.8, "Analytical", color="r", fontsize=20)

plt.grid(True)
plt.axis([120, 200, -7, -3])
plt.xlabel("X [m]")
plt.ylabel("Bed Level [m]")

if model_name == "MUSTANG":
    plt.title("ANA-DUNE Test Case (MUSTANG) - Bed evolution")
    output_suffix = "mustang"
else:
    plt.title("ANA-DUNE Test Case (USGS) - Bed evolution")
    output_suffix = "usgs"

plt.tick_params(labelsize=15)

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, f"ana_dune_{output_suffix}.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, f"ana_dune_{output_suffix}.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
print(f"Processed {len(t)} time steps with analytical dune migration model")
print(f"Final plot shows evolution from t=0 to t={t[n_steps - 1] / 3600:.1f} hours")
