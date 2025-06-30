#!/usr/bin/env python3
"""
Sed toy resuspension Test Case

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

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the Sed_toy Resuspension test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="sed_toy_resusp_his.nc", help="Path to the NetCDF file"
)
parser.add_argument("--idy", type=int, default=1, help="Y grid index (default: 1)")
parser.add_argument("--idx", type=int, default=1, help="X grid index (default: 1)")
parser.add_argument(
    "--depth",
    type=float,
    default=20.0,
    help="Water column depth in meters (default: 20.0)",
)
parser.add_argument(
    "--beddepth",
    type=float,
    default=0.041,
    help="Sediment bed depth in meters (default: 0.041)",
)
parser.add_argument(
    "--it",
    type=int,
    default=61,
    help="Time index for stratigraphy (default: 61 - 5 days)",
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

print("=== Sed_toy Resuspension Test Case ===")
print("Analyzing erosion-deposition events and bed stratigraphy")

if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# Parameters
fname = args.file
idy = args.idy - 1  # Convert to Python 0-based indexing
idx = args.idx - 1  # Convert to Python 0-based indexing
depth = args.depth
beddepth = args.beddepth
it = args.it - 1  # Convert to Python 0-based indexing

print(f"Using file: {fname}")
print(f"Grid point: idy={idy + 1}, idx={idx + 1} (Matlab indexing)")
print(f"Water depth: {depth} m")
print(f"Bed depth: {beddepth} m ({beddepth * 100:.1f} cm)")
print(f"Stratigraphy time index: {it + 1} (Matlab indexing)")

# Read data
print("Reading data from CROCO file...")
try:
    with Dataset(fname, "r") as nc:
        # Basic variables
        h = np.squeeze(nc.variables["h"][idy, idx])
        T = np.squeeze(nc.variables["scrum_time"][:]) / 86400  # Convert to days
        nt = len(T)

        # Bottom stress
        bostr = np.squeeze(nc.variables["bostr"][:, idy, idx])

        # Active layer thickness
        ALT = np.squeeze(nc.variables["act_thick"][:, idy, idx])

        # Sediment concentrations in water column
        s1 = np.squeeze(nc.variables["sand_01"][:, idy, idx])
        s2 = np.squeeze(nc.variables["sand_02"][:, idy, idx])
        m1 = np.squeeze(nc.variables["mud_01"][:, idy, idx])
        m2 = np.squeeze(nc.variables["mud_02"][:, idy, idx])

        # Convert masked arrays if needed and ensure 1D
        for var_name, var in [("s1", s1), ("s2", s2), ("m1", m1), ("m2", m2)]:
            if hasattr(var, "mask"):
                var = np.asarray(var)
            # Ensure 1D array
            if var.ndim > 1:
                var = var.flatten()
            # Update the variable
            if var_name == "s1":
                s1 = var
            elif var_name == "s2":
                s2 = var
            elif var_name == "m1":
                m1 = var
            elif var_name == "m2":
                m2 = var

        print(
            f"Sediment array shapes: s1={s1.shape}, s2={s2.shape}, m1={m1.shape}, m2={m2.shape}"
        )

        # Convert concentrations to mass (Kg/m²)
        sand1 = s1 * depth
        sand2 = s2 * depth
        mud1 = m1 * depth
        mud2 = m2 * depth

        # Ensure valid time index
        it = min(it, nt - 1)
        print(f"Using time index: {it + 1}, corresponding to {T[it]:.1f} days")

        # Bed fractions for stratigraphy
        fracs1 = np.squeeze(nc.variables["bed_frac_sand_01"][it, :, idy, idx])
        fracs2 = np.squeeze(nc.variables["bed_frac_sand_02"][it, :, idy, idx])
        fracm1 = np.squeeze(nc.variables["bed_frac_mud_01"][it, :, idy, idx])
        fracm2 = np.squeeze(nc.variables["bed_frac_mud_02"][it, :, idy, idx])

        # Bed thickness
        bedthick = np.squeeze(nc.variables["bed_thick"][it, :, idy, idx])
        NL = len(bedthick)

        # Convert remaining masked arrays if needed
        variables_to_check = [
            ("T", T),
            ("bostr", bostr),
            ("ALT", ALT),
            ("fracs1", fracs1),
            ("fracs2", fracs2),
            ("fracm1", fracm1),
            ("fracm2", fracm2),
            ("bedthick", bedthick),
        ]

        for var_name, var in variables_to_check:
            if hasattr(var, "mask"):
                var = np.asarray(var)
                # Update the variable in the local scope
                if var_name == "T":
                    T = var
                elif var_name == "bostr":
                    bostr = var
                elif var_name == "ALT":
                    ALT = var
                elif var_name == "fracs1":
                    fracs1 = var
                elif var_name == "fracs2":
                    fracs2 = var
                elif var_name == "fracm1":
                    fracm1 = var
                elif var_name == "fracm2":
                    fracm2 = var
                elif var_name == "bedthick":
                    bedthick = var

        print(f"Time series length: {nt} steps")
        print(f"Bed layers: {NL}")
        print(f"Time range: {T[0]:.2f} to {T[-1]:.2f} days")

except FileNotFoundError:
    print(f"Error: File '{fname}' not found.")
    exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    exit(1)

print("Processing stratigraphy data...")
# Compute bed levels
# MATLAB: zbed=[-beddepth;-beddepth+cumsum(flipud(bedthick(2:end,:)),1)];
# Start from bottom (-beddepth) and add cumulative thickness upward

print(f"Bed layers: {NL}, bedthick shape: {bedthick.shape}")

# Create zbed array with correct size
zbed = np.zeros(NL)
zbed[0] = -beddepth

# Add cumulative thickness (skip first layer like MATLAB)
if NL > 1:
    cumthick = np.cumsum(np.flip(bedthick[1:]))  # Skip first element, flip, then cumsum
    zbed[1 : len(cumthick) + 1] = -beddepth + cumthick

# Convert to cm and flip to match MATLAB orientation
zbed = np.flip(zbed) * 100

print("Organizing sediment data...")
# Organize concentration data for plotting
Yconc = np.zeros((4, nt))  # 4: 2 sands + 2 muds
for i in range(nt):
    Yconc[:, i] = [sand1[i], sand2[i], mud1[i], mud2[i]]

# Organize fraction data for stratigraphy
Yfrac = np.zeros((4, NL))  # 4: 2 sands + 2 muds
for k in range(NL):
    Yfrac[:, k] = [fracs1[k], fracs2[k], fracm1[k], fracm2[k]]

print("Creating plots...")
# Create figure
fig = plt.figure(figsize=(15, 11))  # 1050x800 equivalent
fig.suptitle("Two successive erosion–deposition events lasting 5 days", fontsize=14)

# Sediment labels and colors
sediment_labels = ["Sand1", "Sand2", "Mud1", "Mud2"]
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]  # Default matplotlib colors

# Subplot 1: Mass of sediment in suspension
ax1 = plt.subplot(3, 2, 1)
plt.stackplot(
    T,
    Yconc[0],
    Yconc[1],
    Yconc[2],
    Yconc[3],
    labels=sediment_labels,
    colors=colors,
    alpha=0.8,
)
plt.ylabel("Mass of sediment (Kg/m²)", fontsize=12)
plt.legend(loc="upper right")
plt.grid(True)
plt.tick_params(labelsize=10)

# Subplot 2: Bottom stress
ax3 = plt.subplot(3, 2, 3)
plt.plot(T, bostr, "b-", linewidth=1)
plt.ylabel("Bottom stress (Pa)", fontsize=12)
plt.grid(True)
plt.tick_params(labelsize=10)

# Subplot 3: Active Layer Thickness
ax5 = plt.subplot(3, 2, 5)
plt.plot(T, ALT * 100, "g-", linewidth=1)  # Convert to cm
plt.ylabel("Active Layer thickness (cm)", fontsize=12)
plt.xlabel("Time (days)", fontsize=12)
plt.grid(True)
plt.tick_params(labelsize=10)

# Subplot 4: Stratigraphy (spanning right side)
ax_strat = plt.subplot(3, 2, (2, 6))

# Create horizontal stacked area plot for stratigraphy
# We need to rotate this to be vertical like MATLAB's camroll(-270)
bottom = np.zeros(len(zbed))
for i in range(4):
    if i < len(Yfrac):
        # Interpolate fractions to bed levels if needed
        if len(Yfrac[i]) == len(zbed):
            values = Yfrac[i]
        else:
            # Interpolate to match zbed length
            from scipy.interpolate import interp1d

            f = interp1d(
                np.linspace(0, 1, len(Yfrac[i])),
                Yfrac[i],
                kind="linear",
                bounds_error=False,
                fill_value=0,
            )
            values = f(np.linspace(0, 1, len(zbed)))

        plt.barh(
            zbed,
            values,
            left=bottom,
            height=np.diff(zbed, prepend=zbed[0] - 0.1),
            color=colors[i],
            label=sediment_labels[i],
            alpha=0.8,
        )
        bottom += values

# Add horizontal lines for layer boundaries
for z in zbed:
    plt.axhline(y=z, color="black", linestyle="--", alpha=0.5, linewidth=0.8)

plt.xlim([0, 1])
plt.ylim([-4.1, 0])
plt.xlabel("Fraction of sediment", fontsize=12)
plt.ylabel(f"Z (cm) - Day {T[it]:.0f}", fontsize=12)
plt.legend(loc="upper right")
plt.grid(True)
plt.tick_params(labelsize=10)

# Adjust layout
plt.tight_layout()

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "sed_toy_resusp.pdf")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "sed_toy_resusp.png")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

# Show or close plot
if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")

# Print summary statistics
print("\n=== Resuspension Analysis Summary ===")
print(f"Maximum bottom stress: {np.max(bostr):.3f} Pa")
print(f"Maximum active layer thickness: {np.max(ALT) * 100:.2f} cm")
print(f"Total sediment in suspension (final): {np.sum(Yconc[:, -1]):.3f} Kg/m²")

print("\nSediment mass evolution:")
for i, label in enumerate(sediment_labels):
    initial = Yconc[i, 0]
    final = Yconc[i, -1]
    max_val = np.max(Yconc[i, :])
    print(
        f"  {label}: Initial={initial:.4f}, Final={final:.4f}, Max={max_val:.4f} Kg/m²"
    )

print(f"\nBed stratigraphy at day {T[it]:.0f}:")
for i, label in enumerate(sediment_labels):
    mean_frac = np.mean(Yfrac[i, :])
    print(f"  {label}: Mean fraction={mean_frac:.3f}")
