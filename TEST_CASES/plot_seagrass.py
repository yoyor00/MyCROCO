#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr
import xarray


def LOAD_MODEL(filein, i0, j0):
    # get the mean flow over the last 30s of the run (results are outputed every 0.25s)
    data = xarray.open_dataset(filein)
    ptu = data.u.isel(eta_rho=j0, xi_u=i0, time=slice(-31 * 4, -1)).mean(axis=0)
    sig = data.s_rho
    sigw = data.s_w
    Cs = data.Cs_rho
    Csw = data.Cs_w
    hc = data.hc
    h = data.h.isel(eta_rho=j0, xi_rho=i0)
    zeta = data.zeta.isel(eta_rho=j0, xi_rho=i0, time=slice(-31 * 4, -1)).mean(axis=0)
    z = zeta + (zeta + h) * ((hc * sig + h * Cs) / (hc + h)) + h
    return z, ptu


def plot_pt(ax, filemodel, pt, labelpt):
    zmodel, Pmodel = LOAD_MODEL(filemodel, pt, 3)
    ax.plot(Pmodel.values, zmodel.values, ".-", label=labelpt)
    ax.set_xlim([-0.05, 0.5])
    ax.set_ylim([0.0, 0.25])
    ax.legend()


# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot and analyze SEAGRASS test case output from a CROCO NetCDF file.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="seagrass_his.nc",
    help="Path to the NetCDF file (default: seagrass_his.nc)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate PDF files for plots"
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate PNG files for plots"
)
parser.add_argument("--no-show", action="store_true", help="Suppress displaying plots")
parser.add_argument(
    "--output-dir",
    type=str,
    default="./",
    help="Directory to save output files (default: current directory)",
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Parameters
makepdf = args.makepdf
makepng = args.makepng
no_show = args.no_show
nc_file = args.file

# Read NetCDF file
try:
    nc = Dataset(nc_file, "r")
except FileNotFoundError:
    print(f"Error: File '{nc_file}' not found.")
    exit(1)


dx = 0.05
ptlist = [11, 17, 23, 29]
fig, [ax1, ax2, ax3, ax4] = plt.subplots(
    nrows=1, ncols=4, figsize=(15, 4), constrained_layout=True
)
fig.suptitle("seagrass")
plot_pt(ax1, nc_file, ptlist[0], "P1")
plot_pt(ax2, nc_file, ptlist[1], "P2")
plot_pt(ax3, nc_file, ptlist[2], "P3")
plot_pt(ax4, nc_file, ptlist[3], "P4")

# Save outputs
if args.makepng:
    path = os.path.join(args.output_dir, "profile_seagrass.png")
    fig.savefig(path)
    print(f"PNG saved: {path}")

if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "profile_seagrass.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved: {pdf_path}")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()
