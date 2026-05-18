#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse
import croco_utils as cr

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the SANDBAR test case.",
    epilog="Example usage:\n"
    + "  python plot_sandbar.py --mycase 1B --makepng\n"
    + "  python plot_sandbar.py --mycase 1C --makepdf --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="sandbar_his.nc", help="Path to the NetCDF file"
)
parser.add_argument(
    "--mycase",
    type=str,
    choices=["1B", "1C"],
    default="1B",
    help="LIP experiment case: '1B' (erosion) or '1C' (accretion)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Save the plot as a PDF file"
)
parser.add_argument(
    "--makepng", action="store_true", help="Save the plot as a PNG file"
)
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

# Open NetCDF file
try:
    nc = Dataset(args.file, "r")
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# =====================================================================
# Input parameters
# =====================================================================
mycase = args.mycase
morph_fac = 18 if mycase == "1B" else 13
morph_cpl = True

# =====================================================================
# Read grid from numerical model
# =====================================================================
yindex = 1  # Python indexing (MATLAB yindex = 2)
tindex = len(nc.variables["scrum_time"][:]) - 1  # Last time record

time = morph_fac * nc.variables["scrum_time"][:] / 3600  # Time in hours
tindex0 = np.argmin(np.abs(time - (4 if mycase == "1B" else 3)))

# Horizontal grid
hr = np.squeeze(nc.variables["h"][yindex, :])
xr = np.squeeze(nc.variables["x_rho"][yindex, :])
hu = 0.5 * (hr[:-1] + hr[1:])
xu = 0.5 * (xr[:-1] + xr[1:])
L = len(hr)

# New bathymetry from last record
if morph_cpl:
    hnew = np.squeeze(nc.variables["hmorph"][tindex, yindex, :])
    h = hnew.copy()
    h0 = np.squeeze(nc.variables["hmorph"][tindex0, yindex, :])
    if hnew is None or np.all(np.isnan(hnew)):
        h = hr.copy()
        hnew = hr.copy()
        h0 = hr.copy()
else:
    h = hr.copy()
    hnew = hr.copy()
    h0 = hr.copy()

# Vertical grid parameters
N = len(nc.dimensions["s_rho"])
theta_s = nc.theta_s
theta_b = nc.theta_b
hc = nc.hc
zeta = np.squeeze(nc.variables["zeta"][tindex, yindex, :])
Dcrit = 1.1 * nc.variables["Dcrit"][:]

# Add land topography
zeta[h < Dcrit] = zeta[h < Dcrit] - h[h < Dcrit]

# Compute vertical levels
zr = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "r", 2)
zw = cr.zlevs(h, zeta, theta_s, theta_b, hc, N, "w", 2)
zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
zwu = 0.5 * (zw[:, :-1] + zw[:, 1:])

# 2D grids for plotting
xu2d = np.tile(xu, (N, 1))
xr2d = np.tile(xr, (N, 1))
D = zw[N, :] - zw[0, :]
Du = zwu[N, :] - zwu[0, :]
Du2d = np.tile(Du, (N, 1))

# =====================================================================
# Read/compute model fields at tindex
# =====================================================================
time_days = (
    morph_fac
    / 86400
    * (nc.variables["scrum_time"][tindex] - nc.variables["scrum_time"][0])
)
thour = int(np.floor(time_days * 24))

# Zonal velocity
u = np.squeeze(nc.variables["u"][tindex, :, yindex, :])

# Sediment concentration
C = 2 * np.squeeze(nc.variables["sand_01"][tindex, :, yindex, :])

# =====================================================================
# Read/compute fields at tindex0 (mid-experiment time)
# =====================================================================
# Hrms
hrms = np.squeeze(nc.variables["hrm"][tindex0, yindex, :])

# U undertow at reference height
zref = 0.1
u0 = np.squeeze(nc.variables["u"][tindex0, :, yindex, :])
hu_0 = 0.5 * (h0[:-1] + h0[1:])
Lu = len(hu_0)
hu2d = np.tile(hu_0, (N, 1))
z0 = np.squeeze(nc.variables["zeta"][tindex0, yindex, :])
z0[h0 < Dcrit] = z0[h0 < Dcrit] - h0[h0 < Dcrit]  # Add land topography
zr0 = cr.zlevs(h0, z0, theta_s, theta_b, hc, N, "r", 2)
zru0 = 0.5 * (zr0[:, :-1] + zr0[:, 1:])
zzu = hu2d + zru0

ubot = np.zeros(Lu)
for ix in range(Lu):
    indices = np.where(zzu[:, ix] > zref)[0]
    if len(indices) > 0:
        nn = indices[0]
    else:
        nn = N - 1
    ubot[ix] = u0[nn, ix]

# Sediment concentration at reference height
C0 = 2 * np.squeeze(nc.variables["sand_01"][tindex0, :, yindex, :])
zref_c = 0.05
h2d = np.tile(h0, (N, 1))
zz = h2d + zr0

Cbot = np.zeros(L)
for ix in range(L):
    indices = np.where(zz[:, ix] > zref_c)[0]
    if len(indices) > 0:
        nn = indices[0]
    else:
        nn = N - 1
    Cbot[ix] = C0[nn, ix]

nc.close()

# Apply masks
zeta_plot = zeta.copy()
zeta_plot[D < Dcrit] = np.nan
u_plot = u.copy()
u_plot[Du2d < Dcrit] = np.nan
hrms_plot = hrms.copy()
hrms_plot[D <= max(0.1, Dcrit)] = np.nan
hrms_plot[hr < 0] = np.nan

# =====================================================================
# FLUME DATA (LIP: Roelvink & Reniers 1995)
# =====================================================================
xd = np.linspace(0, 200, 64)

# LIP-1B bathymetry at 18h on xd grid
hd1B = -np.array(
    [
        -4.0935,
        -4.0931,
        -4.0926,
        -4.0921,
        -4.0917,
        -4.0912,
        -4.0217,
        -3.9788,
        -3.8282,
        -3.6638,
        -3.5145,
        -3.3757,
        -3.2330,
        -3.0776,
        -2.9134,
        -2.7556,
        -2.6207,
        -2.4878,
        -2.3855,
        -2.3182,
        -2.2756,
        -2.2480,
        -2.2272,
        -2.2064,
        -2.1806,
        -2.1467,
        -2.1033,
        -2.0508,
        -1.9910,
        -1.9269,
        -1.8621,
        -1.8006,
        -1.7456,
        -1.6995,
        -1.6625,
        -1.6319,
        -1.6011,
        -1.5583,
        -1.4985,
        -1.4276,
        -1.3233,
        -1.0958,
        -0.8628,
        -0.8772,
        -0.9906,
        -1.0213,
        -1.0502,
        -1.0420,
        -0.9964,
        -0.8379,
        -0.5619,
        -0.4863,
        -0.4529,
        -0.4259,
        -0.3967,
        -0.3527,
        -0.2812,
        -0.1732,
        -0.0251,
        0.1585,
        0.3641,
        0.5696,
        0.7438,
        0.8490,
    ]
)

x_hrms_d1B = np.array([20, 65, 100, 115, 130, 138, 145, 152, 160, 170])
hrms_d1B = np.array([0.85, 0.80, 0.72, 0.65, 0.58, 0.52, 0.39, 0.36, 0.36, 0.25])

x_ubot_d1B = np.array([65, 102, 115, 130, 138, 145, 152, 160, 170])
ubot_d1B = -np.array([12.3, 13.0, 12.9, 17.2, 30.3, 30.4, 17.4, 15.2, 13.6])

x_Cbot_d1B = np.array([65, 102, 115, 130, 138, 145, 152, 170])
Cbot_d1B = np.array([0.40, 0.23, 0.32, 0.77, 3.08, 1.73, 0.64, 0.87])

# LIP-1C bathymetry at 13h on xd grid
hd1C = -np.array(
    [
        -4.0935,
        -4.0931,
        -4.0926,
        -4.0921,
        -4.0917,
        -4.0912,
        -4.0733,
        -3.9439,
        -3.8105,
        -3.6717,
        -3.5289,
        -3.3835,
        -3.2361,
        -3.0864,
        -2.9344,
        -2.7828,
        -2.6401,
        -2.4939,
        -2.3813,
        -2.3146,
        -2.2765,
        -2.2534,
        -2.2351,
        -2.2143,
        -2.1864,
        -2.1494,
        -2.1028,
        -2.0481,
        -1.9878,
        -1.9250,
        -1.8632,
        -1.8055,
        -1.7545,
        -1.7112,
        -1.6749,
        -1.6427,
        -1.6085,
        -1.5627,
        -1.4999,
        -1.4240,
        -1.3120,
        -1.1622,
        -0.9755,
        -0.7739,
        -0.7177,
        -1.0002,
        -1.0464,
        -1.0356,
        -0.9975,
        -0.8333,
        -0.5603,
        -0.4902,
        -0.4558,
        -0.4270,
        -0.4016,
        -0.3659,
        -0.3035,
        -0.2013,
        -0.0536,
        0.1341,
        0.3454,
        0.5532,
        0.7225,
        0.8159,
    ]
)

x_hrms_d1C = np.array([20, 40, 65, 100, 115, 130, 132, 138, 145, 152, 160, 170])
hrms_d1C = np.array(
    [0.4, 0.41, 0.43, 0.44, 0.43, 0.43, 0.46, 0.43, 0.35, 0.33, 0.32, 0.22]
)

x_ubot_d1C = np.array([65, 102, 115, 125, 130, 134, 152, 160])
ubot_d1C = -np.array([1, 1, 1, 2, 2, 3, 13, 11])

x_Cbot_d1C = np.array([65, 102, 115, 125, 130, 134, 152, 160])
Cbot_d1C = np.array([0.1, 0.1, 0.3, 0.2, 0.35, 0.5, 0.3, 0.8])

# Select data based on case
if mycase == "1B":
    hd = hd1B
    x_hrms_d, hrms_d = x_hrms_d1B, hrms_d1B
    x_ubot_d, ubot_d = x_ubot_d1B, ubot_d1B
    x_Cbot_d, Cbot_d = x_Cbot_d1B, Cbot_d1B
    case_title = "SANDBAR EROSION LIP-1B"
else:
    hd = hd1C
    x_hrms_d, hrms_d = x_hrms_d1C, hrms_d1C
    x_ubot_d, ubot_d = x_ubot_d1C, ubot_d1C
    x_Cbot_d, Cbot_d = x_Cbot_d1C, Cbot_d1C
    case_title = "SANDBAR ACCRETION LIP-1C"

# =====================================================================
# Plot
# =====================================================================
fig, axs = plt.subplots(4, 1, figsize=(10, 12))
xmin, xmax = 60, 190

# --- Subplot 1: U velocity section ---
ax1 = axs[0]
cmin, cmax, nbcol = -0.5, 0.5, 20
levels = np.linspace(cmin, cmax, nbcol + 1)
cf1 = ax1.contourf(xu2d, zru, u_plot, levels=levels, cmap="jet", extend="both")
fig.colorbar(cf1, ax=ax1, orientation="horizontal", pad=0.02, aspect=40)
ax1.plot(xr, -hr, "k:", linewidth=3, label="Initial")
ax1.plot(xr, -hnew, "k", linewidth=3, label="Final Model")
ax1.plot(xd, -hd, "r", linewidth=3, label="Final Data")
ax1.plot(xr, zeta_plot, "g", linewidth=3)
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(-2.5, 0.5)
ax1.set_ylabel("Depth [m]", fontsize=12)
ax1.legend(loc="lower right", fontsize=10)
ax1.grid(True)
ax1.set_title(f"{case_title} - U [m/s] at Time {thour} hour", fontsize=12)
ax1.tick_params(labelsize=12)

# --- Subplot 2: Hrms ---
ax2 = axs[1]
ax2.plot(xr, hrms_plot, "b", linewidth=2, label="Model")
ax2.plot(x_hrms_d, hrms_d, "b*", markersize=10, label="Flume")
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(0, 1)
ax2.set_ylabel("Hrms [m]", fontsize=12)
ax2.legend(loc="lower left", fontsize=10)
ax2.grid(True)
ax2.tick_params(labelsize=12)

# --- Subplot 3: Undertow ---
ax3 = axs[2]
ax3.plot(xu, 100 * ubot, "b", linewidth=2, label="Model")
ax3.plot(x_ubot_d, ubot_d, "b*", markersize=10, label="Flume")
ax3.set_xlim(xmin, xmax)
ax3.set_ylim(-40, 0)
ax3.set_ylabel("U [cm/s]", fontsize=12)
ax3.legend(loc="lower left", fontsize=10)
ax3.grid(True)
ax3.tick_params(labelsize=12)

# --- Subplot 4: Sand concentration ---
ax4 = axs[3]
ax4.plot(xr, Cbot, "b", linewidth=2, label="Model")
ax4.plot(x_Cbot_d, Cbot_d, "b*", markersize=10, label="Flume")
ax4.set_xlim(xmin, xmax)
ax4.set_ylim(0, 4)
ax4.set_xlabel("X [m]", fontsize=12)
ax4.set_ylabel("C [g/l]", fontsize=12)
ax4.legend(loc="lower left", fontsize=10)
ax4.grid(True)
ax4.tick_params(labelsize=12)

plt.tight_layout()

# Save outputs
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, f"sandbar_LIP_{mycase}.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF saved: {pdf_path}")

if args.makepng:
    png_path = os.path.join(args.output_dir, f"sandbar_LIP_{mycase}.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG saved: {png_path}")

# Show or suppress plots
if not args.no_show:
    plt.show()
else:
    plt.close()
    print("Plot display suppressed (--no-show used).")
