#!/usr/bin/env python3
"""
INNERSHELF Test Case
Compare Ekman 2D theoretical model solution with CROCO numerical solution:
vg, psi/|Uek| & v

References:
- Estrade P., P. Marchesiello, A. Colin de Verdiere, C. Roy, 2008:
  Cross-shelf structure of coastal upwelling: a two-dimensional
  expansion of Ekman's theory and a mechanism for innershelf upwelling
  shut down. Journal of Marine Research, 66, 589-616.

- Marchesiello P., and P. Estrade, 2010: Upwelling limitation by
  geostrophic onshore flow. Journal of Marine Research, 68, 37-62.

Translated from Matlab by Claude (2025)
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import croco_utils as cr
import argparse
import os

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Plot results from the INNERSHELF test case.",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file", type=str, default="inner_avg.nc", help="Path to the NetCDF file"
)
parser.add_argument("--makepdf", action="store_true", help="Save plot as a PDF file")
parser.add_argument("--makepng", action="store_true", help="Save plot as PNG file")
parser.add_argument("--no-show", action="store_true", help="Suppress plot display")
parser.add_argument(
    "--output-dir", type=str, default=".", help="Directory to save output files"
)
args = parser.parse_args()

# Ensure output directory exists
os.makedirs(args.output_dir, exist_ok=True)

print("=== INNERSHELF Test Case ===")
print("Comparing Ekman 2D theoretical vs CROCO numerical solutions")
if args.makepdf:
    print("PDF output enabled")
if args.makepng:
    print("PNG output enabled")
if args.no_show:
    print("Display disabled")

# User defined parameters
# Model params
fname = args.file  # croco file name from argument
tauy = 0.07  # alongshore wind-stress (Pa)
taux = 0  # onshore wind-stress (Pa)
ug0 = 2  # onshore geostrophic flow (cm/s)
g = 9.81  # gravity acceleration (m²/s)
rho0 = 1000  # reference density (kg/m³)

# Figure params
hoDmax = 2.5  # index of offshore limit h/D
tick = np.arange(-2.4, -0.1, 0.2)
lx = 18
ly = 24
dxleft = 2
dxright = 0.25
dyup = 2.5
dybot = 1.5
dyint = 0.5
subfigratio = 2 / 3

print("Reading grid and parameters from CROCO file...")

# Get grid, Av & f from numerical model
with Dataset(fname, "r") as nc:
    tindex = len(nc.variables["scrum_time"][:]) - 1  # reads last record (0-based)
    yindex = 2  # y index (Python 0-based, so 3->2)

    Av = nc.Akv_bak  # viscosity (global attribute)
    f = nc.variables["f"][yindex, 0]  # Coriolis frequency
    D = np.pi * np.sqrt(2 * Av / abs(f))  # Ekman depth

    print(f"Coriolis parameter f = {f:.2e} s⁻¹")
    print(f"Ekman depth D = {D:.1f} m")

    # Horizontal grid
    hr = np.squeeze(nc.variables["h"][yindex, :])
    tmp = np.where(hr / D > hoDmax)[0]  # reduce offshore grid
    xindex = tmp[-1] if len(tmp) > 0 else 0
    hr = hr[xindex:]
    L = len(hr)
    hu = 0.5 * (hr[:-1] + hr[1:])
    xr = np.squeeze(nc.variables["x_rho"][yindex, xindex:])
    dx = xr[1] - xr[0]
    xu = 0.5 * (xr[:-1] + xr[1:])

    # Vertical grid
    N = len(nc.dimensions["s_rho"])
    theta_s = nc.theta_s  # S-coordinate surface control parameter (attribute)
    theta_b = nc.theta_b  # S-coordinate bottom control parameter (attribute)
    hc = nc.hc  # S-coordinate critical depth (attribute)
    zeta = np.squeeze(nc.variables["zeta"][tindex, yindex, xindex:])

    # Compute vertical coordinates
    zr = cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "r", 1)
    dzr = zr[1:, :] - zr[:-1, :]  # --> zw(2:N,:)
    zru = 0.5 * (zr[:, :-1] + zr[:, 1:])
    dzru = zru[1:, :] - zru[:-1, :]  # --> zwu(2:N,:)
    zw = cr.zlevs(hr, zeta, theta_s, theta_b, hc, N, "w", 1)
    dzw = zw[1:, :] - zw[:-1, :]  # --> zr
    zwu = 0.5 * (zw[:, :-1] + zw[:, 1:])
    dzwu = zwu[1:, :] - zwu[:-1, :]  # --> zru

    # Topographic slope
    hx = np.zeros_like(hr)
    hx[1:-1] = (hu[1:] - hu[:-1]) / (xu[1:] - xu[:-1])
    hx = np.tile(hx, (N + 1, 1))

    print("Reading numerical model fields...")

    # Numerical model fields
    # Geostrophic velocity vg1 = g/f dzeta --> xu
    vg1 = g / f / dx * (zeta[1:] - zeta[:-1])

    # Zonal velocity --> xu,zru
    u1 = np.squeeze(nc.variables["u"][tindex, :, yindex, xindex:])

    # Meridional velocity --> xr,zr
    v1 = np.squeeze(nc.variables["v"][tindex, :, yindex, xindex:])

    # Vertical velocity --> xr,zw
    w1 = np.squeeze(nc.variables["w"][tindex, :, yindex, xindex:])

    # Stream function --> xu,zwu
    psi1 = np.zeros((N + 1, L - 1))
    for k in range(N):
        psi1[k + 1, :] = psi1[k, :] - u1[k, :] * dzwu[k, :]

print("Computing theoretical fields...")

# Compute theoretical fields on rho and u grids
v0 = np.pi * tauy / (rho0 * abs(f) * D)
ug = ug0 * 0.01
gamma = f / abs(f)
c = (1 + 1j * gamma) * np.pi / D
cu0 = (1 - gamma * 1j) * np.pi * complex(taux, tauy) / (rho0 * abs(f) * D)

# Compute geostrophic velocities --> xr
qsir = np.pi * hr / D
alphar = (np.cosh(qsir) ** 2) * (np.cos(qsir) ** 2) + (np.sinh(qsir) ** 2) * (
    np.sin(qsir) ** 2
)
s1r = np.cosh(qsir) * np.cos(qsir) / alphar
s2r = np.sinh(qsir) * np.sin(gamma * qsir) / alphar
t1r = np.cosh(qsir) * np.sinh(qsir) / alphar
t2r = np.cos(qsir) * np.sin(gamma * qsir) / alphar
f1r = np.sinh(qsir) * np.cos(qsir)
f2r = np.cosh(qsir) * np.sin(gamma * qsir)

vgr = 2 * np.pi / (rho0 * f * D) * ((1 - s1r) * tauy + s2r * taux) - ug * (
    t1r + gamma * t2r - 2 * qsir
)
vgr = vgr / (gamma * t1r - t2r)
ugr = np.zeros_like(hr) + ug
cugr = ugr + 1j * vgr

func = (
    2
    + (np.cosh(qsir) * np.sin(qsir) - np.sinh(qsir) * np.cos(qsir))
    / (t1r - t2r)
    / alphar
    - 2
    * (1 - s1r)
    * (np.cosh(qsir) ** 2 - np.cos(qsir) ** 2)
    / alphar
    / ((t1r - t2r) ** 2)
)

# --> xu
qsiu = np.pi * hu / D
alphau = (np.cosh(qsiu) ** 2) * (np.cos(qsiu) ** 2) + (np.sinh(qsiu) ** 2) * (
    np.sin(qsiu) ** 2
)
s1u = np.cosh(qsiu) * np.cos(qsiu) / alphau
s2u = np.sinh(qsiu) * np.sin(gamma * qsiu) / alphau
t1u = np.cosh(qsiu) * np.sinh(qsiu) / alphau
t2u = np.cos(qsiu) * np.sin(gamma * qsiu) / alphau
f1u = np.sinh(qsiu) * np.cos(qsiu)
f2u = np.cosh(qsiu) * np.sin(gamma * qsiu)

vgu = 2 * np.pi / (rho0 * f * D) * ((1 - s1u) * tauy + s2u * taux) - ug * (
    t1u + gamma * t2u - 2 * qsiu
)
vgu = vgu / (gamma * t1u - t2u)
ugu = np.zeros_like(hu) + ug
cugu = ugu + 1j * vgu

# Compute total horizontal velocity --> zr
wer = np.zeros((N, L), dtype=complex)
for j in range(L):
    wer[:, j] = cu0 * np.sinh(c * (zr[:, j] + hr[j])) / np.cosh(c * hr[j]) - cugr[
        j
    ] * np.cosh(c * zr[:, j]) / np.cosh(c * hr[j])

# --> zu
weu = np.zeros((N, L - 1), dtype=complex)
for j in range(L - 1):
    weu[:, j] = cu0 * np.sinh(c * (zru[:, j] + hu[j])) / np.cosh(c * hu[j]) - cugu[
        j
    ] * np.cosh(c * zru[:, j]) / np.cosh(c * hu[j])

# Compute stream function --> zwu
psi2 = np.zeros((N + 1, L - 1), dtype=complex)
for j in range(L - 1):
    psi2[:, j] = cu0 * (
        1 - np.cosh(c * (zwu[:, j] + hu[j])) / np.cosh(c * hu[j])
    ) + cugu[j] * np.sinh(c * zwu[:, j]) / np.cosh(c * hu[j])
    psi2[:, j] = psi2[:, j] / c - cugu[j] * zwu[:, j]

psi2 = np.real(psi2)

# Compute vertical velocity --> zw
w2 = np.zeros((N + 1, L), dtype=complex)
for j in range(L):
    w2[:, j] = (
        -np.sinh(c * zw[:, j])
        * (s1r[j] - 1j * s2r[j])
        * (
            (s1r[j] - 1j * s2r[j])
            * (1 + 1j + 2j * vgr[j] / 2 / v0 * (f1r[j] + 1j * f2r[j]))
            - (1 + 1j) * func[j]
        )
    )

w2 = v0 * hx * np.real(w2)

print("Creating plots...")

# Create figure
fig = plt.figure(figsize=(lx / 2.54, ly / 2.54))  # Convert cm to inches
dy2 = (ly - 2 * dyint - dyup - dybot) / (2 + subfigratio)
dy1 = dy2 * subfigratio
dx = lx - dxleft - dxright

# Geostrophic velocity
ax1 = plt.subplot(3, 1, 1)
plt.plot(-hu / D, vgu, color=[0.7, 0.7, 0.7], linewidth=3, label="analytical")
plt.plot(-hu / D, vg1, "k--", linewidth=2, label="numerical")
plt.axis([-hoDmax, -hr[-1] / D, -0.2, 0.2])
plt.gca().tick_params(labelsize=16)
plt.yticks(np.arange(-0.2, 0.3, 0.1))
plt.gca().set_yticklabels(["", "-0.1", "0", "0.1", ""])
plt.xticks(tick)
plt.gca().set_xticklabels([])
plt.grid(True)
plt.text(-2.3, 0, r"$v_G$ [m/s]", fontsize=20)
plt.legend(loc="lower right")

# Normalized stream function
ax2 = plt.subplot(3, 1, 2)
tmpx = np.tile(-hu / D, (N + 1, 1))
Uek = abs(tauy / rho0 / f)
psirange1 = np.arange(10, 160, 10)
psirange2 = np.arange(-150, -5, 10)

# Analytical contours
cs1 = plt.contour(
    tmpx,
    zwu / D,
    -100 * psi2 / Uek,
    levels=psirange1,
    colors=[[0.7, 0.7, 0.7]],
    linewidths=3,
)
cs2 = plt.contour(
    tmpx,
    zwu / D,
    -100 * psi2 / Uek,
    levels=psirange2,
    colors=[[0.7, 0.7, 0.7]],
    linewidths=3,
    linestyles=":",
)

# Numerical contours
cs3 = plt.contour(
    tmpx, zwu / D, -100 * psi1 / Uek, levels=psirange1, colors="k", linewidths=2
)
plt.clabel(cs3)
cs4 = plt.contour(
    tmpx,
    zwu / D,
    -100 * psi1 / Uek,
    levels=psirange2,
    colors="k",
    linewidths=2,
    linestyles=":",
)
plt.clabel(cs4)

# Topography
plt.plot(-hr / D, -hr / D, "k")

# Legend (artifice)
plt.plot(
    [-hoDmax, -hu[-1] / D],
    [100, 100],
    color=[0.7, 0.7, 0.7],
    linewidth=3,
    label="analytical",
)
plt.plot([-hoDmax, -hu[-1] / D], [101, 101], "k:", linewidth=2, label="numerical")

plt.text(-1, -1.5, r"$\psi / U_{ek}$ [%]", fontsize=20)
plt.axis([-hoDmax, -hu[-1] / D, -hoDmax, zr[N - 1, 0] / D])
plt.xticks(tick)
plt.gca().set_xticklabels([])
plt.yticks(tick)
plt.gca().set_yticklabels(
    ["", "-2.2", "", "-1.8", "", "-1.4", "", "-1", "", "-0.6", "", "-0.2"]
)
plt.ylabel("z/D")
plt.grid(True)
plt.legend(loc="lower right")

# Normalized meridional Ekman velocity
ax3 = plt.subplot(3, 1, 3)
tmp = np.tile(vgr, (N, 1))
v2 = np.imag(wer) + tmp
vrange1 = np.arange(-100, 0, 5)
vrange2 = np.arange(5, 105, 5)
tmpx = np.tile(-hr / D, (N, 1))

# Analytical contours
cs5 = plt.contour(
    tmpx, zr / D, 100 * v2, levels=vrange1, colors="k", linewidths=2, linestyles=":"
)
cs6 = plt.contour(tmpx, zr / D, 100 * v2, levels=vrange2, colors="k", linewidths=2)
cs7 = plt.contour(
    tmpx, zr / D, 100 * v2, levels=[0], colors=[[0.7, 0.7, 0.7]], linewidths=2
)

# Numerical contours
cs8 = plt.contour(
    tmpx,
    zr / D,
    100 * v1,
    levels=vrange1,
    colors=[[0.7, 0.7, 0.7]],
    linewidths=2,
    linestyles=":",
)
plt.clabel(cs8)
cs9 = plt.contour(
    tmpx, zr / D, 100 * v1, levels=vrange2, colors=[[0.7, 0.7, 0.7]], linewidths=2
)
plt.clabel(cs9)
cs10 = plt.contour(
    tmpx, zr / D, 100 * v1, levels=[0], colors=[[0.7, 0.7, 0.7]], linewidths=2
)

# Topography
plt.plot(-hr / D, -hr / D, "k")

# Legend (artifice)
plt.plot(
    [-hoDmax, -hu[-1] / D],
    [100, 100],
    color=[0.7, 0.7, 0.7],
    linewidth=3,
    label="analytical",
)
plt.plot([-hoDmax, -hu[-1] / D], [101, 101], "k:", linewidth=2, label="numerical")

plt.text(-0.8, -1.5, r"$v$ [cm/s]", fontsize=20)
plt.axis([-hoDmax, -hu[-1] / D, -hoDmax, zr[N - 1, 0] / D])
plt.gca().tick_params(labelsize=16)
plt.xticks(tick)
plt.gca().set_xticklabels(
    ["", "-2.2", "", "-1.8", "", "-1.4", "", "-1", "", "-0.6", "", "-0.2"]
)
plt.yticks(tick)
plt.gca().set_yticklabels(
    ["", "-2.2", "", "-1.8", "", "-1.4", "", "-1", "", "-0.6", "", "-0.2"]
)
plt.xlabel("h/D")
plt.ylabel("z/D")
plt.grid(True)
plt.legend(loc="lower right")

# Set subplot positions (convert cm to figure coordinates)
pos1 = [dxleft / lx, (dybot + 2 * (dy2 + dyint)) / ly, dx / lx, dy1 / ly]
pos2 = [dxleft / lx, (dybot + dy2 + dyint) / ly, dx / lx, dy2 / ly]
pos3 = [dxleft / lx, dybot / ly, dx / lx, dy2 / ly]

ax1.set_position(pos1)
ax2.set_position(pos2)
ax3.set_position(pos3)

plt.tight_layout()

# Save and/or show figure
if args.makepdf:
    pdf_path = os.path.join(args.output_dir, "valid_innershelf.pdf")
    print(f"Saving figure as {pdf_path}...")
    plt.savefig(
        pdf_path, format="pdf", bbox_inches="tight", facecolor="white", transparent=True
    )
    print(f"PDF file '{pdf_path}' has been created.")

if args.makepng:
    png_path = os.path.join(args.output_dir, "valid_innershelf.png")
    print(f"Saving figure as {png_path}...")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight", facecolor="white")
    print(f"PNG file '{png_path}' has been created.")

if not args.no_show:
    plt.show()
else:
    plt.close()

print("Script completed successfully!")
