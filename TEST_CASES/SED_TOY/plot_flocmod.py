#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot results from the FLOCMOD flocculation test cases.

Two modes are available:
  - 0D : zero-dimensional box test case (sed_toy_floc_his.nc).
         Produces a d50 time series + stacked particle-size distribution.
  - 1DV: one-dimensional vertical water-column test case.
         Produces time-depth colormaps of shear rate, SSC and d50,
         plus time-series at three depths and stacked PSD panels.

Default sediment classes (15 classes, MUSTANG/FLOCMOD naming):
  diameters : 4 µm → 1500 µm  (see --diam)
  variables : SED1 … SED15    (see --listvar)
"""

import os
import numpy as np
import argparse
import matplotlib
import matplotlib.pyplot as plt
import xarray


# ── Helper functions ──────────────────────────────────────────────────────────


def compute_dperc(mes, diam, perc):
    proportion = 0.0
    dperc = 0.0
    nbsed = len(diam)
    n = 0
    while (proportion < perc) and (n < nbsed):
        proportion = proportion + mes[n] / mes.sum()
        n = n + 1
    n = n - 1
    if n == nbsed - 1:
        dperc = diam[n]
    else:
        if n == 0:
            dperc = diam[n]
        else:
            a = (diam[n] - diam[n - 1]) / (mes[n] / mes.sum())
            dperc = diam[n - 1] + a * (perc - (proportion - mes[n] / mes.sum()))
    return dperc * 1000000.0


def comp_d50(data, listvar, diam):
    d10 = data[listvar[0]] * 0.0
    d50 = data[listvar[0]] * 0.0
    d90 = data[listvar[0]] * 0.0
    ssc = data[listvar[0]] * 0.0
    mes = np.zeros((data[listvar[0]].shape[0], data[listvar[0]].shape[1], len(diam)))
    for i, var in enumerate(listvar):
        mes[:, :, i] = data[var]
        ssc = ssc + data[var]
    for ik, k in enumerate(data.s_rho.values):
        for it, t in enumerate(data.time.values):
            d50[it, ik] = compute_dperc(mes[it, ik, :], diam, 0.5)
            d10[it, ik] = compute_dperc(mes[it, ik, :], diam, 0.1)
            d90[it, ik] = compute_dperc(mes[it, ik, :], diam, 0.9)
    return d10, d50, d90, ssc


def plot_d50(ax, data, listvar, diam, inth=1, k=0, xlab="time", ylab=""):
    if listvar[0] == "SED1":
        times = matplotlib.dates.date2num(data.time.values)
    else:
        times = data.time.values
    d10, d50, d90, ssc = comp_d50(data, listvar, diam)
    ax.plot(times, d50[:, k], color="black", label="FLOCMOD", linewidth=1)
    ax.set_xlim([times[0], times[-1]])
    if listvar[0] == "SED1":
        ax.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / inth) * inth)
        )
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    return ax, times, d10, d50, d90, ssc


def plot_class_distri(
    ax,
    data,
    listvar,
    diam,
    inth=1,
    k=0,
    xlab="time",
    ylab="particle size distribution",
    add_legend=True,
    add_label=True,
):
    if listvar[0] == "SED1":
        times = matplotlib.dates.date2num(data.time.values)
    else:
        times = data.time.values
    dt = times[1] - times[0]
    datak1 = data.isel(s_rho=k)
    mes = np.zeros((len(diam), len(times)))
    prop_mes = np.zeros((len(diam), len(times)))
    idiamt = np.zeros((len(diam), len(times)))
    for it, t in enumerate(times):
        datait = datak1.isel(time=it)
        for i, var in enumerate(listvar):
            mes[i, it] = datait[var]
            idiamt[i, it] = i
        for i, var in enumerate(listvar):
            prop_mes[i, it] = mes[i, it] / mes[:, it].sum()
    colnorm = matplotlib.colors.Normalize(vmin=0.0, vmax=len(diam))
    varcolor = matplotlib.cm.ScalarMappable(norm=colnorm, cmap="jet")
    varcolor.set_array(idiamt)
    col = varcolor.to_rgba(idiamt)
    offset_class = 0.0 * times
    for i in range(len(diam)):
        if i > 0:
            offset_class = offset_class + prop_mes[i - 1, :]
        if add_label:
            ax.bar(
                times,
                prop_mes[i, :],
                bottom=offset_class,
                color=col[i, :, :],
                width=dt,
                linewidth=0.5,
                edgecolor="grey",
                label="%i" % int(diam[i] * 1e6),
            )
        else:
            ax.bar(
                times,
                prop_mes[i, :],
                bottom=offset_class,
                color=col[i, :, :],
                width=dt,
                linewidth=0.5,
                edgecolor="grey",
            )
    ax.set_xlim([times[0], times[-1]])
    if listvar[0] == "SED1":
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        ax.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / inth) * inth)
        )
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if add_legend:
        ax.legend(
            loc="upper right",
            bbox_to_anchor=(1.12, 1.1),
            fontsize="small",
            labelspacing=0.2,
        )
    return ax


# ── Plot functions ────────────────────────────────────────────────────────────


def trace0D(diam, listvar, namecas, ncfic, output_dir="."):
    datafic = xarray.open_dataset(ncfic)
    data = datafic.isel(eta_rho=3, xi_rho=3).squeeze()

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=[10, 5])
    ax1, times, d10, d50, d90, ssc = plot_d50(
        ax1, data, listvar, diam, xlab="", ylab="Median floc diameter (\u03bcm)"
    )
    ax1.set_xlim([times[0], times[-1]])
    ax1.set_ylim([0, 300.0])
    if listvar[0] == "SED1":
        ax1.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / 3) * 3)
        )
        ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax1.grid()
    ax2 = plot_class_distri(
        ax2,
        data,
        listvar,
        diam,
        inth=3,
        xlab="",
        ylab="Stacked Floc size\ndistribution (fraction)",
    )
    ax2.set_ylim([0.0, 1.0])
    out_path = os.path.join(output_dir, "%s.png" % namecas)
    plt.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f"  saved: {out_path}")


def trace1DV(
    diam,
    listvar,
    namecas,
    ncfic,
    shear_rate_mod=True,
    eta0=3,
    xi0=3,
    srholist=[0, 25, 49],
    maxssc=0.5,
    maxd50=700,
    limz=[-5, 0],
    output_dir=".",
):
    datafic = xarray.open_dataset(ncfic)
    data = datafic.isel(eta_rho=eta0, xi_rho=xi0).squeeze()
    d10, d50, d90, ssc = comp_d50(data, listvar, diam)
    z = data.s_rho * (data.h + data.zeta) + data.zeta

    if listvar[0] == "SED1":
        ustar = np.sqrt(data.TAUSKIN / 1025.0)
        times = matplotlib.dates.date2num(data.time.values)
    else:
        ustar = 0.0 * data.wstr
        times = data.time.values
    dist_surf_on_bottom = (data.s_w[-1] - data.s_rho) / (data.s_rho - data.s_w[0])
    shear_rate = np.sqrt(
        ustar**3.0
        / (0.4 * (data.h + data.zeta) * (1.0e-6 + 1.0e-09))
        * dist_surf_on_bottom
    )
    if shear_rate_mod:
        N = len(data.s_rho)
        Eps_gls = data.AKt
        for k in range(N):
            if k == 0:
                diss = Eps_gls.isel(s_w=k + 1)
            elif k == N - 1:
                diss = Eps_gls.isel(s_w=k)
            else:
                diss = 0.5 * (Eps_gls.isel(s_w=k + 1) + Eps_gls.isel(s_w=k))
            shear_rate[:, k] = np.sqrt(diss / 1.0e-6)

    # Figure 1 : time-depth colormaps
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        4, 1, sharex=True, figsize=[9, 6], constrained_layout=True
    )
    ax1.plot(times, ustar)
    ax1.set_ylabel("u* (m/s)")
    ax1.set_ylim([0.001, 0.1])
    ax1.set_yscale("log")
    ax1.set_title("u* (m/s)")
    p2 = ax2.pcolormesh(
        times,
        z,
        np.transpose(shear_rate.values),
        norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=100.0),
        cmap="jet",
        shading="auto",
    )
    ax2.set_ylabel("z (m)")
    ax2.set_yticks(limz)
    ax2.set_title("Shear rate (/s)")
    fig.colorbar(p2, ax=ax2, extend="max", aspect=10, pad=0.01)
    p3 = ax3.pcolormesh(
        times,
        z,
        np.transpose(ssc.values),
        vmin=0.0,
        vmax=maxssc,
        cmap="jet",
        shading="nearest",
    )
    ax3.set_ylabel("z (m)")
    ax3.set_yticks(limz)
    ax3.set_title("SSC (g/L)")
    fig.colorbar(p3, ax=ax3, extend="max", aspect=10, pad=0.01)
    p4 = ax4.pcolormesh(
        times,
        z,
        np.transpose(d50.values),
        vmin=0.0,
        vmax=maxd50,
        cmap="jet",
        shading="nearest",
    )
    ax4.set_xlabel("time")
    ax4.set_ylabel("z (m)")
    ax4.set_yticks(limz)
    ax4.set_title("Floc d50 (\u03bcm)")
    fig.colorbar(p4, ax=ax4, extend="max", aspect=10, pad=0.01)
    ax4.set_xlim([times[-144], times[-1]])
    if listvar[0] == "SED1":
        ax4.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / 3) * 3)
        )
        ax4.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    for ax in (ax1, ax2, ax3, ax4):
        ax.grid()
    out1 = os.path.join(output_dir, "%s.png" % namecas)
    plt.savefig(out1, dpi=200)
    plt.close(fig)
    print(f"  saved: {out1}")

    # Figure 2 : time-series at 3 depths
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, sharex=True, figsize=[9, 5], constrained_layout=True
    )
    ax1.plot(times, shear_rate.isel(s_rho=srholist[0]), label="Bottom", color="red")
    ax1.plot(
        times, shear_rate.isel(s_rho=srholist[1]), label="Mid depth", color="green"
    )
    ax1.plot(times, shear_rate.isel(s_rho=srholist[2]), label="Surface", color="blue")
    ax1.set_ylabel("Shear rate (/s)")
    ax1.set_ylim([0.2, 200])
    ax1.set_yscale("log")
    ax2.plot(times, ssc.isel(s_rho=srholist[0]), label="Bottom", color="red")
    ax2.plot(times, ssc.isel(s_rho=srholist[1]), label="Mid depth", color="green")
    ax2.plot(times, ssc.isel(s_rho=srholist[2]), label="Surface", color="blue")
    ax2.set_ylabel("SSC (g/L)")
    ax2.set_ylim([0.0, maxssc])
    ax3.plot(times, d50.isel(s_rho=srholist[0]), label="Bottom", color="red")
    ax3.plot(times, d50.isel(s_rho=srholist[1]), label="Mid depth", color="green")
    ax3.plot(times, d50.isel(s_rho=srholist[2]), label="Surface", color="blue")
    ax3.set_ylabel("d50 (\u03bcm)")
    ax3.set_ylim([0.0, maxd50])
    ax3.set_xlim([times[-144], times[-1]])
    if listvar[0] == "SED1":
        ax3.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / 3) * 3)
        )
        ax3.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax1.grid()
    ax1.legend()
    ax2.grid()
    ax3.grid()
    out2 = os.path.join(output_dir, "%s_2.png" % namecas)
    plt.savefig(out2, dpi=200)
    plt.close(fig)
    print(f"  saved: {out2}")

    # Figure 3 : stacked PSD at 3 depths
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=[9, 5])
    plot_class_distri(
        ax1,
        data,
        listvar,
        diam,
        inth=1,
        k=srholist[2],
        xlab="",
        ylab="surface\nPSD",
        add_legend=False,
        add_label=False,
    )
    plot_class_distri(
        ax2,
        data,
        listvar,
        diam,
        inth=1,
        k=srholist[1],
        xlab="",
        ylab="mid depth\nPSD",
        add_legend=False,
        add_label=False,
    )
    plot_class_distri(
        ax3,
        data,
        listvar,
        diam,
        inth=1,
        k=srholist[0],
        xlab="",
        ylab="bottom\nPSD",
        add_legend=False,
        add_label=True,
    )
    ax3.set_xlim([times[-144], times[-1]])
    if listvar[0] == "SED1":
        ax3.xaxis.set_major_locator(
            matplotlib.dates.HourLocator(byhour=np.arange(24 / 3) * 3)
        )
        ax3.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    plt.subplots_adjust(right=0.8)
    fig.legend(loc="center right", bbox_to_anchor=(0.95, 0.5), labelspacing=0.2)
    out3 = os.path.join(output_dir, "%s_3.png" % namecas)
    plt.savefig(out3, dpi=200)
    plt.close(fig)
    print(f"  saved: {out3}")


# ── CLI ───────────────────────────────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Plot FLOCMOD results (0D or 1DV test case).",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--mode",
    type=str,
    choices=["0D", "1DV"],
    required=True,
    help="Test-case mode: '0D' (box) or '1DV' (vertical column).",
)
parser.add_argument(
    "--file",
    type=str,
    default=None,
    help=(
        "Path to the NetCDF history file.\n"
        "Default: '../cases/0D/sed_toy_floc_his.nc'  (0D)\n"
        "         '../cases/1DV/sed_toy_floc_his.nc' (1DV)"
    ),
)
parser.add_argument(
    "--name",
    type=str,
    default=None,
    help="Case name used as prefix for output PNG files.\n"
    "Default: 'flocmod_0D' (0D) or 'flocmod_1DV' (1DV).",
)
parser.add_argument(
    "--diam",
    type=float,
    nargs="+",
    default=None,
    help="Sediment class diameters in metres (space-separated).\n"
    "Default: 15-class standard FLOCMOD distribution.",
)
parser.add_argument(
    "--listvar",
    type=str,
    nargs="+",
    default=None,
    help="Sediment variable names in the NetCDF file. Default: SED1 … SED15.",
)
parser.add_argument(
    "--gls",
    action="store_true",
    help="(1DV only) Use GLS turbulent dissipation (AKt) for shear rate.",
)
parser.add_argument(
    "--eta0",
    type=int,
    default=3,
    help="(1DV only) eta_rho extraction index (default: 3).",
)
parser.add_argument(
    "--xi0",
    type=int,
    default=3,
    help="(1DV only) xi_rho extraction index (default: 3).",
)
parser.add_argument(
    "--srholist",
    type=int,
    nargs=3,
    default=[0, 25, 49],
    metavar=("BOTTOM", "MID", "SURFACE"),
    help="(1DV only) s_rho layer indices: bottom / mid / surface. Default: 0 25 49.",
)
parser.add_argument(
    "--maxssc",
    type=float,
    default=0.5,
    help="(1DV only) SSC colormap max (g/L). Default: 0.5.",
)
parser.add_argument(
    "--maxd50",
    type=float,
    default=700.0,
    help="(1DV only) d50 colormap max (µm). Default: 700.",
)
parser.add_argument(
    "--limz",
    type=float,
    nargs=2,
    default=[-5.0, 0.0],
    metavar=("ZMIN", "ZMAX"),
    help="(1DV only) Vertical axis limits in metres. Default: -5 0.",
)
parser.add_argument(
    "--output-dir",
    type=str,
    default=".",
    help="Directory where all PNG files are saved (default: '.').",
)
parser.add_argument(
    "--makepdf", action="store_true", help="(kept for template compatibility, not used)"
)
parser.add_argument(
    "--makepng", action="store_true", help="(kept for template compatibility, not used)"
)
parser.add_argument(
    "--no-show", action="store_true", help="(kept for template compatibility, not used)"
)

args = parser.parse_args()

# ── Defaults ──────────────────────────────────────────────────────────────────

DIAM_DEFAULT = {
    "0D": np.array(
        [
            0.000004,
            0.000006,
            0.000009,
            0.000014,
            0.000021,
            0.000033,
            0.000050,
            0.000077,
            0.000118,
            0.000180,
            0.000275,
            0.000421,
            0.000643,
            0.000982,
            0.001500,
        ]
    ),
    "1DV": np.array(
        [
            0.000050,
            0.000077,
            0.000118,
            0.000180,
            0.000275,
            0.000421,
            0.000643,
            0.000982,
        ]
    ),
}

LISTVAR_DEFAULT = {
    "0D": [
        "SED1",
        "SED2",
        "SED3",
        "SED4",
        "SED5",
        "SED6",
        "SED7",
        "SED8",
        "SED9",
        "SED10",
        "SED11",
        "SED12",
        "SED13",
        "SED14",
        "SED15",
    ],
    "1DV": [
        "SED1",
        "SED2",
        "SED3",
        "SED4",
        "SED5",
        "SED6",
        "SED7",
        "SED8",
    ],
}
FILE_DEFAULT = {
    "0D": "../cases/0D/sed_toy_floc_his.nc",
    "1DV": "../cases/1DV/sed_toy_floc_his.nc",
}
NAME_DEFAULT = {
    "0D": "flocmod_0D",
    "1DV": "flocmod_1DV",
}

diam = np.array(args.diam) if args.diam else DIAM_DEFAULT[args.mode]
listvar = args.listvar if args.listvar else LISTVAR_DEFAULT[args.mode]
ncfic = args.file if args.file else FILE_DEFAULT[args.mode]
namecas = args.name if args.name else NAME_DEFAULT[args.mode]
output_dir = args.output_dir

# ── Validate ──────────────────────────────────────────────────────────────────

os.makedirs(output_dir, exist_ok=True)

if not os.path.isfile(ncfic):
    print(f"Error: NetCDF file '{ncfic}' not found.")
    raise SystemExit(1)

if len(diam) != len(listvar):
    print(
        f"Error: --diam ({len(diam)}) and --listvar ({len(listvar)}) "
        "must have the same length."
    )
    raise SystemExit(1)

# ── Run ───────────────────────────────────────────────────────────────────────

if args.mode == "0D":
    print(
        f"[0D] case='{namecas}'  file='{ncfic}'  "
        f"nclasses={len(diam)}  output='{output_dir}'"
    )
    trace0D(diam, listvar, namecas, ncfic, output_dir=output_dir)
    print("[0D] Done.")

elif args.mode == "1DV":
    shear_rate_mod = args.gls
    print(
        f"[1DV] case='{namecas}'  file='{ncfic}'  "
        f"nclasses={len(diam)}  output='{output_dir}'\n"
        f"      shear_rate_mod={shear_rate_mod}  eta0={args.eta0}  xi0={args.xi0}\n"
        f"      srholist={args.srholist}  maxssc={args.maxssc}  "
        f"maxd50={args.maxd50}  limz={args.limz}"
    )
    trace1DV(
        diam,
        listvar,
        namecas,
        ncfic,
        shear_rate_mod=shear_rate_mod,
        eta0=args.eta0,
        xi0=args.xi0,
        srholist=args.srholist,
        maxssc=args.maxssc,
        maxd50=args.maxd50,
        limz=args.limz,
        output_dir=output_dir,
    )
    print("[1DV] Done.")
