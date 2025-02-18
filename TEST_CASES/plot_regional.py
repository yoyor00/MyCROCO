#!/usr/bin/env python3

import os
import math
import argparse
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
import croco_utils as cr

####### MAIN
# Parse command-line arguments
parser = argparse.ArgumentParser(
    description="Plot data from a NetCDF file and optionally save as PDF or PNG.",
    epilog="Example usage:\n"
    + "  python script.py --makepdf --file croco_his.nc\n"
    + "  python script.py --makepng --no-show",
    formatter_class=argparse.RawTextHelpFormatter,
)
parser.add_argument(
    "--file",
    type=str,
    default="croco_his.nc",
    help="Path to the NetCDF file (default: basin.nc)",
)
parser.add_argument(
    "--makepdf", action="store_true", help="Generate a PDF of the plots"
)
parser.add_argument(
    "--makepng", action="store_true", help="Generate a PNG of the plots"
)
parser.add_argument(
    "--no-show",
    action="store_true",
    help="Do not display the plots on the screen (default: display)",
)
parser.add_argument(
    "--output-dir",
    type=str,
    default=".",
    help="Directory to save the output files (default: current directory)",
)
args = parser.parse_args()

# Read data from NetCDF file
try:
    hisfile = args.file
except FileNotFoundError:
    print(f"Error: File '{args.file}' not found.")
    exit(1)

# # Dev
# data_path = "/local/tmp/3/gcambon/DEV_CROCO/2024_CI_BENCH/croco/BENCH/rundir/BENGUELA_LR_TIDES_BULK/mpi-4/CROCO_FILES/"
# hisfile = data_path + "croco_his.nc"

#################################################@
# read lat,lon,mask
lat, lon, mask = cr.read_latlonmask(hisfile, "r")
latu, lonu, masku = cr.read_latlonmask(hisfile, "u")
latv, lonv, maskv = cr.read_latlonmask(hisfile, "v")

# get the number of vertical level
nc = Dataset(hisfile)
nc = Dataset(hisfile)
temp = nc.variables["temp"][:]
T = np.shape(temp)[0]
N = np.shape(temp)[1]

# choose type of plots
plot_horiz = 1
plot_section = 1

# Parameters for vertical sections
gridfile = hisfile  # Grid file
tindex = -1  # Last time step
ind_sec = 17  # Column index for the vertical slice

# FOR THE HORIZONTAL PLOT
if plot_horiz == 1:
    #
    # read the horizslice : last time step , surface level
    tndx = -1
    vlev = N
    # Read var
    ssh = cr.get_hslice(hisfile, hisfile, "zeta", tndx, 0, "r")
    sst = cr.get_hslice(hisfile, hisfile, "temp", tndx, vlev, "r")
    sss = cr.get_hslice(hisfile, hisfile, "salt", tndx, vlev, "r")
    ubar = cr.get_hslice(hisfile, hisfile, "ubar", tndx, 0, "u")
    vbar = cr.get_hslice(hisfile, hisfile, "vbar", tndx, 0, "v")
    usurf = cr.get_hslice(hisfile, hisfile, "u", tndx, vlev, "u")
    vsurf = cr.get_hslice(hisfile, hisfile, "v", tndx, vlev, "v")
    sustr = cr.get_hslice(hisfile, hisfile, "sustr", tndx, 0, "u")
    svstr = cr.get_hslice(hisfile, hisfile, "svstr", tndx, 0, "v")

    # Do the subplot
    # Liste des variables à tracer
    variables = {
        "SSH": (lon, lat, ssh, "Sea Surface Height", "m"),
        "SST": (lon, lat, sst, "Sea Surface Temperature", "°C"),
        "SSS": (lon, lat, sss, "Sea Surface Salinity", "PSU"),
        "U-Bar": (lonu, latu, ubar, "Barotropic U", "m/s"),
        "V-Bar": (lonv, latv, vbar, "Barotropic V", "m/s"),
        "U-Surf": (lonu, latu, usurf, "Surface U", "m/s"),
        "V-Surf": (lonv, latv, vsurf, "Surface V", "m/s"),
        "U-Wstress": (lonu, latu, usurf, "U Wind Stress", "N/m^{2}"),
        "V-Wstress": (lonv, latv, vsurf, "V Wind Stress", "N/m^{2}"),
    }
    # Création des subplots
    n_variables = len(variables)
    ncols = 3  # Nombre de colonnes
    nrows = (n_variables + ncols - 1) // ncols  # Calcul du nombre de lignes nécessaires
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(15, 5 * nrows),
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    axes = axes.flatten()

    # Parcours des variables
    for idx, (key, (lon, lat, data, title, units)) in enumerate(variables.items()):
        ax = axes[idx]
        ax.set_title(rf"{title} ($\mathrm{{{units}}}$)")

        # Ajout des données
        mesh = ax.pcolormesh(
            lon,
            lat,
            data,
            shading="auto",
            cmap="jet",
            norm=Normalize(vmin=np.nanmin(data), vmax=np.nanmax(data)),
        )

        # # Ajout des lignes de côte
        # ax.coastlines()

        # Ajout des sections sur la cartes
        ax.plot(
            lon[ind_sec, :],
            lat[ind_sec, :],
            color="black",
            linestyle="--",
            linewidth=1.5,
        )
        ax.plot(
            lon[:, ind_sec],
            lat[:, ind_sec],
            color="black",
            linestyle="--",
            linewidth=1.5,
        )

        # Ajout des longitudes et latitudes
        ax.set_xticks(
            np.arange(np.floor(lon.min()), np.ceil(lon.max()) + 1, 5),
            crs=ccrs.PlateCarree(),
        )
        ax.set_yticks(
            np.arange(np.floor(lat.min()), np.ceil(lat.max()) + 1, 5),
            crs=ccrs.PlateCarree(),
        )
        ax.xaxis.set_major_formatter(cticker.LongitudeFormatter())
        ax.yaxis.set_major_formatter(cticker.LatitudeFormatter())
        ax.tick_params(labelsize=8)

        # Ajout de la colorbar
        plt.colorbar(mesh, ax=ax, orientation="vertical", fraction=0.046, pad=0.04)

    # Supprimer les axes inutilisés
    for idx in range(n_variables, len(axes)):
        fig.delaxes(axes[idx])

    # Ajuster l'espacement des subplots
    plt.tight_layout()
    # plt.show()
    # plt.savefig('REGIONAL_maps.png')

#################################
# Determine output file paths
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(output_dir, "REGIONAL_maps.pdf")
    plt.savefig(pdf_path, transparent=True)
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(output_dir, "REGIONAL_maps.png")
    plt.savefig(png_path, dpi=300)
    print(f"PNG file '{png_path}' has been created.")

# Show plots if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --show to enable).")

############################################################################
# FOR THE VERTICAL SECTION PLOT
if plot_section == 1:
    # For temp, salt, u, v do merdian and zonal section crossing the mask
    gridfile = hisfile

    tindex = -1  # Dernier pas de temps
    direction = "x"  # Coupe zonale
    # ind_sec = 15  # Indice de la colonne pour la coupe

    # Extraire les données
    var_sec = "temp"
    section = cr.get_vertical_section(
        hisfile, var_sec, tindex, direction, ind_sec, gridfile
    )

    # Tracer la section verticale
    cr.plot_vertical_section(
        section, title=f"Zonal Vertical section -  {var_sec} - j={ind_sec}", unit="°C"
    )

# Exemple de dictionnaire avec variables et leurs directions
variables = {
    "temp": ["x", "y"],  # Coupe zonale et méridienne pour la température
    "salt": ["x", "y"],  # Coupe zonale et méridienne pour la salinité
    "u": ["x", "y"],  # Zonal et méridien pour la vitesse zonale
    "v": ["x", "y"],  # Zonal et méridien pour la vitesse méridienne
}

titles = {
    "temp": "Temperature (°C)",
    "salt": "Salinity (PSU)",
    "u": "Zonal Velocity (m/s)",
    "v": "Meridional Velocity (m/s)",
}

units = {
    "temp": "°C",
    "salt": "PSU",
    "u": "m/s",
    "v": "m/s",
}

# Calculer la disposition des sous-graphiques (grille optimale)
n_vars = len(variables)
n_dirs = sum(
    len(directions) for directions in variables.values()
)  # Nombre total de sections
n_cols = 2  # Nombre de colonnes fixes
n_rows = math.ceil(n_dirs / n_cols)  # Nombre de lignes nécessaires

# Initialisation des subplots
fig, axes = plt.subplots(
    n_rows, n_cols, figsize=(15, n_rows * 3)
)  # , constrained_layout=True

# Aplatir les axes pour gestion facile (en cas de grille)
axes = axes.flatten()

# Initialisation de l'indice pour les axes
axis_index = 0

# Boucle pour parcourir chaque variable et ses directions
for var, directions in variables.items():
    for direction in directions:
        # Extraire les données de la section
        print(f"Processing {var} in direction {direction}")
        section = cr.get_vertical_section(
            hisfile, var, tindex, direction, ind_sec, gridfile
        )
        data = section["variable"]
        m_data = np.ma.masked_equal(data, 0)
        # Configurer les axes
        ax = axes[axis_index]
        im = ax.pcolor(
            section["distance"],
            section["depth"],
            m_data,
            cmap="viridis",
            shading="auto",
        )
        ax.grid(True)

        # # Add x-ticks every 1 unit
        # xticks = np.arange(section["distance"].min(), section["distance"].max() + 1, 1)
        # ax.set_xticks(xticks)

        # # Orientation verticale des étiquettes
        # ax.tick_params(axis="x", labelrotation=90)

        # Add a black line for topo section
        topo_profile = section["topo"]
        distance_profile = section["distance"]
        # Plot the depth profile as a black dashed line
        ax.plot(
            distance_profile,
            topo_profile,
            color="black",
            linestyle="--",
            linewidth=1.5,
        )

        # Titre et labels
        ax.set_title(f"{titles.get(var, var)} - Along {direction}")
        if direction == "x":
            ax.set_xlabel("Along Longitude")
            ax.set_ylabel("Depth (m)")
        elif direction == "y":
            ax.set_xlabel("Along Latitude")
            ax.set_ylabel("Depth (m)")
        else:
            ax.set_xlabel("Longitude/Latitude")
            ax.set_ylabel("Depth (m)")

        # Ajouter une barre de couleur
        fig.colorbar(im, ax=ax, orientation="vertical", label=units.get(var, ""))

        # Mettre à jour l'indice de l'axe
        axis_index += 1

# Supprimer les axes inutilisés si le nombre de variables est inférieur à la grille
for j in range(n_dirs, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()

# Titre général
plt.suptitle(
    f"Vertical Sections for Variables at J/I index {ind_sec}", fontsize=16, y=1.02
)

# plt.show()
# plt.savefig('REGIONAL_sections.png', bbox_inches='tight')

#################################
# Determine output file paths
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(output_dir, "REGIONAL_sections.pdf")
    plt.savefig(pdf_path, transparent=True, bbox_inches="tight")
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(output_dir, "REGIONAL_sections.png")
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"PNG file '{png_path}' has been created.")

# Show plots if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --show to enable).")
