#!/usr/bin/env python3

from abc import update_abstractmethods
from asyncio.subprocess import SubprocessStreamProtocol
import os
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import argparse

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker

from matplotlib.colors import Normalize
from matplotlib.ticker import MaxNLocator

import math


#################################################################
def read_latlonmask(fname, var_type):
    """
    Lit les latitudes, longitudes et le masque depuis un fichier NetCDF CROCO.

    Args:
        fname (str): Nom du fichier NetCDF CROCO.
        var_type (str): Type de variable ('r' pour rho, 'w', 'u', ou 'v').

    Returns:
        lat (np.ndarray): Latitude (matrice 2D).
        lon (np.ndarray): Longitude (matrice 2D).
        mask (np.ndarray): Masque (1 pour la mer, NaN pour la terre).
    """
    # Ouverture du fichier NetCDF
    nc = Dataset(fname)
    lat = nc.variables["lat_rho"][:]
    lon = nc.variables["lon_rho"][:]
    mask = nc.variables.get("mask_rho", np.ones_like(lon))

    # Gestion des types spécifiques
    if var_type == "u":
        lat = rho2u_2d(lat)
        lon = rho2u_2d(lon)
        mask = mask[:, :-1] * mask[:, 1:]
    elif var_type == "v":
        lat = rho2v_2d(lat)
        lon = rho2v_2d(lon)
        mask = mask[:-1, :] * mask[1:, :]

    # Conversion du masque : 0 -> NaN
    mask = np.where(mask == 0, np.nan, mask)

    return lat, lon, mask


#################################################################
def rho2u_2d(var):
    """
    Convertit une variable rho-grid vers u-grid.
    """
    return 0.5 * (var[:, :-1] + var[:, 1:])


#################################################################
def rho2v_2d(var):
    """
    Convertit une variable rho-grid vers v-grid.
    """
    return 0.5 * (var[:-1, :] + var[1:, :])


#################################################################
def rho2u_3d(var):
    return 0.5 * (var[:, :, :-1] + var[:, :, 1:])


#################################################################
def rho2v_3d(var):
    return 0.5 * (var[:, :-1, :] + var[:, 1:, :])


#################################################################
def get_csf(sc, theta_s, theta_b):
    """
    function CSF = get_csf(sc,theta_s,theta_b) Get CS-Curves for the new s-ccordinate system
    NOTE: Mathematical limits of CSF,csrf for theta_s, theta_b --> 0 match that under "else" logical branches.
    """

    if theta_s > 0.0:
        csrf = (1.0 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - 1.0)
    else:
        csrf = -(sc**2)

    if theta_b > 0.0:
        CSF = (np.exp(theta_b * csrf) - 1.0) / (1.0 - np.exp(-theta_b))
    else:
        CSF = csrf

    return CSF


###################################################
def scoordinate(theta_s, theta_b, N, hc, vtransform):
    """
    function [sc_r,Cs_r,sc_w,Cs_w] = scoordinate(theta_s,theta_b,N,hc,vtransform)
    Define S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    """

    #  sc_r = np.zeros(N)
    #  Cs_r = np.zeros(N)
    #  sc_w = np.zeros(N+1)
    #  Cs_w = np.zeros(N+1)

    #
    if vtransform == 2:
        # print("NEW_S_COORD")

        ds = 1.0 / N

        sc_r = ds * (np.arange(1, N + 1) - N - 0.5)
        Cs_r = get_csf(sc_r, theta_s, theta_b)
        #
        #    sc_w[0]   = -1.
        #    sc_w[N]   =  0.
        #    Cs_w[0]   = -1.
        #    Cs_w[N]   =  0.

        sc_w = ds * (np.arange(0, N + 1) - N)
        Cs_w = get_csf(sc_w, theta_s, theta_b)

    else:

        print("OLD_S_COORD")

        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)

        sc_w = (np.arange(0, N + 1) - N) / N
        Cs_w = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc_w) + theta_b * (
            cff2 * np.tanh(theta_s * (sc_w + 0.5)) - 0.5
        )

        sc_r = (np.arange(1, N + 1) - N - 0.5) / N
        Cs_r = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc_r) + theta_b * (
            cff2 * np.tanh(theta_s * (sc_r + 0.5)) - 0.5
        )

    return sc_r, Cs_r, sc_w, Cs_w


###################################################
def zlevs(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    function z=zlevs(h,zeta,theta_s,theta_b,hc,N,type,vtransform)

     this function compute the depth of rho or w points for CROCO

     On Input:

       type    'r': rho point 'w': w point
       vtransform  1=> old v transform (Song, 1994);
                   2=> new v transform (Shcheptekin, 2006)
     On Output:

       z       Depths (m) of RHO- or W-points (3D matrix).
    """

    #
    # Test the number of dimension for h
    #

    Ndim = np.size(np.shape(h))
    if Ndim == 2:
        M, L = np.squeeze(np.shape(h))
    elif Ndim == 1:
        L = np.squeeze(np.shape(h))
    else:
        print("zlevs: error - incorrect dimension for h")
        return

    hmin = h.min()
    hmax = h.max()

    #
    # Set S-Curves in domain [-1 < sc < 0] at vertical W- and RHO-points.
    #

    ds = 1.0 / float(N)

    if type == "w":

        sc = ds * (np.arange(0, N + 1) - N)
        Nmax = N + 1

    elif type == "r":

        sc = ds * (np.arange(1, N + 1) - N - 0.5)
        Nmax = N

    else:

        print("Problem with type = ", type)
        sys.exit()

    if vtransform == 1:

        print("OLD_S_COORD")

        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5
        )

    elif vtransform == 2:
        # print("NEW_S_COORD")
        Cs = get_csf(sc, theta_s, theta_b)

    else:
        print("Problem with vtransform = ", vtransform)
        sys.exit()

    #
    #  Set vertical grid
    #

    if Ndim == 2:
        z = np.zeros([Nmax, M, L])
    elif Ndim == 1:
        z = np.zeros([Nmax, L])

    for k in np.arange(0, Nmax):

        if vtransform == 1:

            cff = hc * (sc[k] - Cs[k])
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / h
            if Ndim == 2:
                z[k, :, :] = z0 + zeta * (1.0 + z0 * hinv)
            elif Ndim == 1:
                z[k, :] = z0 + zeta * (1.0 + z0 * hinv)

        elif vtransform == 2:

            cff = hc * sc[k]
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / (h + hc)
            if Ndim == 2:
                z[k, :, :] = z0 * h * hinv + zeta * (1.0 + z0 * hinv)
            elif Ndim == 1:
                z[k, :] = z0 * h * hinv + zeta * (1.0 + z0 * hinv)

    return z


#################################################################
def rempoints(var, npts):
    """
    Réduit une matrice 2D en supprimant des points aux frontières.

    Args:
        var (np.ndarray): Données à traiter (matrice 2D).
        npts (list[int]): Nombre de points à supprimer sous la forme [west, east, south, north].

    Returns:
        np.ndarray: Matrice 2D après suppression des points de bordure.
    """
    # Dimensions initiales
    M, L = var.shape

    # Extraire les indices pour chaque direction
    west, east, south, north = npts

    # Réduire les dimensions en fonction des points à supprimer
    var = var[south : M - north, west : L - east]

    return var


#################################################################
def get_type(fname, vname, vlevin):
    """
    Détermine le "type" d'une variable CROCO : rho, u ou v.

    Args:
        fname (str): Nom du fichier NetCDF.
        vname (str): Nom de la variable à examiner.
        vlevin (int): Niveau vertical d'entrée.

    Returns:
        tuple:
            - type (str): Type de la variable ('r', 'u', 'v', 'w', ou '').
            - vlevout (int): Niveau vertical de sortie (inchangé ou modifié).
    """
    vlevout = vlevin
    var_type = "r"  # Par défaut, type est 'r'

    try:
        # Ouverture du fichier NetCDF
        with Dataset(fname, "r") as nc:
            if vname not in nc.variables:
                return "", vlevout

            var = nc.variables[vname]
            dims = var.dimensions
            ndim = len(dims)

            if ndim == 1:  # Cas d'une variable 1D
                return "", vlevout

            i = 0
            name = dims[i]

            # Vérifier si la dimension est temporelle
            if not (name.endswith("time") or name.startswith("time")):
                print("Avertissement : Pas dépendant du temps.")
            else:
                i += 1

            # Vérifier si la dimension est verticale
            name = dims[i]
            if name.startswith("s"):
                if name.endswith("w"):
                    return "w", vlevout
                else:
                    i += 1
            else:
                vlevout = 0

            # Vérifier si la dimension est spatiale (lat/y ou lon/x)
            name = dims[i]
            if not (name.startswith("e") or name.startswith("y")):
                return "", vlevout
            else:
                if name.endswith("v"):
                    return "v", vlevout
                if name.endswith("u"):
                    return "u", vlevout

            # Vérifier une deuxième dimension spatiale
            if i + 1 < ndim:
                name = dims[i + 1]
                if name.startswith("x"):
                    if name.endswith("u"):
                        return "u", vlevout
                    if name.endswith("v"):
                        return "v", vlevout

    except Exception as e:
        print(f"Erreur lors de la lecture du fichier NetCDF : {e}")
        return "", vlevout

    return var_type, vlevout


#################################################################
def vinterp(var, z, depth):
    """
    Interpole une variable 3D sur un niveau horizontal de profondeur constante.

    Args:
        var (numpy.ndarray): Variable à traiter (matrice 3D de taille [N, Mp, Lp]).
        z (numpy.ndarray): Profondeurs (m) des points RHO ou W (matrice 3D de taille [N, Mp, Lp]).
        depth (float): Profondeur de la coupe (scalaire; mètres, valeur négative).

    Returns:
        numpy.ndarray: Coupe horizontale interpolée (matrice 2D de taille [Mp, Lp]).
    """
    N, Mp, Lp = z.shape

    # Trouver la position des niveaux verticaux les plus proches
    below = z < depth
    levs = np.sum(below, axis=0)  # Nombre de niveaux sous la profondeur spécifiée
    levs[levs == N] = (
        N - 1
    )  # Limiter à l'avant-dernier niveau pour éviter les débordements

    # Créer un masque pour identifier les indices valides
    mask = np.where(levs > 0, 1, np.nan)

    # Générer les positions linéaires pour un accès rapide
    imat, jmat = np.meshgrid(np.arange(Lp), np.arange(Mp), indexing="ij")
    pos = N * Mp * imat + N * jmat + levs
    pos = pos.astype(int)  # Assurer que les indices sont des entiers
    pos[np.isnan(mask)] = 0  # Éviter les erreurs sur les valeurs NaN

    # Effectuer l'interpolation
    z1 = z.ravel()[pos + 1]
    z2 = z.ravel()[pos]
    v1 = var.ravel()[pos + 1]
    v2 = var.ravel()[pos]

    # Calculer la variable interpolée
    vnew = mask * (((v1 - v2) * depth + v2 * z1 - v1 * z2) / (z1 - z2))

    return vnew.reshape(Mp, Lp)


#################################################################
def get_depths(fname, gname, tindex, point_type):
    """
    Calcule les profondeurs des niveaux sigma.

    Args:
        fname (str): Chemin du fichier historique.
        gname (str): Chemin du fichier grille.
        tindex (int): Indice temporel.
        point_type (str): Type de point ('rho', 'u', 'v', 'w').
    Returns:
        numpy.ndarray: Profondeurs des niveaux sigma (matrice 3D).
    """

    # Lecture des données de la grille
    with Dataset(gname, mode="r") as nc:
        h = nc.variables["h"][:]

    # Lecture des données historiques
    with Dataset(fname, mode="r") as nc:
        zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
        hmorph = nc.variables.get("hmorph", None)
        if hmorph is not None:
            h = np.squeeze(hmorph)

        # Lecture des paramètres de verticalité
        try:
            theta_s = nc.variables["theta_s"][:]
            theta_b = nc.variables["theta_b"][:]
            Tcline = nc.variables["Tcline"][:]
            hc = nc.variables.get("hc", None)
        except KeyError:
            theta_s = nc.theta_s
            theta_b = nc.theta_b
            Tcline = nc.Tcline
            hc = getattr(nc, "hc", None)

        if Tcline is not None:
            hmin = np.min(h)
            hc = min(hmin, Tcline)

        N = len(nc.dimensions["s_rho"])
        VertCoordType = getattr(nc, "VertCoordType", "")
        vtrans = np.squeeze(nc.variables.get("Vtransform", None))
        s_coord = 2 if VertCoordType == "NEW" or vtrans == 2 else 1

        if s_coord == 2:
            hc = Tcline

    if zeta is None:
        zeta = np.zeros_like(h)

    # Calcul des profondeurs des niveaux sigma
    vtype = point_type
    if point_type in ('u', 'v'):
        vtype = 'r'
        
    z = zlevs(h, zeta, theta_s, theta_b, hc, N, vtype, vtrans)

    # Convertir les profondeurs si type est 'u' ou 'v'
    if point_type == "u":
        z = rho2u_3d(z)
    elif point_type == "v":
        z = rho2v_3d(z)

    return z


#################################################################
def get_hslice(fname, gname, vname, tindex, level, var_type):
    """
    Extrait une coupe horizontale d'une variable CROCO.

    Args:
        fname (str): Chemin du fichier CROCO (average ou history).
        gname (str): Chemin du fichier de grille CROCO.
        vname (str): Nom de la variable.
        tindex (int): Indice temporel.
        level (float): Niveau vertical pour la coupe :
                       - level = 0 : variable 2D horizontale (ex. zeta)
                       - level > 0 : coupe au niveau sigma (entier >= 1)
                       - level < 0 : interpolation à z = level.
        var_type (str): Type de la variable ('r', 'u', 'v', ou 'w').

    Returns:
        numpy.ndarray: Coupe horizontale (matrice 2D).
    """

    with Dataset(fname, mode="r") as nc:
        if level == 0:
            # Variable 2D
            var = np.squeeze(nc.variables[vname][tindex, :, :])
            var[var == 0] = np.nan

            # Correction pour le mouillage/séchage
            with Dataset(gname, mode="r") as ng:
                h = ng.variables["h"][:]
            hmorph = nc.variables.get("hmorph", None)
            if hmorph is not None:
                h = np.squeeze(hmorph[tindex, :, :])
            zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
            D = zeta + h

            Dcrit = nc.variables.get("Dcrit", None)
            if Dcrit is not None:
                Dcrit = Dcrit[:] + 1e-5
            else:
                Dcrit = 0.2 + 1e-5

            if var_type == "u":
                D = rho2u_2d(D)
            elif var_type == "v":
                D = rho2v_2d(D)

            var[D <= Dcrit] = np.nan

        elif level > 0:
            # Coupe à un niveau sigma spécifique
            var = np.squeeze(nc.variables[vname][tindex, level - 1, :, :])
            var[var == 0] = np.nan

        else:
            # Coupe horizontale à un niveau de profondeur donné
            # Obtenir les profondeurs des niveaux sigma
            z = get_depths(fname, gname, tindex, var_type)

            # Lire la variable 3D et interpoler
            var_sigma = np.squeeze(nc.variables[vname][tindex, :, :, :])
            var = vinterp(var_sigma, z, level)

    return var


#################################################################
def get_var(hisfile, gridfile, vname, tindex, vlevel, coef, rempts):
    """
    Récupère une variable depuis un fichier CROCO NetCDF.

    Args:
        hisfile (str): Nom du fichier historique NetCDF.
        gridfile (str): Nom du fichier de grille NetCDF.
        vname (str): Nom de la variable à récupérer.
        tindex (int): Index temporel de la variable.
        vlevel (float): Niveau vertical.
        coef (float): Coefficient multiplicatif.
        rempts (int): Points à supprimer aux bords.

    Returns:
        tuple: (lat, lon, mask, var) où:
            lat (np.ndarray): Latitudes.
            lon (np.ndarray): Longitudes.
            mask (np.ndarray): Masque (1 = mer, NaN = terre).
            var (np.ndarray): Variable demandée.
    """
    print(f"-- {vname} --")

    if not gridfile:
        gridfile = hisfile
    if vname in ["h", "pm", "pn"]:
        lat, lon, mask = read_latlonmask(gridfile, "r")
        with Dataset(gridfile, "r") as nc:
            var = nc.variables[vname][:]

        mask = np.where(mask == 0, np.nan, mask)
        var = mask * var

    else:
        var_type, vlev = get_type(hisfile, vname, vlevel)
        lat, lon, mask = read_latlonmask(gridfile, var_type)
        var = coef * mask * get_hslice(hisfile, gridfile, vname, tindex, vlev, var_type)

    return lat, lon, mask, var


#################################################################
def do_plot(lon, lat, var, colmin, colmax, ncol, hisfile, gridfile, tindex, vlevel):
    """
    Do the plot
    """
    # Pseudocolor plot
    plt.figure()

    ax = plt.axes(projection=ccrs.PlateCarree())
    mesh = ax.pcolormesh(
        lon,
        lat,
        var,
        shading="flat",
        cmap="jet",
        norm=Normalize(vmin=colmin, vmax=colmax),
    )
    plt.colorbar(mesh, ax=ax)

    plt.show()


#################################################################
def get_vertical_section(fname, vname, tindex, direction, idx, gname):
    """
    Extrait une section verticale d'une variable 3D.

    Args:
        fname (str): Nom du fichier historique NetCDF.
        vname (str): Nom de la variable (par exemple, "temp").
        tindex (int): Index temporel.
        direction (str): Direction de la section ("x" pour zonale ou "y" pour méridienne).
        idx (int): Indice de la ligne ou colonne pour la coupe.
        gname (str): Nom du fichier de grille NetCDF.

    Returns:
        dict: Contient les profondeurs, 
                la distance (number of grid point), 
                et les valeurs interpolées.
    """

    with Dataset(fname, "r") as nc:
        # Déterminer le type de grille pour obtenir les profondeurs
        if vname == "w":
            point_type = "w"  # Grille régulière pour u
        elif vname == "u": 
            point_type = "u"  # Grille régulière pour u
        elif vname == "v": 
            point_type = "v"  # Grille régulière pour v
        else:
            point_type = "r"  # Grille régulière pour les scalaires (temp, salt, etc.)
            
        # Lire les profondeurs        
        z = get_depths(fname, gname, tindex, point_type)
        temp = np.squeeze(nc.variables[vname][tindex, :, :, :])

    if direction == "x":
        section_data = temp[:, :, idx]
        depth_section = z[:, :, idx]
        distance = np.arange(section_data.shape[1])  # Distance en points de grille
    elif direction == "y":
        section_data = temp[:, idx, :]
        depth_section = z[:, idx, :]
        distance = np.arange(section_data.shape[1])  # Distance en points de grille
    else:
        raise ValueError("Direction invalide. Choisissez 'x' ou 'y'.")

    return {
        "depth": depth_section,
        "distance": distance,
        "variable": section_data,
    }


#################################################################
def plot_vertical_section(section_data, title="Vertical Section", unit="°C"):
    """
    Trace une section verticale.

    Args:
        section_data (dict): Contient les profondeurs, la distance et les valeurs.
        title (str): Titre du graphique.
        unit (str): Unité de la variable.
    """
    depth = section_data["depth"]
    distance = section_data["distance"]
    variable = section_data["variable"]

    plt.figure(figsize=(12, 6))
    # plt.contourf(distance, -depth, variable, levels=50, cmap="viridis")
    plt.pcolormesh(distance, depth, variable, cmap="viridis",shading='auto')
    plt.colorbar(label=unit)
    plt.title(title)
    plt.xlabel("Distance (grille)")
    plt.ylabel("Profondeur (m)")
    plt.grid(True)
    #plt.show()


#################################################################
####### MAIN MAIN #######
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
lat, lon, mask = read_latlonmask(hisfile, "r")
latu, lonu, masku = read_latlonmask(hisfile, "u")
latv, lonv, maskv = read_latlonmask(hisfile, "v")

# get the number of vertical level
nc = Dataset(hisfile)
temp = nc.variables["temp"][:]
T = np.shape(temp)[0]
N = np.shape(temp)[1]

# choose type of plots
plot_horiz = 1
plot_section = 1

# FOR THE HORIZONTAL PLOT
if plot_horiz == 1:
    #
    # read the horizslice : last time step , surface level
    tndx = -1
    vlev = N
    # Read var
    ssh = get_hslice(hisfile, hisfile, "zeta", tndx, 0, "r")
    sst = get_hslice(hisfile, hisfile, "temp", tndx, vlev, "r")
    sss = get_hslice(hisfile, hisfile, "salt", tndx, vlev, "r")
    ubar = get_hslice(hisfile, hisfile, "ubar", tndx, 0, "u")
    vbar = get_hslice(hisfile, hisfile, "vbar", tndx, 0, "v")
    usurf = get_hslice(hisfile, hisfile, "u", tndx, vlev, "u")
    vsurf = get_hslice(hisfile, hisfile, "v", tndx, vlev, "v")
    sustr = get_hslice(hisfile, hisfile, "sustr", tndx, 0, "u")
    svstr = get_hslice(hisfile, hisfile, "svstr", tndx, 0, "v")

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
    #plt.show()
    #plt.savefig('REGIONAL_maps.png')
    
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
    idx = 20  # Indice de la colonne pour la coupe

    # Extraire les données
    var_sec = "temp"
    section = get_vertical_section(hisfile, var_sec, tindex, direction, idx, gridfile)

    # Tracer la section verticale
    plot_vertical_section(
        section, title=f"Zonal Vertical section -  {var_sec} - j={idx}", unit="°C"
    )

# Exemple de dictionnaire avec variables et leurs directions
variables = {
    "temp": ["x", "y"],  # Coupe zonale et méridienne pour la température
    "salt": ["x", "y"],  # Coupe zonale et méridienne pour la salinité
    "u": ["x", "y"],     # Zonal et méridien pour la vitesse zonale
    "v": ["x", "y"],     # Zonal et méridien pour la vitesse méridienne
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

# Paramètres pour les sections verticales
gridfile = hisfile  # Fichier de la grille
tindex = -1  # Dernier pas de temps
idx = 20  # Indice de la colonne pour la coupe

# Calculer la disposition des sous-graphiques (grille optimale)
n_vars = len(variables)
n_dirs = sum(len(directions) for directions in variables.values())  # Nombre total de sections
n_cols = 2  # Nombre de colonnes fixes
n_rows = math.ceil(n_dirs / n_cols)  # Nombre de lignes nécessaires

# Initialisation des subplots
fig, axes = plt.subplots(
    n_rows, n_cols, figsize=(15, n_rows * 3))#, constrained_layout=True

# Aplatir les axes pour gestion facile (en cas de grille)
axes = axes.flatten()

# Initialisation de l'indice pour les axes
axis_index = 0

# Boucle pour parcourir chaque variable et ses directions
for var, directions in variables.items():
    for direction in directions:
        # Extraire les données de la section
        print(f"Processing {var} in direction {direction}")
        section = get_vertical_section(hisfile, var, tindex, direction, idx, gridfile)

        # Configurer les axes
        ax = axes[axis_index]
        im = ax.pcolor(
            section["distance"],
            section["depth"],
            section["variable"],
            cmap="viridis",
            shading="auto",
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
plt.suptitle(f"Vertical Sections for Variables at J/I index {idx}", fontsize=16, y=1.02)

#plt.show()
#plt.savefig('REGIONAL_sections.png', bbox_inches='tight')

#################################
# Determine output file paths
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

# Save to PDF if requested
if args.makepdf:
    pdf_path = os.path.join(output_dir, "REGIONAL_sections.pdf")
    plt.savefig(pdf_path, transparent=True, bbox_inches='tight')
    print(f"PDF file '{pdf_path}' has been created.")

# Save to PNG if requested
if args.makepng:
    png_path = os.path.join(output_dir, "REGIONAL_sections.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"PNG file '{png_path}' has been created.")

# Show plots if not suppressed
if not args.no_show:
    plt.show()
else:
    print("Plot display suppressed (use --show to enable).")