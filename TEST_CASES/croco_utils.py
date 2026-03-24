#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


# ──────────────────────────────────────────────────────────
#  Internal building blocks
# ──────────────────────────────────────────────────────────

def _compute_sc(N, point_type):
    """
    Compute raw S-coordinate values in [-1, 0].

    Parameters
    ----------
    N : int
        Number of vertical rho-levels.
    point_type : str
        'r' for rho points, 'w' for w points.

    Returns
    -------
    sc : np.ndarray
        S-coordinate values.  Shape (N,) for rho, (N+1,) for w.
    """
    ds = 1.0 / N
    if point_type == "w":
        return ds * (np.arange(0, N + 1) - N)
    elif point_type == "r":
        return ds * (np.arange(1, N + 1) - N - 0.5)
    else:
        raise ValueError(f"Unknown point_type '{point_type}', expected 'r' or 'w'")


def _compute_cs(sc, theta_s, theta_b, vtransform):
    """
    Compute the CS stretching function.

    Parameters
    ----------
    sc : np.ndarray
        S-coordinate values (from _compute_sc).
    theta_s : float
        Surface control parameter.
    theta_b : float
        Bottom control parameter.
    vtransform : int
        1 = old (Song & Haidvogel 1994), 2 = new (Shchepetkin 2005).

    Returns
    -------
    Cs : np.ndarray
        Stretching function, same shape as *sc*.
    """
    if vtransform == 2:
        # New S-coordinate: delegate to the two-step Cs formula
        if theta_s > 0.0:
            csrf = (1.0 - np.cosh(theta_s * sc)) / (np.cosh(theta_s) - 1.0)
        else:
            csrf = -(sc ** 2)

        if theta_b > 0.0:
            Cs = (np.exp(theta_b * csrf) - 1.0) / (1.0 - np.exp(-theta_b))
        else:
            Cs = csrf

    elif vtransform == 1:
        # Old S-coordinate: sinh/tanh stretching
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = ((1.0 - theta_b) * cff1 * np.sinh(theta_s * sc)
              + theta_b * (cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5))

    else:
        raise ValueError(f"Unknown vtransform={vtransform}, expected 1 or 2")

    return Cs


def _compute_z(h, zeta, sc, Cs, hc, vtransform):
    """
    Compute 3D depth field from S-coordinate parameters.

    Parameters
    ----------
    h : np.ndarray
        Bathymetry, shape (M, L).
    zeta : np.ndarray
        Free surface elevation, shape (M, L).
    sc : np.ndarray
        S-coordinate values, shape (Nk,).
    Cs : np.ndarray
        Stretching function, shape (Nk,).
    hc : float
        Critical depth parameter.
    vtransform : int
        Vertical transformation equation (1 or 2).

    Returns
    -------
    z : np.ndarray
        Depths, shape (Nk, M, L).
    """
    Nk = len(sc)
    M, L = h.shape
    z = np.empty((Nk, M, L))

    if vtransform == 1:
        hinv = 1.0 / h
        for k in range(Nk):
            z0 = hc * (sc[k] - Cs[k]) + Cs[k] * h
            z[k, :, :] = z0 + zeta * (1.0 + z0 * hinv)

    elif vtransform == 2:
        hinv = 1.0 / (h + hc)
        for k in range(Nk):
            z0 = hc * sc[k] + Cs[k] * h
            z[k, :, :] = z0 * h * hinv + zeta * (1.0 + z0 * hinv)

    return z


# ──────────────────────────────────────────────────────────
#  Public API (signatures preserved for backward compat)
# ──────────────────────────────────────────────────────────

def get_csf(sc, theta_s, theta_b):
    """
    CS stretching function for the new S-coordinate (vtransform=2).

    Kept for backward compatibility.  Equivalent to
    ``_compute_cs(sc, theta_s, theta_b, vtransform=2)``.
    """
    return _compute_cs(sc, theta_s, theta_b, vtransform=2)


def scoordinate(theta_s, theta_b, N, hc, vtransform):
    """
    Compute S-coordinate values and stretching at both rho and w points.

    Parameters
    ----------
    theta_s, theta_b : float
        Stretching parameters.
    N : int
        Number of rho-levels.
    hc : float
        Critical depth (unused here, kept for interface compatibility).
    vtransform : int
        1 (old) or 2 (new).

    Returns
    -------
    sc_r, Cs_r : np.ndarray
        S-coordinate and stretching at rho points, shape (N,).
    sc_w, Cs_w : np.ndarray
        S-coordinate and stretching at w points, shape (N+1,).
    """
    sc_r = _compute_sc(N, "r")
    sc_w = _compute_sc(N, "w")
    Cs_r = _compute_cs(sc_r, theta_s, theta_b, vtransform)
    Cs_w = _compute_cs(sc_w, theta_s, theta_b, vtransform)
    return sc_r, Cs_r, sc_w, Cs_w


def zlevs(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    Compute the depths of rho or w points for CROCO.

    Parameters
    ----------
    h : np.ndarray
        Bathymetry (1D or 2D).
    zeta : np.ndarray
        Free surface elevation (same shape as h).
    theta_s, theta_b : float
        Stretching parameters.
    hc : float
        Critical depth parameter.
    N : int
        Number of rho-levels.
    type : str
        'r' for rho points, 'w' for w points.
    vtransform : int
        Vertical transformation equation (1 or 2).

    Returns
    -------
    z : np.ndarray
        Depths (3D squeezed if h was 1D).
    """
    # Handle 1D input
    ndim = np.ndim(h)
    if ndim == 1:
        h = h[:, np.newaxis]
        zeta = zeta[:, np.newaxis]
    elif ndim != 2:
        raise ValueError("zlevs: h must be 1D or 2D")

    sc = _compute_sc(N, type)
    Cs = _compute_cs(sc, theta_s, theta_b, vtransform)
    z = _compute_z(h, zeta, sc, Cs, hc, vtransform)

    return z.squeeze()


################################################################
def rho2u_2d(rho_in):
    """
    Convert a 2D field at rho points to a field at u points
    """

    def _r2u(rho_in, Lp):
        u_out = 0.5 * (rho_in[:, :Lp - 1] + rho_in[:, 1:Lp])
        return u_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2u(rho_in, Lshp)


#################################################################
def rho2u_3d(var):
    """
    Convert 3D variable from rho-grid to u-grid.
    """
    return 0.5 * (var[:, :, :-1] + var[:, :, 1:])


################################################################
def rho2v_2d(rho_in):
    """
    Convert a 2D field at rho points to a field at v points
    """

    def _r2v(rho_in, Mp):
        v_out = 0.5 * (rho_in[:Mp - 1] + rho_in[1:Mp])
        return v_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2v(rho_in, Mshp)


################################################################
def u2rho_2d(u_in):
    """
    Convert a 2D field at u points to a field at rho points
    """

    def _uu2ur(uu_in, Mp, Lp):
        L, Lm = Lp - 1, Lp - 2
        u_out = np.zeros((Mp, Lp))
        u_out[:, 1:L] = 0.5 * (u_in[:, 0:Lm] + u_in[:, 1:L])
        u_out[:, 0] = u_out[:, 1]
        u_out[:, L] = u_out[:, Lm]
        return u_out.squeeze()

    assert u_in.ndim == 2, "u_in must be 2d"
    Mp, Lp = u_in.shape
    return _uu2ur(u_in, Mp, Lp + 1)


################################################################
def v2rho_2d(v_in):
    """
    Convert a 2D field at v points to a field at rho points
    """

    def _vv2vr(v_in, Mp, Lp):
        M, Mm = Mp - 1, Mp - 2
        v_out = np.zeros((Mp, Lp))
        v_out[1:M] = 0.5 * (v_in[:Mm] + v_in[1:M])
        v_out[0] = v_out[1]
        v_out[M] = v_out[Mm]
        return v_out.squeeze()

    assert v_in.ndim == 2, "v_in must be 2d"
    Mp, Lp = v_in.shape
    return _vv2vr(v_in, Mp + 1, Lp)


#################################################################
def rho2v_3d(var):
    """
    Convert 3D variable from rho-grid to v-grid.
    """
    return 0.5 * (var[:, :-1, :] + var[:, 1:, :])


################################################################
def rho2uvp(rfield):
    """
    Convert 3D variable from rho-grid to uvp-grid.
    """
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1

    vfield = 0.5 * (rfield[np.arange(0, M), :] + rfield[np.arange(1, Mp), :])
    ufield = 0.5 * (rfield[:, np.arange(0, L)] + rfield[:, np.arange(1, Lp)])
    pfield = 0.5 * (ufield[np.arange(0, M), :] + ufield[np.arange(1, Mp), :])

    return ufield, vfield, pfield


################################################################
def read_latlonmask(fname, var_type):
    """
    Read latitudes and longitudes and mask from CROCO netcdf file.

    Args:
        fname (str): CROCO netcdf file
        var_type (str): variable type ('r' for rho, 'w', 'u', or 'v').

    Returns:
        lat (np.ndarray): Latitude (matrix 2D).
        lon (np.ndarray): Longitude (matrix 2D).
        mask (np.ndarray): Mask (1 for sea, NaN for land).
    """
    with Dataset(fname) as nc:
        lat = nc.variables["lat_rho"][:]
        lon = nc.variables["lon_rho"][:]
        mask = (nc.variables["mask_rho"][:] if "mask_rho" in nc.variables
                else np.ones_like(lon))

    # Manage specific types
    if var_type == "u":
        lat = rho2u_2d(lat)
        lon = rho2u_2d(lon)
        mask = mask[:, :-1] * mask[:, 1:]
    elif var_type == "v":
        lat = rho2v_2d(lat)
        lon = rho2v_2d(lon)
        mask = mask[:-1, :] * mask[1:, :]

    # Convert mask : 0 -> NaN
    mask = np.where(mask == 0, np.nan, mask)

    return lat, lon, mask


#################################################################
def rempoints(var, npts):
    """
    Reduce a 2D matrix by removing border points

    Args:
        var (np.ndarray): input 2D matrix
        npts (list[int]): number of points to remove at each border
                          with format [west, east, south, north].

    Returns:
        np.ndarray: 2D matrix after removing border points
    """
    # Initial dimensions
    M, L = var.shape

    # Extract index from input
    west, east, south, north = npts

    # Reduce dimensions from input number of point to remove
    var = var[south : M - north, west : L - east]

    return var


#################################################################
def get_type(fname, vname, vlevin):
    """
    Determines the grid type of a CROCO variable from its dimension names.

    The grid type is inferred from the CROCO dimension naming convention:
        eta_rho / xi_rho -> 'r',  eta_u / xi_u -> 'u',
        eta_v / xi_v -> 'v',      s_w -> 'w'.

    Args:
        fname (str): Name of the NetCDF file.
        vname (str): Name of the variable to examine.
        vlevin (int): Input vertical level.

    Returns:
        tuple:
            - var_type (str): Grid type ('r', 'u', 'v', 'w', or '').
            - vlevout (int): Output vertical level (set to 0 if no vertical dim).
    """
    vlevout = vlevin

    try:
        with Dataset(fname, "r") as nc:
            if vname not in nc.variables:
                return "", vlevout

            dims = nc.variables[vname].dimensions

            if len(dims) <= 1:
                return "", vlevout

            # Detect vertical dimension
            has_vertical = any(d.startswith("s_") for d in dims)
            if not has_vertical:
                vlevout = 0

            # Check for w-level vertical coordinate
            if any(d == "s_w" for d in dims):
                return "w", vlevout

            # Detect horizontal grid type from dimension names
            # CROCO uses: eta_u/xi_u for u-grid, eta_v/xi_v for v-grid
            dim_set = set(dims)
            if "xi_u" in dim_set or "eta_u" in dim_set:
                return "u", vlevout
            if "xi_v" in dim_set or "eta_v" in dim_set:
                return "v", vlevout

            # Fallback: check for rho-grid or generic dimensions
            if any(d.startswith("eta") or d.startswith("xi") or
                   d.startswith("x") or d.startswith("y") or
                   d.startswith("e") for d in dims):
                return "r", vlevout

            return "", vlevout

    except Exception as e:
        print(f"Error reading NetCDF file '{fname}': {e}")
        return "", vlevout


#################################################################
def vinterp(var, z, depth):
    """
    Interpolates a 3D variable onto a horizontal level of constant depth.

    Args:
        var (numpy.ndarray): Variable to process (3D matrix of size [N, Mp, Lp]).
        z (numpy.ndarray): Depths (m) of RHO or W points (3D matrix of size [N, Mp, Lp]).
        depth (float): Depth of the slice (scalar; meters, negative value).

    Returns:
        numpy.ndarray: Interpolated horizontal slice (2D matrix of size [Mp, Lp]).
    """

    N, Mp, Lp = z.shape

    # Find the number of levels below the requested depth
    below = z < depth
    levs = np.sum(below, axis=0)  # shape (Mp, Lp)

    # Clamp to [1, N-1] so that both levs-1 and levs are valid indices in [0, N-1]
    levs = np.clip(levs, 1, N - 1)

    # Create a mask: NaN where the depth is entirely above or below the grid
    mask = np.where((np.sum(below, axis=0) > 0) & (np.sum(below, axis=0) < N),
                    1.0, np.nan)

    # Build index arrays for the two bracketing levels
    jmat, imat = np.meshgrid(np.arange(Mp), np.arange(Lp), indexing="ij")
    # C-order ravel: z[k, j, i] -> k * Mp * Lp + j * Lp + i
    stride = Mp * Lp
    pos_lo = (levs - 1) * stride + jmat * Lp + imat   # level below
    pos_hi = levs * stride + jmat * Lp + imat          # level above

    z_lo = z.ravel()[pos_lo]
    z_hi = z.ravel()[pos_hi]
    v_lo = var.ravel()[pos_lo]
    v_hi = var.ravel()[pos_hi]

    # Linear interpolation
    vnew = mask * (v_lo + (v_hi - v_lo) * (depth - z_lo) / (z_hi - z_lo))

    return vnew


#################################################################
def get_depths(fname, gname, tindex, point_type):
    """
    Calculates the depths of sigma levels.

    Args:
        fname (str): Path to the history file.
        gname (str): Path to the grid file.
        tindex (int): Time index.
        point_type (str): Point type ('rho', 'u', 'v', 'w').

    Returns:
        numpy.ndarray: Depths of sigma levels (3D matrix).
    """

    # Reading grid data
    with Dataset(gname, mode="r") as nc:
        h = nc.variables["h"][:]

    # Reading variables
    with Dataset(fname, mode="r") as nc:
        zeta = np.squeeze(nc.variables["zeta"][tindex, :, :])
        hmorph = nc.variables.get("hmorph", None)
        if hmorph is not None:
            h = np.squeeze(hmorph[tindex, :, :])

        # Reading sigma grid parameter
        try:
            theta_s = float(nc.variables["theta_s"][:])
            theta_b = float(nc.variables["theta_b"][:])
            Tcline = float(nc.variables["Tcline"][:])
            hc_var = nc.variables.get("hc", None)
            hc = float(hc_var[:]) if hc_var is not None else None
        except KeyError:
            theta_s = float(nc.theta_s)
            theta_b = float(nc.theta_b)
            Tcline = float(nc.Tcline)
            hc = float(nc.hc) if hasattr(nc, "hc") else None

        if hc is None:
            hc = min(float(np.min(h)), Tcline)

        N = len(nc.dimensions["s_rho"])
        VertCoordType = getattr(nc, "VertCoordType", "")
        vtrans = nc.variables.get("Vtransform", None)
        vtrans = int(np.squeeze(vtrans[:])) if vtrans is not None else None
        s_coord = 2 if VertCoordType == "NEW" or vtrans == 2 else 1

        if s_coord == 2:
            hc = Tcline

    if zeta is None:
        zeta = np.zeros_like(h)

    # Calculating depths of sigma levels
    vtype = point_type
    if point_type in ("u", "v"):
        vtype = "r"

    z = zlevs(h, zeta, theta_s, theta_b, hc, N, vtype, s_coord)

    # Swicth to U- or V- grid
    if point_type == "u":
        z = rho2u_3d(z)
    elif point_type == "v":
        z = rho2v_3d(z)

    return z


#################################################################
def get_hslice(fname, gname, vname, tindex, level, var_type):
    """
    Extracts a horizontal slice of a CROCO variable.

    Args:
        fname (str): Path to the CROCO file (average or history).
        gname (str): Path to the CROCO grid file.
        vname (str): Name of the variable.
        tindex (int): Time index.
        level (float): Vertical level for the slice:
                       - level = 0: 2D horizontal variable (e.g., zeta)
                       - level > 0: slice at the sigma level (integer >= 1)
                       - level < 0: interpolation at z = level.
        var_type (str): Type of the variable ('r', 'u', 'v', or 'w').

    Returns:
        numpy.ndarray: Horizontal slice (2D matrix).
    """

    with Dataset(fname, mode="r") as nc:
        if level == 0:
            # Variable 2D
            var = np.squeeze(nc.variables[vname][tindex, :, :])

            # Correction for wetting/drying
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
            # Slice at a specific sigma level
            var = np.squeeze(nc.variables[vname][tindex, level - 1, :, :])

        else:
            # Horizontal slice at a given depth level
            # Get the depths of sigma levels
            z = get_depths(fname, gname, tindex, var_type)

            # Read the 3D variable and interpolate
            var_sigma = np.squeeze(nc.variables[vname][tindex, :, :, :])
            var = vinterp(var_sigma, z, level)

    return var


#################################################################
def get_var(hisfile, gridfile, vname, tindex, vlevel, coef, rempts):
    """
    Retrieves a variable from a CROCO NetCDF file.

    Args:
        hisfile (str): Name of the NetCDF history file.
        gridfile (str): Name of the NetCDF grid file.
        vname (str): Name of the variable to retrieve.
        tindex (int): Time index of the variable.
        vlevel (float): Vertical level.
        coef (float): Multiplicative coefficient.
        rempts (list[int]): Points to remove at the boundaries
                            [west, east, south, north], or None to skip.

    Returns:
        tuple: (lat, lon, mask, var) where:
            lat (np.ndarray): Latitudes.
            lon (np.ndarray): Longitudes.
            mask (np.ndarray): Mask (1 = sea, NaN = land).
            var (np.ndarray): Requested variable.
    """

    print(f"-- {vname} --")

    if not gridfile:
        gridfile = hisfile
    if vname in ["h", "pm", "pn"]:
        lat, lon, mask = read_latlonmask(gridfile, "r")
        with Dataset(gridfile, "r") as nc:
            var = nc.variables[vname][:]

        var = mask * var

    else:
        var_type, vlev = get_type(hisfile, vname, vlevel)
        lat, lon, mask = read_latlonmask(gridfile, var_type)
        var = coef * mask * get_hslice(hisfile, gridfile, vname, tindex, vlev, var_type)

    # Remove boundary points if requested
    if rempts is not None:
        lat = rempoints(lat, rempts)
        lon = rempoints(lon, rempts)
        mask = rempoints(mask, rempts)
        var = rempoints(var, rempts)

    return lat, lon, mask, var


#################################################################
def do_plot(lon, lat, var, colmin, colmax, ncol, hisfile, gridfile, tindex, vlevel):
    """
    Do the plot (requires cartopy).
    """
    import cartopy.crs as ccrs
    from matplotlib.colors import Normalize

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
def get_vertical_section(fname, vname, tindex, direction, ind_sec, gname):
    """
    Extracts a vertical section of a 3D variable.

    Args:
        fname (str): Name of the NetCDF history file.
        vname (str): Name of the variable (e.g., "temp").
        tindex (int): Time index.
        direction (str): Section direction ("x" for zonal or "y" for meridional).
        ind_sec (int): Index of the row or column for the section.
        gname (str): Name of the NetCDF grid file.

    Returns:
        dict: Contains depths,
              distance (number of grid points),
              and interpolated values.
    """

    with Dataset(fname, "r") as nc:
        topo = np.squeeze(nc.variables["h"][:, :])
        # Déterminer le type de grille pour obtenir les profondeurs
        if vname == "w":
            point_type = "w"
        elif vname == "u":
            point_type = "u"
            topo = rho2u_2d(topo)
        elif vname == "v":
            point_type = "v"
            topo = rho2v_2d(topo)
        else:
            point_type = "r"

        # Get depths
        z = get_depths(fname, gname, tindex, point_type)
        var_3d = np.squeeze(nc.variables[vname][tindex, :, :, :])

    if direction == "y":
        section_data = var_3d[:, :, ind_sec]
        depth_section = z[:, :, ind_sec]
        distance = np.arange(section_data.shape[1])  # Distance in grid points
        topo_section = -topo[:, ind_sec]
    elif direction == "x":
        section_data = var_3d[:, ind_sec, :]
        depth_section = z[:, ind_sec, :]
        distance = np.arange(section_data.shape[1])  # Distance in grid points
        topo_section = -topo[ind_sec, :]
    else:
        raise ValueError("Direction invalide. Choisissez 'x' ou 'y'.")

    return {
        "topo": topo_section,
        "depth": depth_section,
        "distance": distance,
        "variable": section_data,
    }


#################################################################
def plot_vertical_section(section_data, title="Vertical Section", unit="°C",
                          show=True):
    """
    Plots a vertical section.

    Args:
        section_data (dict): Contains depths, distance, and values.
        title (str): Title of the plot.
        unit (str): Unit of the variable.
        show (bool): If True, call plt.show().

    Returns:
        matplotlib.figure.Figure: The created figure.
    """

    depth = section_data["depth"]
    topo = section_data["topo"]
    distance = section_data["distance"]
    variable = section_data["variable"]

    fig = plt.figure(figsize=(12, 6))
    plt.pcolormesh(distance, depth, variable, cmap="viridis", shading="auto")
    plt.plot(
        distance,
        topo,
        color="black",
        linestyle="--",
        linewidth=1.5,
    )
    plt.colorbar(label=unit)
    plt.title(title)
    plt.xlabel("Distance (grid points)")
    plt.ylabel("Depth (m)")
    plt.grid(True)

    if show:
        plt.show()

    return fig


################################################################
def vorticity(ubar, vbar, pm, pn):
    """
    Compute the relative vorticity.

    Translated from Matlab vort.m function.

    Args:
        ubar (np.ndarray): Zonal velocity on u-grid (Mp, L) or rho-grid (Mp, Lp)
        vbar (np.ndarray): Meridional velocity on v-grid (M, Lp) or rho-grid (Mp, Lp)
        pm (np.ndarray): Curvilinear coordinate metric in XI (1/dx) on rho-grid
        pn (np.ndarray): Curvilinear coordinate metric in ETA (1/dy) on rho-grid

    Returns:
        np.ndarray: Relative vorticity on psi-grid (M, L)
    """
    Mp, Lp = pm.shape
    L = Lp - 1
    M = Mp - 1

    # Compute uom = u / dx on the appropriate grid
    if ubar.shape == (Mp, L):
        # ubar already on u-grid
        uom = 2 * ubar / (pm[:, :L] + pm[:, 1:Lp])
    elif ubar.shape == (Mp, Lp):
        # ubar on rho-grid: interpolate to u-grid first
        ubar_u = rho2u_2d(ubar)  # shape (Mp, L)
        uom = 2 * ubar_u / (pm[:, :L] + pm[:, 1:Lp])
    else:
        raise ValueError(f"Unexpected ubar shape {ubar.shape} for pm shape {pm.shape}")

    # Compute von = v / dy on the appropriate grid
    if vbar.shape == (M, Lp):
        # vbar already on v-grid
        von = 2 * vbar / (pn[:M, :] + pn[1:Mp, :])
    elif vbar.shape == (Mp, Lp):
        # vbar on rho-grid: interpolate to v-grid first
        vbar_v = rho2v_2d(vbar)  # shape (M, Lp)
        von = 2 * vbar_v / (pn[:M, :] + pn[1:Mp, :])
    else:
        raise ValueError(f"Unexpected vbar shape {vbar.shape} for pn shape {pn.shape}")

    # Compute metric products at psi points
    mn = pm * pn
    mn_p = (mn[:M, :L] + mn[:M, 1:Lp] + mn[1:Mp, 1:Lp] + mn[1:Mp, :L]) / 4

    # Compute vorticity: ∂v/∂x - ∂u/∂y
    xi = mn_p * (von[:, 1:Lp] - von[:, :L] - uom[1:Mp, :] + uom[:M, :])

    return xi


################################################################
def tridim(var_2d, N):
    """
    Create a 3D array by replicating a 2D array N times along the first dimension.
    """
    return np.tile(var_2d[np.newaxis, :, :], (N, 1, 1))
