# croco_utils.py

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.colors import Normalize
import cartopy.crs as ccrs


################################################################
def zlevs(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    Compute the depth of rho or w points for CROCO.
    """
    Ndim = np.ndim(h)
    if Ndim == 2:
        M, L = h.shape
    elif Ndim == 1:
        M = len(h)
        L = 1
        h = h[:, np.newaxis]
        zeta = zeta[:, np.newaxis]
    else:
        raise ValueError("zlevs: error - incorrect dimension for h")

    ds = 1.0 / float(N)

    if type == "w":
        sc = ds * (np.arange(0, N + 1) - N)
        Nmax = N + 1
    elif type == "r":
        sc = ds * (np.arange(1, N + 1) - N - 0.5)
        Nmax = N
    else:
        raise ValueError(f"Problem with type = {type}")

    if vtransform == 1:
        cff1 = 1.0 / np.sinh(theta_s)
        cff2 = 0.5 / np.tanh(0.5 * theta_s)
        Cs = (1.0 - theta_b) * cff1 * np.sinh(theta_s * sc) + theta_b * (
            cff2 * np.tanh(theta_s * (sc + 0.5)) - 0.5
        )
    elif vtransform == 2:
        Cs = get_csf(sc, theta_s, theta_b)
    else:
        raise ValueError(f"Problem with vtransform = {vtransform}")

    z = np.zeros([Nmax, M, L])

    for k in np.arange(0, Nmax):
        if vtransform == 1:
            cff = hc * (sc[k] - Cs[k])
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / h
            z[k, :, :] = z0 + zeta * (1.0 + z0 * hinv)
        elif vtransform == 2:
            cff = hc * sc[k]
            cff1 = Cs[k]
            z0 = cff + cff1 * h
            hinv = 1.0 / (h + hc)
            z[k, :, :] = z0 * h * hinv + zeta * (1.0 + z0 * hinv)

    return z.squeeze()


################################################################
def get_csf(sc, theta_s, theta_b):
    """
    Get CS-Curves for the new s-coordinate system
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


################################################################
def rho2u_2d(rho_in):
    """
    Convert a 2D field at rho points to a field at u points
    """

    def _r2u(rho_in, Lp):
        u_out = rho_in[:, : Lp - 1]
        u_out += rho_in[:, 1:Lp]
        u_out *= 0.5
        return u_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2u(rho_in, Lshp)


################################################################
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
        v_out = rho_in[: Mp - 1]
        v_out += rho_in[1:Mp]
        v_out *= 0.5
        return v_out.squeeze()

    assert rho_in.ndim == 2, "rho_in must be 2d"
    Mshp, Lshp = rho_in.shape
    return _r2v(rho_in, Mshp)


################################################################
def rho2v_3d(var):
    """
    Convert 3D variable from rho-grid to v-grid.
    """
    return 0.5 * (var[:, :-1, :] + var[:, 1:, :])


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


################################################################
def rho2uvp(rfield):
    """
    Compute the values at u, v and psi points from rho points.

    Translated from Matlab rho2uvp.m function.
    """
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1

    # Interpolate to v points (average in eta direction)
    vfield = 0.5 * (rfield[:M, :] + rfield[1:Mp, :])

    # Interpolate to u points (average in xi direction)
    ufield = 0.5 * (rfield[:, :L] + rfield[:, 1:Lp])

    # Interpolate to psi points (average u field in eta direction)
    pfield = 0.5 * (ufield[:M, :] + ufield[1:Mp, :])

    return ufield, vfield, pfield


################################################################
def vorticity(ubar, vbar, pm, pn):
    """
    Compute the relative vorticity.

    Translated from Matlab vort.m function.

    Args:
        ubar (np.ndarray): Zonal velocity on u-grid
        vbar (np.ndarray): Meridional velocity on v-grid
        pm (np.ndarray): Curvilinear coordinate metric in XI (1/dx) on rho-grid
        pn (np.ndarray): Curvilinear coordinate metric in ETA (1/dy) on rho-grid

    Returns:
        np.ndarray: Relative vorticity on psi-grid
    """
    Mp, Lp = pm.shape
    L = Lp - 1
    M = Mp - 1

    # Initialize arrays
    xi = np.zeros((M, L))
    mn_p = np.zeros((M, L))
    uom = np.zeros((M, Lp))
    von = np.zeros((Mp, L))

    # Check if ubar dimensions match what Matlab expects
    if ubar.shape == (Mp, L):  # ubar on u-grid as expected by Matlab
        uom = 2 * ubar / (pm[:, :L] + pm[:, 1:Lp])
    else:
        # Interpolate ubar to rho grid first, then extract appropriate slice
        ubar_rho = u2rho_2d(ubar)
        uom = 2 * ubar_rho[:M, :] / (pm[:M, :L] + pm[:M, 1:Lp])

    # Check if vbar dimensions match what Matlab expects
    if vbar.shape == (M, Lp):  # vbar on v-grid as expected by Matlab
        von = 2 * vbar / (pn[:M, :] + pn[1:Mp, :])
    else:
        # Interpolate vbar to rho grid first, then extract appropriate slice
        vbar_rho = v2rho_2d(vbar)
        von = 2 * vbar_rho[:, :L] / (pn[:M, :L] + pn[1:Mp, :L])

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


################################################################
def read_latlonmask(fname, var_type):
    """
    Read latitudes and longitudes and mask from CROCO netcdf file.
    """
    nc = Dataset(fname)
    lat = nc.variables["lat_rho"][:]
    lon = nc.variables["lon_rho"][:]
    mask = nc.variables.get("mask_rho", np.ones_like(lon))

    if var_type == "u":
        lat = rho2u_2d(lat)
        lon = rho2u_2d(lon)
        mask = mask[:, :-1] * mask[:, 1:]
    elif var_type == "v":
        lat = rho2v_2d(lat)
        lon = rho2v_2d(lon)
        mask = mask[:-1, :] * mask[1:, :]

    mask = np.where(mask == 0, np.nan, mask)
    return lat, lon, mask
