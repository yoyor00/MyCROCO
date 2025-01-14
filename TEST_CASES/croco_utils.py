# croco_utils.py

import numpy as np

# Add other utility functions here as needed.


def zlevs(h, zeta, theta_s, theta_b, hc, N, type, vtransform):
    """
    Compute the depth of rho or w points for CROCO.

    Parameters:
        h (np.ndarray): Bathymetry (depths of the ocean bottom).
        zeta (np.ndarray): Free surface elevation.
        theta_s (float): S-coordinate surface control parameter.
        theta_b (float): S-coordinate bottom control parameter.
        hc (float): Width (m) of the surface or bottom boundary layer where higher resolution is required.
        N (int): Number of vertical levels.
        type (str): 'r' for rho points, 'w' for w points.
        vtransform (int): Vertical transformation equation (1 or 2).

    Returns:
        np.ndarray: Depths of RHO- or W-points (3D matrix).
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


def get_csf(sc, theta_s, theta_b):
    """
    function CSF = get_csf(sc,theta_s,theta_b)
    Get CS-Curves for the new s-coordinate system
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


def rho2u_3d(rho_in):
    """
    Convert a 3D field at rho points to a field at u points
    Calls rho2u_2d
    """

    def levloop(rho_in):
        Nlevs, Mshp, Lshp = rho_in.shape
        rho_out = np.zeros((Nlevs, Mshp, Lshp - 1))
        for k in xrange(Nlevs):
            rho_out[k] = ROMS.rho2u_2d(rho_in[k])
        return rho_out

    assert rho_in.ndim == 3, "rho_in must be 3d"
    return levloop(rho_in)


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


def rho2v_3d(rho_in):
    """
    Convert a 3D field at rho points to a field at v points
    Calls rho2v_2d
    """

    def levloop(rho_in):
        Nlevs, Mshp, Lshp = rho_in.shape
        rho_out = np.zeros((Nlevs, Mshp - 1, Lshp))
        for k in xrange(Nlevs):
            rho_out[k] = ROMS.rho2v_2d(rho_in[k])
        return rho_out

    assert rho_in.ndim == 3, "rho_in must be 3d"
    return levloop(rho_in)


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


def u2rho_3d(u_in):
    """
    Convert a 3D field at u points to a field at rho points
    Calls u2rho_2d
    """

    def _levloop(u_in):
        Nlevs, Mshp, Lshp = u_in.shape
        u_out = np.zeros((Nlevs, Mshp, Lshp + 1))
        for Nlev in xrange(Nlevs):
            u_out[Nlev] = ROMS.u2rho_2d(u_in[Nlev])
        return u_out

    assert u_in.ndim == 3, "u_in must be 3d"
    return _levloop(u_in)


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


def v2rho_3d(v_in):
    """
    Convert a 3D field at v points to a field at rho points
    Calls v2rho_2d
    """

    def levloop(v_in):
        Nlevs, Mshp, Lshp = v_in.shape
        v_out = np.zeros((Nlevs, Mshp + 1, Lshp))
        for Nlev in xrange(Nlevs):
            v_out[Nlev] = ROMS.v2rho_2d(v_in[Nlev])
        return v_out

    assert v_in.ndim == 3, "v_in must be 3d"
    return levloop(v_in)


def rho2uvp(rfield):
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1

    vfield = 0.5 * (rfield[np.arange(0, M), :] + rfield[np.arange(1, Mp), :])
    ufield = 0.5 * (rfield[:, np.arange(0, L)] + rfield[:, np.arange(1, Lp)])
    pfield = 0.5 * (ufield[np.arange(0, M), :] + ufield[np.arange(1, Mp), :])

    return ufield, vfield, pfield
