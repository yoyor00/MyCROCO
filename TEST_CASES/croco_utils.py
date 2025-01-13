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
