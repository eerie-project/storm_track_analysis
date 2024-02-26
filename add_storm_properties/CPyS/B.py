import pandas as pd
import numpy as np

def right_left(field, th):
    """
    Separate geopotential field into left and right of the th line.

    Parameters
    ----------
    field (xr.DataArray): The geopotential field
    th: The direction (in degrees)

    Returns
    -------
    left, right (2 xr.DataArray): The left and right side of the geopt. field.
    """
    if th <= 180:
        return field.where((field.az <= th) | (field.az > 180 + th)), field.where(
            (field.az > th) & (field.az <= 180 + th)
        )
    else:
        return field.where((field.az <= th) & (field.az > th - 180)), field.where(
            (field.az > th) | (field.az <= th - 180)
        )


def right_left_vector(z, th):
    """
    Separate geopotential field into left and right of the th line.

    Parameters
    ----------
    z (xr.DataArray): The geopotential field
    th: The direction (in degrees)

    Returns
    -------
    left, right (2 xr.DataArray): The left and right side of the z. field.
    """

    A = pd.DataFrame([list(z.az.values)] * len(z.snapshot))  # matrix of az x snapshot
    mask = np.array((A.lt(th, 0) & A.ge(th - 180, 0)) | A.ge(th + 180, 0))  # Mask in 2D (az, snapshot)
    mask = np.array([mask] * len(z.r))  # Mask in 3D (r, az, snapshot)
    mask = np.swapaxes(mask, 0, 1)  # Mask in 3D (az, r, snapshot)
    R, L = z.where(mask), z.where(
        ~mask
    )
    return R, L


def area_weights(field):
    """
    Computes the weights needed for the weighted mean of polar field.

    Parameters
    ----------
    field (xr.DataArray): The geopotential field

    Returns
    -------
    w (xr.DataArray): The weights corresponding to the area wrt the radius.
    """
    δ = (field.r[1] - field.r[0]) / 2
    w = (field.r + δ) ** 2 - (field.r - δ) ** 2
    return w


def B(th, geopt, SH=False, names=["snap_z900", "snap_z600"]):
    """
    Computes the B parameter for a point, with the corresponding snapshot of geopt at 600hPa and 900hPa

    Parameters
    ----------
    th: The direction (in degrees)
    geopt (xr.DataSet): The snapshots at both levels
    SH (bool): Set to True if the point is in the southern hemisphere
    names: names of the 900hPa and 600hPa geopt. variables in geopt.

    Returns
    -------
    B, the Hart phase space parameter for symetry.
    """
    if type(names) == str:
        z900 = geopt[names].sel(plev=900e2, method="nearest")
        print("Level " + str(z900.plev.values) + " is taken for 900hPa")
        z600 = geopt[names].sel(plev=600e2, method="nearest")
        print("Level " + str(z600.plev.values) + " is taken for 600hPa")
    else:
        z900 = geopt[names[0]]
        z600 = geopt[names[1]]

    ΔZ = z600 - z900
    ΔZ_R, ΔZ_L = right_left_vector(ΔZ, th)
    if SH:
        h = -1
    else:
        h = 1
    return h * (
            ΔZ_R.weighted(area_weights(ΔZ_R)).mean(["r", "az"])
            - ΔZ_L.weighted(area_weights(ΔZ_L)).mean(["r", "az"])
    ).values


def B_vector(th_vec, z900, z600, lat):
    """
    Computes the B parameter for a vector of points, with the corresponding snapshot of geopt at 600hPa and 900hPa

    Parameters
    ----------
    th_vec : The theta parameter for each point
    z900 : The z900 field for each point
    z600 : The z600 field for each point
    lat : The latitude of each point

    Returns
    -------
    B, the Hart phase space parameter for symetry.
    """
    ΔZ = z600 - z900
    ΔZ_R, ΔZ_L = right_left_vector(ΔZ, th_vec)
    h = np.where(lat < 0, -1, 1)
    return h * (
        ΔZ_R.weighted(area_weights(ΔZ_R)).mean(["az", "r"])
        - ΔZ_L.weighted(area_weights(ΔZ_L)).mean(["az", "r"])
    )
