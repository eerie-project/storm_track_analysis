from .theta import theta_multitrack
from .B import B_vector
from .VT import VT

import numpy as np
import xarray as xr

def compute_CPS_parameters(
    tracks, geopt, levels, geopt_name = "snap_zg", plev_name="level",
):
    """
    Computes the three (+ theta) Hart parameters for all the points in tracks.

    Parameters
    ----------
    tracks (pd.DataFrame): The set of TC points
    geopt (xr.DataSet): The geopotential snapshots associated with the tracks
        level coordinate must be named plev, in Pa.
    names (str) : Provide the name of the 3D (plev, r, az) geopt snapshots variables as a string.

    Returns
    -------
    tracks (pd.DataFrame): The set of TC points with four new columns corresponding to the parameters
    """

    #old_settings = np.seterr(divide='ignore', invalid='ignore')

    #geopt = geopt.rename({plev_name:"plev"})

    # 1/ B computation
    mask_value = 0.0
    ## Select 900 & 600 hPa levels
    fname = geopt.format(level='925')
    geo = xr.open_dataset(fname)
    z900 = geo['snap_zg_925']
    z900m = z900.where((z900 > 0.0))
    fname = geopt.format(level='600')
    geo1 = xr.open_dataset(fname)
    z600 = geo1['snap_zg_600']
    z600m = z600.where((z600 > 0.0))
    fname = geopt.format(level='250')
    geo2 = xr.open_dataset(fname)
    z200 = geo2['snap_zg_250']
    z200m = z200.where((z200 > 0.0))
    #print('z900 ',z900)
    #z900, z600 = geopt[geopt_name].sel(plev = 900e2, method = "nearest"), \
    #                   geopt[geopt_name].sel(plev = 600e2, method = "nearest"),
    #print("Level "+str(z900.plev.values)+" is taken for 900hPa"+"\n"+
    #      "Level "+str(z600.plev.values)+" is taken for 600hPa"+"\n")

    ## theta computation
    if "theta" not in tracks.columns :
        tracks = tracks.assign(theta=theta_multitrack(tracks))

    ## B computation
    tracks = tracks.assign(
        B=B_vector(tracks.theta.values, z900m, z600m, tracks.lat.values)
    )

    # 2/ VTL & VTU computation
    #geopt = geopt.sortby("plev", ascending = False)
    geopt_list = [z900, z600, z200]
    VTL, VTU = VT(geopt_list, levels, name = geopt_name)

    # Output
    tracks = tracks.assign(VTL=VTL, VTU=VTU)
    #np.seterr(**old_settings)
    print('tracks ',tracks)

    return tracks


if __name__ == "__main__":
    from CPyS import *
    import numpy as np
    import pandas as pd
    import xarray as xr

    # Test theta_multitrack
    tracks = pd.read_csv("/home/h06/hadom/workspace/tenten/storms/composite_cps_phasespace/da746_trackprofile_19980801_19980901_tc_psl_subset.txt", index_col=False)
    geopt = xr.open_dataset("/scratch/hadom/snaps/zg_da746.nc")
    print('geopt ',geopt['snap_zg'])

    df = compute_CPS_parameters(tracks, geopt)
