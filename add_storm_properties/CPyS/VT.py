import numpy as np

def find_max_min(xarr, il):
    '''
    Find the median of the highest and lowest 5 values of the data
    to test if odd values cause noise in Hart values
    '''
    size = xarr.sizes['snapshot']
    array_max = np.ma.zeros((size))
    array_min = np.ma.zeros((size))
    for it in range(size):
        values = xarr[it,:,:].values.flatten()
        #print('find vals ',it, il, values)
        if values.min() == 0.:
            array_min[it] = 0.0
            array_max[it] = 0.0
        else:
            values_sort_ascend = np.sort(values[values > 0.])
            val_min = np.ma.median(values_sort_ascend[:5])
            val_max = np.ma.median(values_sort_ascend[-6:])
            #print('find vals ',it, il, values_sort_ascend[:5])
            #val_min = np.ma.min(values_sort_ascend[:5])
            #val_max = np.ma.max(values_sort_ascend[-6:])
            array_min[it] = val_min
            array_max[it] = val_max
    return array_max, array_min

def VT(geopt, levels, name = "snap_zg"):
    """
    Parameters
    ----------
    geopt (xr.DataArray) : The Geopotential snapshots DataArray.
        plev must be decreasing
    name (str) : Name of the geopotential snapshots variable.

    Returns
    -------
    VTL, VTU : The Hart Phase Space parameters for upper and lower thermal wind respectively.
    """
    #from sklearn.linear_model import LinearRegression
    from scipy.stats import linregress
    Zmax1 = {}; Zmin1 = {}
    ΔZ = {}
    for il, level in enumerate(levels):
        #Zmax1[level] = geopt[il].max(["az", "r"])
        #Zmin1[level] = geopt[il].min(["az", "r"])
        Zmax1[level], Zmin1[level]  = find_max_min(geopt[il], il)
        ΔZ[level] = Zmax1[level] - Zmin1[level]

    print('VT, levels ',levels)
    # for each level, calculate the maximum and minimum value
    #Z_max = geopt[name].max(["az", "r"]) # Maximum of Z at each level for each snapshot
    #Z_min = geopt[name].min(["az", "r"]) # Minimum of ...
    #ΔZ = Z_max - Z_min  # Fonction of snapshot & plev
    #ΔZ_bottom = ΔZ.sel(plev=slice(950e2, 600e2)) # Lower troposphere
    ΔZ_bottom = []; ΔZ_top = []
    for i in range(len(ΔZ[levels[0]])):
        #ΔZ_bottom.append([ΔZ[levels[0]].values[i], ΔZ[levels[1]].values[i]])
        #ΔZ_top.append([ΔZ[levels[1]].values[i], ΔZ[levels[2]].values[i]])
        ΔZ_bottom.append([ΔZ[levels[0]][i], ΔZ[levels[1]][i]])
        ΔZ_top.append([ΔZ[levels[1]][i], ΔZ[levels[2]][i]])

    #ΔZ_bottom = ΔZ[900] - ΔZ[600]
    #print('ΔZ_bottom ',ΔZ_bottom)
    #ΔZ_top = ΔZ.sel(plev=slice(600e2, 250e2))    # Upper tropo
    #ΔZ_top = ΔZ[600] - ΔZ[200]
    # X is a list of the log of the pressure levels
    #X = np.log(ΔZ_bottom.plev).values.reshape(-1, 1).flatten()
    X = np.log(np.array([levels[0], levels[1]])).reshape(-1, 1).flatten()

    #VTL = [LinearRegression().fit(X, y).coef_[0] if not np.isnan(y).any() else np.nan for y in ΔZ_bottom.values]
    VTL = [linregress(X, y).slope if not np.isnan(y).any() else np.nan for y in ΔZ_bottom]
    print('len(VTL) ',len(VTL), np.ma.amax(VTL), np.ma.amin(VTL))
    large = np.where(np.array(VTL) > 160.0)[0]
    print('No large pts ',len(large), large)
    for i in range(len(large)):
        print('VTL ',i, np.array(VTL)[large[i]], np.array(ΔZ_bottom)[large[i]], Zmax1[925][large[i]], Zmin1[925][large[i]], np.array(X))
        #print('VTL ',i, np.array(VTL)[large[i]+1], np.array(ΔZ_bottom)[large[i]+1], np.array(X))
        print('')
    # upper level
    #X = np.log(ΔZ_top.plev).values.reshape(-1, 1).flatten()
    X = np.log(np.array([levels[1], levels[2]])).reshape(-1, 1).flatten()
    #VTU = [LinearRegression().fit(X, y).coef_[0] for y in ΔZ_top.values]
    VTU = [linregress(X, y).slope for y in ΔZ_top]
    #print('VTU ',VTU)

    for i in range(len(VTL)):
        if (57.5 < VTL[i] < 57.9) and (4.2 < VTU[i] < 4.3):
            print('VTL ',i, VTL[i], VTU[i], Zmax1[925][i], Zmin1[925][i])

    return VTL, VTU
