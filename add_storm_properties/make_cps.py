"""
Calculates Hart phase space parameters and STJ, adds these parameters to existing tracks, and archives a new version of the track file to MASS

Read the wind speed and abs(u) running meaned files, locate the STJ (sub-tropical jet).
Add this information into the track files (csv format at this point), lat_STJ_{NH/SH} is the latitude, and ET 0/1 is an indicator of whether the storm is extra-tropical (0=no).

Use the zg composite snapshots from the make_composites.py step, and use these to calculate the Hart phase space parameters VTL, VTU and B.
Add these to the csv track file

Use the original netcdf track file to make a copy, and then append these extra variables into that new netcdf file.

Archive this file to MASS is possible.

"""

from netCDF4 import Dataset
import os, glob, sys
import subprocess
import time
import cftime
import datetime
import pandas as pd
import shutil
import csv
import xarray as xr
import numpy as np
from copy import copy

from CPyS import compute_CPS_parameters
from CPyS import STJ

# Extra variables to add to track file. 
extra_variables = ['lat_STJ_NH', 'lat_STJ_SH', 'ET', 'theta', 'B', 'VTL', 'VTU']
fname_out = '{}a.{}{}_{}.{}'

def run_cmd(
        cmd,
    check=True
):
    sts = subprocess.run(
        cmd,
        shell=True,
        universal_newlines=True,
        check=check,
        capture_output=True
    )
    print('cmd ', sts.stdout, sts.returncode)
    if sts.stderr:
        if 'Warning' in sts.stderr:            
            msg = (
                "Warning found in output ", sts.stderr
            )
            print(msg)
        else:
            msg = (
                "Error found in output ", sts.stderr
            )
            if 'TSSC_FILE_DOES_NOT_EXIST' in sts.stderr:
                raise Exception('File does not exist')
            else:
                raise RuntimeError(msg)
    return sts

def get_environment_variables():
    """
    Get the required environment variables from the suite. A list and
    explanation of the required environment variables is included in the
    documentation.
    """
    global um_suiteid, um_runid, cylc_task_cycle_time, time_cycle, startdate, enddate, year, month, day, cycleperiod, dir_out_base, variables, CPS_LEVELS, TRACK_TYPE, TEMPESTEXTREMES_ALGOS, FILEPATTERN_TE, FILEPATTERN_TRACK, resol, tempest_path, moo_base_path, backup_dir

    try:
        um_suiteid = os.environ["SUITEID_OVERRIDE"]
    except:
        um_suiteid = os.environ["CYLC_SUITE_NAME"]

    um_runid = um_suiteid.split('-')[1]
    cylc_task_cycle_time = os.environ["CYLC_TASK_CYCLE_TIME"]
    time_cycle = os.environ["TIME_CYCLE"]
    startdate = os.environ["STARTDATE"]
    enddate = os.environ["ENDDATE"]
    cycleperiod = os.environ["CYCLEPERIOD"]
    dir_out_base = os.environ["DIR_OUT"]
    calendar = os.environ["CALENDAR"]
    variables = os.environ["VARIABLES"]
    TRACK_TYPE = os.environ["TRACK_TYPE"]
    CPS_LEVELS = os.environ["CPS_LEVELS"]
    TEMPESTEXTREMES_ALGOS = os.environ["TEMPESTEXTREMES_ALGOS"]
    FILEPATTERN_TE = os.environ["FILEPATTERN_TE"]
    FILEPATTERN_TRACK = os.environ["FILEPATTERN_TRACK"]
    resol = os.environ["RESOL_ATM"]
    tempest_path = os.environ["TEMPESTEXTREMES_PATH"]
    moo_base_path = os.environ["MOO_BASE_PATH"]
    backup_dir = os.environ["BACKUP_DIR"]
    if not os.path.exists(backup_dir):
        os.makedirs(backup_dir)

    year = time_cycle[0:4]
    month = time_cycle[4:6]
    day = time_cycle[6:8]

def define_netcdf_metadata(var):
    """
    Define potential metadata for the netcdf variables
    """
    long_name = "unknown"
    description = "unknown"
    units = "1"

    if var == "VTU":
        standard_name = var
        long_name = "Upper level thermal wind indicator"
        description = "VTU component from Hart phase space"
        units = "1"
    elif var == "VTL":
        standard_name = var
        long_name = "Lower level thermal wind indicator"
        description = "VTL component from Hart phase space"
        units = "1"
    elif "lat" in var:
        standard_name = var
        long_name = "Latitude of sub-tropical jet"
        description = "Bourdin et al. (2022) treatment"
        units = "degrees"
    elif var == "ET":
        standard_name = var
        long_name = "Extra-tropical transition indication"
        description = "Indicator (0/1) whether this point is beyond STJ"
        units = "1"
    elif var == "B":
        standard_name = var
        long_name = "Extra-tropical transition indication"
        description = "Hart (2003) asymmetry parameter"
        units = "m"
    elif var == "theta":
        standard_name = var
        long_name = "Angular direction between track points"
        description = "Angular direction between track points"
        units = "m"

    return standard_name, long_name, description, units

def rewrite_nc_file(ncfile, track_csv_df, variables):
    '''
    ncfile - string - netcdf file to add new variables to
    track_csv_df - dataframe - containing the new variables to write (not necessarily in correct order)
    variables - list - the new variables to be added to the nc file
    '''
    track_dims = ['FIRST_PT', 'NUM_PTS', 'TRACK_ID']
    
    index_new = track_csv_df.index.values
    track_id_new = track_csv_df.track_id.values
    lon_values_new = track_csv_df.lon.values
        
    with Dataset(ncfile, 'r+') as nc:
        variables = {}

        record = nc.dimensions['record'].name
        tracks = nc.dimensions['tracks'].name
        record_len = nc.dimensions['record'].size
        tracks_len = nc.dimensions['tracks'].size

        lon = nc.variables['lon'][:]
        track_id = nc.variables['TRACK_ID'][:]
        first_pt = nc.variables['FIRST_PT'][:]
        num_pts = nc.variables['NUM_PTS'][:]

        variables_to_write = {}
        for var in extra_variables:
            nc.createVariable(var,'f8', ('record'))

            standard_name, long_name, description, v_units = define_netcdf_metadata(var)
            nc.variables[var].standard_name = standard_name
            nc.variables[var].long_name = long_name
            nc.variables[var].description = description
            nc.variables[var].units = str(v_units)

            variables_to_write[var] = []

            #VTU = track_csv_df.VTU.values
            var_new = track_csv_df[var].values
            print('var ',var, var_new)

            # loop round storms
            # loop round all points in storm
            # check that the extra values being written are from the same storm (i.e. same order) as in this file
            track_index = 0
            threshold = 0.01

            for tid in range(tracks_len):
                this_track = np.where(track_id_new == tid)[0]
                var_values = var_new[this_track]
                var_lon_vals = lon_values_new[this_track]

                track_firstpt = first_pt[tid]
                track_len = num_pts[tid]
                this_track_old = np.arange(track_firstpt, track_firstpt+track_len)
                lon_vals = lon[this_track_old]
                lon_vals_new = np.where(lon_vals > 180, lon_vals - 360, lon_vals)
                #lon_vals_new = lon_vals
                for i in range(len(this_track)):
                    if np.abs(lon_vals_new[i] - var_lon_vals[i]) > threshold:
                        print('Not same track ',tid, i, lon_vals_new[i], var_lon_vals[i])
                        raise Exception('Likely not same track ')
                    variables_to_write[var].append(var_values[i])
                track_index += track_len

            nc.variables[var][:] = variables_to_write[var]

def form_cps(files, fout, var, track_type, levels, output_dir, NH_lim, SH_lim):

    if len(files) == 0:
        raise Exception('No files available')

    for f in files:
        f_copy = f[:-3]+'_STJ_Hart.nc'
        # make copy of nc file 
        if os.path.exists(f_copy):
            os.remove(f_copy)
        shutil.copy(f, f_copy)

        f_csv = f[:-3]+'.csv'
        # now use snapshots and tracks to calculate Hart parameters
        # geopot is currently the ncfile name without the _level.nc
        tracks = pd.read_csv(f_csv, index_col=False, sep=',')
        #print('tracks ',tracks)

        tracks_STJ = STJ.identify_ET(tracks, NH_lim, SH_lim, fill=True)
        track_csv_STJ = f_csv[:-4]+'_STJ.csv'
        tracks_STJ.to_csv(track_csv_STJ, index=False)

        tracks = pd.read_csv(track_csv_STJ, index_col=False, sep=',')

        geopt = fout[:-3]+'_{level}_snaps_'+track_type+'.nc'
        df = compute_CPS_parameters(tracks, geopt, levels['zg'], plev_name='pressure')
        track_csv_hart = f_csv[:-4]+'_STJ_Hart.csv'
        df.to_csv(track_csv_hart, index=False)

        # Now convert tracks back to netcdf format
        rewrite_nc_file(f_copy, df, extra_variables)

        cmd = 'touch '+track_csv_hart+'.arch'
        run_cmd(cmd)
        cmd = 'touch '+f_copy+'.arch'
        run_cmd(cmd)

def archive_to_mass(output_dir):
    '''
    Archive any files that have a .arch extension
    '''
    moo_path_nc = moo_base_path.format(suite=um_suiteid)+'/any.nc.file/'
    moo_path_txt = moo_base_path.format(suite=um_suiteid)+'/ady.file/'
    files_to_archive = glob.glob(os.path.join(output_dir, '*.arch'))
    archived = []
    for f in files_to_archive:
        fname = f[:-5]
        fname_base = os.path.basename(fname)
        nc_file = False
        if fname[-2:] == 'nc':
            nc_file=True

        try:
            if nc_file:
                cmd = 'moo put -F '+fname+' '+moo_path_nc+fname_base
            else:
                cmd = 'moo put -F '+fname+' '+moo_path_txt+fname_base
            run_cmd(cmd)
        except:
            cmd = 'cp '+fname+' '+os.path.join(backup_dir, um_suiteid, fname_base)
            run_cmd(cmd)
        archived.append(f)

    for f in archived:
        os.remove(f)

def main():
    get_environment_variables()

    output_dir = os.path.join(dir_out_base, um_suiteid)
    var = 'zg'

    track_types = TRACK_TYPE.split(',')
    cps_levels = CPS_LEVELS.split(',')
    levels = {}
    levels[var] = list(map(int, cps_levels))
    tempestextremes_algos = TEMPESTEXTREMES_ALGOS.split(',')

    # filename string for the input data to make composite snapshots
    fout = os.path.join(output_dir, fname_out.format(um_suiteid[2:], 'pt', str(year), var, 'nc'))

    # STJ independent of tracks
    outfile = os.path.join(output_dir, fname_out.format(um_suiteid[2:], 'pt', str(year), 'STJ', 'nc'))
    fout_u_abs_runmean = os.path.join(output_dir, fname_out.format(um_suiteid[2:], 'pt', str(year), 'u_250_abs_runmean', 'nc'))
    fout_wspd_runmean = os.path.join(output_dir, fname_out.format(um_suiteid[2:], 'pt', str(year), 'wspd_250_runmean', 'nc'))

    STJ_files = glob.glob(outfile+'*.nc')
    if len(STJ_files) != 2:
        NH_lim, SH_lim = STJ.compute_STJ_latmin(fout_wspd_runmean, fout_u_abs_runmean, outfile)
    else:
        NH_lim = xr.open_dataset(outfile+'_NH.nc').latitude
        SH_lim = xr.open_dataset(outfile+'_SH.nc').latitude

    for track_type in track_types:
        if track_type == 'TE':
            for te_algo in tempestextremes_algos:
                search = FILEPATTERN_TE.format(runid=um_runid, year=str(year), algo_type=te_algo) 
                files = glob.glob(os.path.join(output_dir, search))
                track_type_val = track_type+'_'+te_algo
                form_cps(files, fout, var, track_type_val, levels, output_dir, NH_lim, SH_lim)
        else:
            for hemi in ['NH', 'SH']:
                search = FILEPATTERN_TRACK.format(year=str(year), resol=resol[0].upper()+resol[1:], runid=um_runid, hemi=hemi)
                files = glob.glob(os.path.join(output_dir, search))
                track_type_val = track_type+'_'+hemi
                form_cps(files, fout, var, track_type_val, levels, output_dir, NH_lim, SH_lim)

    archive_to_mass(output_dir)

if __name__ == "__main__":
    main()
