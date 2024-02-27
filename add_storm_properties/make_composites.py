"""
Create composites based on tracks input and variables required.
Based on code from Stella Bourdin and Bourdin et al. (2022), GMD

Reads in (netcdf) format of TC track files, from either TempestExtremes or TRACK (note in netcdf these are essentially the same format), and converts to csv. Note that TRACK has separate NH, SH hemisphere files.

Uses zg at 3 levels, here likely 925, 600, 250 hPa (data already extracted and available in netcdf format).
Uses TempestExtremes NodeFileCompose to calculate composites of zg at these 3 levels, using the TC tracks as the centre of the composite. Calculated on radial grid with a resolution of 0.25 degrees.

The composites are written out with variable names snap_zg_{level}, one file per level with all the composites in order of the tracks in the track file

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
import numpy as np
from copy import copy

fname_out = '{}a.{}{}_{}.{}'

composite_resol = '0.25'
composite_radx = '8'


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
    global um_runid, um_suiteid, cylc_task_cycle_time, time_cycle, startdate, enddate, year, month, day, cycleperiod, dir_out_base, variables, CPS_LEVELS, TRACK_TYPE, TEMPESTEXTREMES_ALGOS, FILEPATTERN_TE, FILEPATTERN_TRACK, resol, tempest_path

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

    year = time_cycle[0:4]
    month = time_cycle[4:6]
    day = time_cycle[6:8]

def write_to_csv(fname_csv, cols_to_keep, variables, index = False, sep=','):
    with open(fname_csv, 'w') as csvfile:
        writer = csv.writer(csvfile)
        csvfile.write(','.join(cols_to_keep)+'\n')
        
        #cols_to_keep = ['index', 'track_id', 'time', 'lon', 'lat', 'year', 'month', 'day', 'hour']
        #for row in zip(variables['index'], variables['track_id'], variables['time'], variables['lon'], variables['lat'], variables['year'], variables['month'], variables['day'], variables['hour'], variables['datetime']):
        for row in range(len(variables['index'])):
            column = []
            for c in cols_to_keep:
                column.append(variables[c][row])
            writer.writerow(column)
       
def convert_track_nc_to_csv(fname_nc, fname_csv):
    """
    Read track netcdf file
    FIRST_PT, NUM_PTS, TRACK_ID are dimensioned by no. tracks
    Rest are dimensioned by record
    Calculate a track_id for each track
    Convert this into a data frame csv file
    """
    track_dims = ['FIRST_PT', 'NUM_PTS', 'TRACK_ID']
    
    with Dataset(fname_nc, 'r') as nc:
        variables = {}
        variable_names = nc.variables.keys()
        print('variables ',variable_names)
        # number of storms in the file
        ntracks = int(nc.dimensions['tracks'].size)
        num_pts = nc.variables['NUM_PTS'][:]
        time_var = nc.variables['time']
        try:
            calendar = time_var.calendar
        except:
            calendar = time_var.time_calendar
        dtime = cftime.num2date(time_var[:], time_var.units, calendar = calendar)
        variables['year'] = []
        variables['month'] = []
        variables['day'] = []
        variables['hour'] = []
        variables['datetime'] = []
        for it, tim in enumerate(dtime):
            time_str = str(tim)
            daytime = time_str.split(' ')
            date = daytime[0]
            year = int(date[0:4])
            month = int(date[5:7])
            day = int(date[8:10])
            hour = int(daytime[1][:2])
            datet = datetime.datetime(year=year, month=month, day=day, hour=hour)
            variables['year'].append(year)
            variables['month'].append(month)
            variables['day'].append(day)
            variables['hour'].append(hour)
            variables['datetime'].append(datet)

        cols = list(variable_names)

        for var in cols:
            variables[var] = nc.variables[var][:]

        variables['track_id'] = []
        for it in range(ntracks):
            pts = num_pts[it]
            for ip in range(pts):
                variables['track_id'].append(it)

        cols_to_keep = ['index', 'track_id', 'time', 'lon', 'lat', 'year', 'month', 'day', 'hour', 'datetime']

        in_vars = ','.join(cols_to_keep)
        print('in_vars ',in_vars)

        list_nc = []
        for c in cols_to_keep:
            list_nc.append(list(variables[c][:]))
        df_nc = pd.DataFrame(list_nc)
        df_nc = df_nc.T
        df_nc.columns = cols_to_keep
        write_to_csv(fname_csv, cols_to_keep, variables, index = False, sep=',')

def tempest_composite(tempest_path, track_csv, in_data, out_data, var_name):
    '''
    Run the tempestExtremes NodeFileCompose command to create composites
    '''

    cmd = tempest_path+'/NodeFileCompose --in_nodefile "'+track_csv+'" --in_nodefile_type "SN" --in_fmt "(auto)" --in_data "'+in_data+'" --out_data "'+out_data+'" --out_grid "RAD" --dx '+composite_resol+' --resx '+composite_radx+' --var "'+var_name+'" --varout "'+var_name+'" --snapshots --regional --latname latitude --lonname longitude'
    print(cmd)
    run_cmd(cmd)

def form_composite(files, levels, var, fout, track_type):
    if len(files) > 0:

        for f in files:
            fname_nc = f
            fname_csv = fname_nc[:-3]+'.csv'
            fname_copy = fname_nc[:-3]+'_copy.nc'
            if not os.path.exists(fname_csv):
                convert_track_nc_to_csv(fname_nc, fname_csv)

            for lev in levels[var]:
                var_name = var+'_'+str(lev)
                in_data_file = fout[:-3]+'_'+str(lev)+'.nc'
                out_file = in_data_file[:-3]+'_snaps_'+track_type+'.nc'
                tempest_composite(tempest_path, fname_csv, in_data_file, out_file, var_name)


def main():

    get_environment_variables()
    output_dir = os.path.join(dir_out_base, um_suiteid)
    var = 'zg'

    track_types = TRACK_TYPE.split(',')
    cps_levels = CPS_LEVELS.split(',')
    levels = {}
    levels[var] = cps_levels
    tempestextremes_algos = TEMPESTEXTREMES_ALGOS.split(',')

    # filename string for the input data to make composite snapshots
    fout = os.path.join(output_dir, fname_out.format(um_suiteid[2:], 'pt', str(year), var, 'nc'))

    for track_type in track_types:
        if track_type == 'TE':
            for te_algo in tempestextremes_algos:
                search = FILEPATTERN_TE.format(runid=um_runid, year=str(year), algo_type=te_algo) 
                #search = '*_trackprofile_year_'+str(year)+'_tc_psl.nc'
                files = glob.glob(os.path.join(output_dir, search))
                form_composite(files, levels, var, fout, track_type+'_'+te_algo)
        else:
            for hemi in ['NH', 'SH']:
                search = FILEPATTERN_TRACK.format(year=str(year), resol=resol[0].upper()+resol[1:], runid=um_runid, hemi=hemi)
                #search = 'TC_'+str(year)+'*.full_track.NH.tracks.nc'
                files = glob.glob(os.path.join(output_dir, search))
                form_composite(files, levels, var, fout, track_type+'_'+hemi)

if __name__ == "__main__":

    main()

