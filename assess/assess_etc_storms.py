""" 
Code to read TC storms into dictionary structure
Call code to calculate TC properties (frequency, density etc) by year and means
Call code to produce plots based on the TC properties
"""
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os, sys, glob
import subprocess
import cf_units
import calendar

TC_SCRIPT_PATH='/data/users/hadom/branches/git/storm_track_analysis/assess'

sys.path.append(TC_SCRIPT_PATH)
import load_data.load_feature_file
import load_data.load_TE_nc as load_TE_nc
import load_data.get_tracks_from_mass
import etc_calculate
import etc_assessment
import tc_intensity
import tc_interannual
import tc_assess_code
import tc_assess_code.calc_interann
import ts_obs.load_data
import ts_obs.track_ibtracs_netcdf as track_ibtracs_netcdf
import make_simple_webpage

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER,LONGITUDE_FORMATTER
import matplotlib.ticker as mticker 
import shapely.geometry as sgeom
import iris
import iris.plot as iplt

import netCDF4, cftime
import numpy as np
import pickle
import time

# some example inputs for the plotting
months = [1,2,3,4,5,6,7,8,9,10,11,12]
months_nh = [11,12,1,2,3]
months_sh = [5,6,7,8,9]

YEARS_OBS = range(1979, 2019)

def calc_data_freq(storms):
    '''
    Calculate the data frequency - 6hr, 1hr
    '''
    hours = []
    medians = []
    for storm in storms[:10]:
        if len(storm.obs) > 6:
            for ob in storm.obs:
                hours.append(ob.date.hour)
        time_diff = np.asarray(hours[1:]) - np.asarray(hours[:-1])
        median = np.median(np.asarray(time_diff))
        medians.append(median)
    median_all = np.median(np.array(medians))
    if median != 1 and median != 6:
        raise Exception('Data frequency not 1 or 6 '+str(median))
    #print('Data freq ',median)
    return median

def work(runid_info, yearstart, yearend, yearsplot, dir_in_base, plot_dir = '', apply_threshold = False, coupled=True, storedir='./', do_meanstate=True, do_variability=True):
    '''
    Read storms from dataset
    Calculate and plot various metrics/assessment aspects
    '''
    time_start = time.perf_counter()

    models = '_'.join(runid_info['model'])
    expt = runid_info['experiment'][0]
    if expt != '':
        expt_end = '_'+expt
    else:
        expt_end = ''

    if plot_dir == '':
        plot_dir = os.path.join('./', models+expt_end)
    else:
        plot_dir = os.path.join(plot_dir, models+expt_end)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    years = np.arange(yearstart, yearend+1)
    if years[0] < YEARS_OBS[0]:
        years_obs_start = YEARS_OBS[0]
    else:
        years_obs_start = years[0]
    if years[-1] > YEARS_OBS[-1]:
        years_obs_end = YEARS_OBS[-1]
    else:
        years_obs_end = years[-1]
    years_obs = range(years_obs_start, years_obs_end)
    years_obs = range(YEARS_OBS[0], YEARS_OBS[-1])
    print('years_obs ',years_obs)

    storms = {}; no_years_with_storms = {}
    keys = []
    for ir, suiteid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            runid = suiteid[2:]
            expt = runid_info['experiment'][ir]
            ens = runid_info['ens'][ir]
            fileformat = runid_info['fileformat'][ir]
            filepattern = runid_info['filepattern'][ir]
            variable_indices = runid_info['variable_indices'][ir]
            calendar = runid_info['calendar'][ir]
            key = etc_calculate.define_key(algo, expt, ens)
            if 'track' in algo or 'TRACK' in algo:
                try:
                    track_method = algo.split('_')[1]
                except:
                    track_method = ''
            else:
                track_method = ''
            resol_in = runid_info['resol'][ir]
            algo_type = runid_info['algo_type'][ir]
            institute = runid_info['institute'][ir]
            rip = runid_info['rip'][ir]
            grid = runid_info['grid'][ir]
            get_from_mass = runid_info['get_from_mass'][ir]
            dir_in = os.path.join(dir_in_base, suiteid, 'tracks', ens)

            print('get mass ',resol_in, algo_type, get_from_mass, dir_in, suiteid, runid)

            if get_from_mass:
                # get data from archive
                load_data.get_tracks_from_mass.work(suiteid, runid, ens, runid_info, algo, algo_type, expt, key, resol_in, track_method, dir_in, years, institute, rip, grid, fileformat='nc', filepattern=filepattern)

            # load storm data from local store into dictionary structure
            storms_runid, no_years_with_storms_runid, year_start, year_end = load_data.load_feature_file.read_storms_with_appropriate_reader(suiteid, runid, ens, runid_info, algo, algo_type, expt, key, resol_in, track_method, dir_in, years, institute, rip, grid, fileformat=fileformat, filepattern=filepattern, variable_indices=variable_indices, calendar=calendar)

            storms[(suiteid, key)] = storms_runid
            storms[(suiteid, key, 'nyears')] = no_years_with_storms_runid
            storms[(suiteid, key, 'years')] = np.arange(year_start, year_end+1)

            data_freq = calc_data_freq(storms_runid)
            storms[(suiteid, key, 'freq')] = data_freq
            keys.append(key)
        
    time_read_data = time.perf_counter()
    print('Reading data took ',(time_read_data - time_start),' seconds')

    # calculate all the TC statistics and write to local disk
    etc_calculate.calculate_etc_properties(runid_info, storms, years, storedir, plot_dir, keys, months_nh, months_sh)

    # plot TC properties
    etc_assessment.plot_etc_properties(runid_info, storms, years, storedir, do_meanstate, do_variability, years_obs, plot_dir, keys, months_nh, months_sh, coupled=coupled)

    # make some web pages from the plots produced above
    algo = runid_info['algorithm']
    algo_fname = runid_info['algo_type'][0]
    make_simple_webpage.make_web_pages(plot_dir, os.path.join(plot_dir, 'html'), runid_info['model'], algo, algo_fname)

if __name__ == '__main__':
    pass
