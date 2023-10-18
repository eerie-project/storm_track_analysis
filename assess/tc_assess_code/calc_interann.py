'''
Calculate the TC frequency and ACE for all basins for a given year, including by storm category
Write out as netcdf file with year, and basins written into metadata
Try and make generic enough to run elsewhere, e.g. JASMIN, so pull out directory assumptions to highest level
Try to make self-contained, so not needing to pull code from other packages

Inputs:
suitename
(directory containing track data)
storms for this year
basins to calculate
where to write netcdf output file
'''
import sys
from netCDF4 import Dataset
from cftime import date2num, datetime
import numpy as np
import os

TC_SCRIPT_PATH='/data/users/hadom/branches/git/track_analysis'
sys.path.append(TC_SCRIPT_PATH)
import tc_assess_code
#sys.path.append('/home/h06/hadom/workspace/tenten/storms/assessment')
#import tc_assess_code

def calculate_model_count_ace(runid, storms, year, basins, months, hemi, nensemble=1, ace=False):
    '''
    Calculate the model tc frequency in each basin
    : param str runid:
    : param collection storms:
    : param int year:
    : param list basins:
    : param list months:
    : param int nensemble:
    : returns: dictionary of storm frequency, with basins as keys
    '''
    no_storms = False
    annual_count = {}
    for ib, basin in enumerate(basins):
        count = None
        count_category = None
        years = np.arange(year, year+1)

        #print('process basin ',ib, basin, year)
        if not ace:
            count = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_count(storms, years, months, basin, nensemble), np.float32)
        else:
            count = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_ace(storms, years, months, basin, nensemble), np.float32)
            #print ('basin, ace ',runid, basin, count)
        annual_count.update({basin: count})

        if not ace:
            count_category = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_count_category(storms, years, months, basin, nensemble), np.float32)
        else:
            count_category = np.ma.asarray(tc_assess_code.get_annual_vmax_mean_ace_category(storms, years, months, basin, nensemble), np.float32)
        for icat in range(6):
            annual_count.update({(basin,icat): count_category[:,icat]})
            #print('annual_count cat ',basin, icat, annual_count[(basin,icat)])
        
    # when the NH (or SH) total is zero, then this is missing, so set the appropriate hemisphere BASINS to missing
    miss = np.where(annual_count[hemi] == 0)[0]
    if len(miss) >= 1:
        no_storms = True
    annual_count[hemi][miss] = np.ma.masked
    
    for icat in range(6):
        annual_count[(hemi,icat)][miss] = np.ma.masked
    for ib, basin in enumerate(basins):
        annual_count[basin][miss] = np.ma.masked
        for icat in range(6):
            annual_count[(basin,icat)][miss] = np.ma.masked

    return annual_count, no_storms

def read_model_count_ace(fname_in, years, runkey, annual_count_models, outdir):
    '''
    Read the netcdf files for the annual count, ace values and store
    : param str runid:
    : param collection storms:
    : param int year:
    
    '''

    for hemi in ['nh', 'sh']:
        years_hem = []
        for year in years:
            fname = fname_in.format(str(year), hemi)
            if os.path.exists(fname):
                years_hem.append(year)
                with Dataset(fname, 'r') as nc:
                    basins = nc.basins
                    count = nc.variables['tc_count'][:]
                    ace = nc.variables['tc_ace'][:]
                    count_by_cat = nc.variables['tc_count_by_cat'][:,:]
                    ace_by_cat = nc.variables['tc_ace_by_cat'][:,:]
                    for ib, bas in enumerate(basins):
                        annual_count_models[runkey][bas].append(count[ib])
                        annual_count_models[runkey]['ace', bas].append(ace[ib])
                        for icat in range(6):
                            annual_count_models[runkey][(bas, icat)].append(count_by_cat[ib, icat])
                            annual_count_models[runkey][('ace', bas, icat)].append(ace_by_cat[ib, icat])
            annual_count_models[runkey][('years', hemi)] = years_hem
            

def write_tc_count_ace_netcdf(fname_out, basins, annual_count, annual_ace, 
                              year, time_units, calendar, 
                              start_date, end_date, algorithm, 
                              runid, hemi='nh'):
    '''
    Write out netcdf file containing the TC count and ACE for each basin
    for this one year
    '''
    nbasins = len(basins)
    print('write out ',fname_out)
    with Dataset(fname_out, 'w', format='NETCDF4') as nc:

        nc.title = "TC count and ACE by basin"
        nc.start_date = start_date
        nc.end_date = end_date
        nc.algorithm = algorithm
        nc.runid = runid
        nc.basins = basins
    
        nc.createDimension('basins', size=nbasins)
        nc.createDimension('category', size=6)
        nc.createDimension('time', size=1)
        nc.createVariable("time", "f8", ("time"))
        nc.createVariable('tc_count', np.int32, ('basins'))
        nc.createVariable('tc_ace', np.float32, ('basins'))
        nc.createVariable('tc_count_by_cat', np.int32, ('basins','category'))
        nc.createVariable('tc_ace_by_cat', np.float32, ('basins','category'))

        nc.variables["time"].units = time_units
        nc.variables["time"].calendar = calendar
        nc.variables["time"].standard_name = "time"
        nc.variables["time"].long_name = "time"

        nc.variables['tc_count'].units = '1'
        nc.variables['tc_count'].standard_name = 'tc_count'
        nc.variables['tc_count'].long_name = 'TC count per year'
        nc.variables['tc_count'].description = 'TC count per year'

        nc.variables['tc_ace'].units = '1'
        nc.variables['tc_ace'].standard_name = 'tc_ace'
        nc.variables['tc_ace'].long_name = 'TC ace per year'
        nc.variables['tc_ace'].description = 'TC ace (from 925hPza winds) per year'

        nc.variables['tc_count_by_cat'].units = '1'
        nc.variables['tc_count_by_cat'].standard_name = 'tc_count_by_cat'
        nc.variables['tc_count_by_cat'].long_name = 'TC count by category per year'
        nc.variables['tc_count_by_cat'].description = 'TC count by category (MSLP intensity) per year'

        nc.variables['tc_ace_by_cat'].units = '1'
        nc.variables['tc_ace_by_cat'].standard_name = 'tc_ace_by_cat'
        nc.variables['tc_ace_by_cat'].long_name = 'TC ace by category per year'
        nc.variables['tc_ace_by_cat'].description = 'TC ace (from 925hPa winds) by category (intensity) per year'

        year = int(start_date[0:4])
        month = int(start_date[4:6])
        day = int(start_date[6:8])
        hour = 0

        t1 = date2num(
                datetime(
                    year,
                    month,
                    day,
                    hour,
                    calendar=calendar,
                ),
                time_units,
                calendar=calendar,
            )

        nc.variables['time'][:] = t1
        tc_count = []
        for basin in basins:
            tc_count.append(annual_count[basin])
        nc.variables['tc_count'][:] = tc_count
        tc_ace = []
        for basin in basins:
            tc_ace.append(annual_ace[basin])
        nc.variables['tc_ace'][:] = tc_ace

        for ib, basin in enumerate(basins):
            tc_count = []
            for icat in range(6):
                tc_count.append(annual_count[(basin, icat)][0])
            #print('tc_count for cat ',basin, tc_count)
            nc.variables['tc_count_by_cat'][ib,:] = tc_count
        for ib, basin in enumerate(basins):
            tc_ace = []
            for icat in range(6):
                tc_ace.append(annual_ace[(basin, icat)])
            nc.variables['tc_ace_by_cat'][ib, :] = tc_ace

def calc_tc_frequency_mean_std(runid, runkey, BASINS, annual_count_obs, annual_count_models, mean_count, std_count, do_obs=False):
    '''
    Calculate some statistics of the model interannual variability
    param:
    '''
    if do_obs:
        for ods in OBS_DATASETS:
            for ib, basin in enumerate(BASINS):
                for otype in OBS_TYPES:
                    mean_count[(ods, otype, basin)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin)].mean(), dtype = np.float32)
                    std_count[(ods, otype, basin)] =  np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype],basin)].std(), dtype = np.float32)
                    for icat in np.arange(6):        
                        mean_count_cat[(ods, otype, basin, icat) ] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin, icat)].mean(), dtype = np.float32)
                        std_count_cat[(ods, otype, basin, icat)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], basin, icat)].std(), dtype = np.float32)

        # do ACE
        for ods in OBS_DATASETS:
            for ib, basin in enumerate(BASINS):
                for otype in OBS_TYPES:
                    mean_count[(ods, otype, 'ace', basin)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin)].mean(), dtype = np.float32)
                    std_count[(ods, otype, 'ace', basin)] =  np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin)].std(), dtype = np.float32)
                    for icat in np.arange(6):        
                        mean_count_cat[(ods, otype, 'ace', basin, icat) ] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin, icat)].mean(), dtype = np.float32)
                        std_count_cat[(ods, otype, 'ace', basin, icat)] = np.ma.asarray(annual_count_obs[ods][(OBS_TYPES[otype], 'ace', basin, icat)].std(), dtype = np.float32)

    for ib, basin in enumerate(BASINS):
        mean_count[(runkey, basin)] = np.ma.asarray(np.array(annual_count_models[runkey][basin]).mean(), dtype = np.float32)
        std_count[(runkey, basin)] = np.ma.asarray(np.array(annual_count_models[runkey][basin]).std(), dtype = np.float32)
        #print('annual_count ', annual_count_models[runkey][(basin)], mean_count[(runkey, basin)], std_count[(runkey, basin)])
        for icat in np.arange(6):        
            mean_count[(runkey, basin, icat)] = np.ma.asarray(np.array(annual_count_models[runkey][(basin,icat)]).mean(), dtype = np.float32)
            std_count[(runkey, basin, icat)] = np.ma.asarray(np.array(annual_count_models[runkey][(basin,icat)]).std(), dtype = np.float32)

    # do ACE
    for ib, basin in enumerate(BASINS):
        #print('calc mean ace ',runid, runkey, basin)

        mean_count[(runkey, 'ace', basin)] = np.ma.asarray(np.array(annual_count_models[runkey][('ace', basin)]).mean(), dtype = np.float32)
        std_count[(runkey, 'ace', basin)] = np.ma.asarray(np.array(annual_count_models[runkey][('ace', basin)]).std(), dtype = np.float32)
        #print('annual_ace ', annual_count_models[runkey][('ace',basin)], mean_count[(runkey, 'ace', basin)], std_count[(runkey, 'ace', basin)])
        for icat in np.arange(6):        
            mean_count[(runkey, 'ace', basin, icat)] = np.ma.asarray(np.array(annual_count_models[runkey][('ace', basin, icat)]).mean(), dtype = np.float32)
            std_count[(runkey, 'ace', basin, icat)] = np.ma.asarray(np.array(annual_count_models[runkey][('ace', basin, icat)]).std(), dtype = np.float32)

