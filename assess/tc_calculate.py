import time
import os, sys
import iris

TC_SCRIPT_PATH='/data/users/hadom/branches/git/eerie-project/storm_track_analysis/assess'
sys.path.append(TC_SCRIPT_PATH)
import tc_assess_code
import tc_assess_code.calc_interann

BASINS_LIST = tc_assess_code.BASINS_LIST
BASINS = tc_assess_code.BASINS

def define_key(algo, expt, ens):
    key = algo
    if expt != '':
        key += '_'+expt
    if ens != '':
        key += ens

    return key

def storm_density_calculation(storms, years, months, basin, genesis=False, nyears=0, data_freq=6):
    """ 
    Plots monthly mean density of tropical storm tracks, genesis locations, 
    lysis locations or location of maximum intensity of model tropical storms
    that form in a desired set of years, months and ocean basin. 
    Default plot is of storm track density.  
    
    To get different plots set:
    Genesis plot: genesis=True
    Lysis plot: lysis=True
    Maximum intensity (location of max wind): 
    max_intensity=True
    
    Basin options:
    None: Whole globe
    'na': North Atlantic
    'ep': Eastern Pacific
    'wp': Western Pacific
    'ni': North Indian Ocean
    'si': Southwest Indian Ocean
    'au': Australian Region
    'sp': South Pacific
    
    Note: months [1,2,3,4,5,6] will obtain storms 
    that *formed* within the time period 1 Jan to
    30 June inclusive. 
    
    Setting months [11,12,1,2,3,4] and years [1996]
    will return storms that formed between 1 Nov
    1996 and 30 April 1997 inclusive.
    
    Setting the basin will return any storms which
    *passed through* (not just formed) in the 
    designated basin region. 
    
    """    
    lats, lons, count = tc_assess_code.storm_lats_lons(storms, years, months, basin, genesis=genesis)
    cube = tc_assess_code._binned_cube(lats, lons)
    if nyears == 0:
        nyears = float(len(years))
    cube /= nyears
    cube /= len(months)
    if data_freq != 6:
        cube *= float(data_freq) / 6
    print('Total number of storms in time period:', years, count)
    return cube, count
        
def calc_and_write_interannual_variability(runid_info, storms, years, scratch_dir, months_nh, months_sh):
    for ir, suiteid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            ens = runid_info['ens'][ir]
            expt = runid_info['experiment'][ir]
            key = define_key(algo, expt, ens)
            if ens == '':
                outdir = scratch_dir.format(suiteid, algo)
            else:
                outdir = scratch_dir.format(suiteid+'/'+ens, algo)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            storms_all = storms[(suiteid, key)]
            years = storms[(suiteid, key, 'years')]
            for year in years[:-1]:
                for hemi in ['nh', 'sh']:
                    fname_out = os.path.join(outdir, 'tc_count_ace_'+str(year)+'_'+hemi+'.nc')
                    #print('do tc_count_ace ',fname_out)
                    if not os.path.exists(fname_out):
                        if hemi == 'nh':
                            months = months_nh
                        else:
                            months = months_sh
                        
                        storms_year = tc_assess_code.get_storms_in_year(storms_all, year, months)
                        basins = []
                        for b in BASINS_LIST:
                            if BASINS[b][0] == hemi:
                                basins.append(b)
                        annual_count, no_storms = tc_assess_code.calc_interann.calculate_model_count_ace(suiteid, storms_year, year, basins, months, hemi, nensemble=1, ace=False)
                        annual_ace, no_storms_ace = tc_assess_code.calc_interann.calculate_model_count_ace(suiteid, storms_year, year, basins, months, hemi, nensemble=1, ace=True)
                        #print('annual count, no storms ',annual_count[hemi], no_storms, no_storms_ace)
                        if not no_storms and not no_storms_ace:
                            # since the year/month is given, the time_units and calendar do not really matter as this is the measure during a period
                            time_units = 'hours since 1950-01-01 00:00:00'
                            calendar = 'gregorian'
                            start_date = str(year)+str(months[0]).zfill(2)+str(1).zfill(2)
                            if hemi == 'nh':
                                end_date = str(year)+str(months[-1]).zfill(2)+str(31).zfill(2)
                            else:
                                end_date = str(year+1)+str(months[-1]).zfill(2)+str(31).zfill(2)
                            algorithm = 'TempestExtremes'
                            tc_assess_code.calc_interann.write_tc_count_ace_netcdf(fname_out, basins, annual_count, annual_ace, year, time_units, calendar, 
                              start_date, end_date, algorithm, 
                                                                               suiteid, hemi=hemi)

def do_track_density_calc(runid_info, storms, years, plot_dir, scratch_dir, months_nh, months_sh, genesis=False, keys = []):
    '''
    Calculate the track density for each year and save
    '''
    if genesis:
        fname = 'genesis'
    else:
        fname = 'density'

    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            ens = runid_info['ens'][ir]
            if ens == '':
                store_dir = scratch_dir.format(runid, algo)
            else:
                store_dir = scratch_dir.format(runid+'/'+ens, algo)
            if not os.path.exists(store_dir):
                os.makedirs(store_dir)
            key = keys[ir]
            storms_runid = storms[(runid, key)]
            nyears = storms[(runid, key, 'nyears')]
            years_model = storms[(runid, key, 'years')]
            data_freq = storms[(runid, key, 'freq')]
            for year in years_model:
                this_year = [year, year]
                fname_out = os.path.join(store_dir, fname+'_'+str(year)+'.nc')
                if not os.path.exists(fname_out):
                    print('do den calc ',runid, algo, year)
                    cube_nh, count_nh = storm_density_calculation(storms_runid, this_year, months_nh, 'nh', genesis=genesis, nyears = 1, data_freq=data_freq)
                    cube_sh, count_sh = storm_density_calculation(storms_runid, this_year, months_sh, 'sh', genesis=genesis, nyears = 1, data_freq=data_freq)
                    if count_nh > 0 and count_sh > 0:
                        cube = cube_nh + cube_sh
                        ncube = tc_assess_code.add_time_coord_to_cube(cube, year, months_nh, months_sh, fname)
                        iris.save(ncube, os.path.join(store_dir, fname+'_'+str(year)+'.nc'))
            #cube_density[(runid, key)] = cube
            #iris.save(cube,'./test_density.nc')

def calculate_tc_properties(runid_info, storms, years, storedir, plot_dir, keys, months_nh, months_sh):
    '''
    Do calculation of TC properties such as density, variability,
and save for later reading/plotting
    '''
    ## main part of calculation
    time_calc_var_start = time.perf_counter()
    calc_and_write_interannual_variability(runid_info, storms, years, storedir, months_nh, months_sh)

    time_calc_var = time.perf_counter()
    print('Calc var took ',(time_calc_var - time_calc_var_start),' seconds')

    do_track_density_calc(runid_info, storms, years, plot_dir, storedir, months_nh, months_sh, genesis = False, keys = keys)
    do_track_density_calc(runid_info, storms, years, plot_dir, storedir, months_nh, months_sh, genesis = True, keys = keys)

    time_calc_den = time.perf_counter()
    print('Calc den took ',(time_calc_den - time_calc_var),' seconds')

