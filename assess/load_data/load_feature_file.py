'''
Routines to read track files
Different formats:
  TRACK ascii - what are the columns?
  TRACK netcdf
  TempestExtremes ascii
  TempestExtremes netcdf
'''

import sys, os, subprocess, glob
import netCDF4
sys.path.append('/data/users/hadom/tracking')
import ts_obs
import ts_obs.load_data

sys.path.append('/data/users/hadom/branches/git/track_analysis')
import load_data.load_TE_nc as load_TE_nc
import load_data.load_TRACK_ascii as load_TRACK_ascii

def get_data_years_from_mass(suiteid, runid, ens, search, dir_final):
    '''
    Extract the track data from MASS
    '''
    if ens == '':
        moodir = 'moose:/crum/'+suiteid+'/any.nc.file/'
    else:
        moodir = 'moose:/ens/'+suiteid+'/'+ens+'/any.nc.file/'
    cmd =  'moo ls '+moodir+search+'  | head -1 | cut -d _ -f 4 | cut -c 1-4'
    print(cmd)
    sts=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
    year_start = sts.stdout
    cmd =  'moo ls '+moodir+search+'  | tail -1 | cut -d _ -f 4 | cut -c 1-4'
    print(cmd)
    sts=subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)
    year_end = sts.stdout
    return int(year_start), int(year_end)

def get_data_from_mass(suiteid, runid, ens, search, dir_final):
    '''
    Extract the track data from MASS
    '''
    if ens == '':
        moo_path = 'moose:/crum/'+suiteid+'/any.nc.file/'+search
    else:
        moo_path = 'moose:/ens/'+suiteid+'/'+ens+'/any.nc.file/'+search
    cmd = 'moo get -i '+moo_path+' '+dir_final
    print(cmd)
    os.system(cmd)

def read_storms_tempest_nc(dir_in, runid_info, suiteid, runid, ens, years, algo, algo_type, filepattern=''):
    '''
    For a given file name format, and directory
    Read the storms from the files
    If vort is in algo_type, then need to read both nh and sh parts
    '''

    if filepattern == '':
        filepattern = "TC_{year}_{resol}_u-{runid}.{track_type}.{hemi}{track_m}.tracks"

    dir_track = dir_in
    if not os.path.exists(dir_track):
        os.makedirs(dir_track)

    if 'vort' in algo or 'ew' in algo:
        algo_search = algo+'_?h'
    else:
        algo_search = algo

    search = filepattern.format(runid=runid, year='*', algo_type=algo_type)

    fnames = []; storms = []
    no_years_with_storms = 0
    print('get/read over years ',years)
    fnames = []
    years_exist = []
    for year in years:
        print('do year ',year)
        search_fname = os.path.join(dir_track, filepattern.format(runid=runid, year=year, algo_type=algo_type))
        print ('search_fname ',search_fname)

        if os.path.exists(search_fname):
            f = search_fname
            storms_list = list(load_TE_nc.load_cmor(f, vort_variable = 'rvT63'))

            storms.extend(storms_list)
            print('no storms, year ',len(storms_list), year)
            if len(storms_list) > 0:
                no_years_with_storms += 1
                fnames.append(f)
                years_exist.append(year)
                
    # check the metadata to discover which algorithm this is, and hence
    # what feature variable is tracked
    with netCDF4.Dataset(fnames[0], 'r') as nc:
        track_algorithm = nc.getncattr('algorithm')
        if track_algorithm == 'TRACK':
            track_extra = nc.getncattr('algorithm_extra')
            if track_extra == 'T63avg':
                feature_variable = 'vortmean_T63'
            else:
                feature_variable = 'rv850_T42'
        elif 'TempestExtremes' in track_algorithm:
            feature_variable = 'psl'
            track_extra = ''
        else:
            raise Exception('Unrecognised algorithm in netcdf file '+fnames[0])

    return storms, feature_variable, track_extra, no_years_with_storms, years_exist[0], years_exist[-1]

def read_storms_track_nc(dir_in_base, runid_info, suiteid, runid, ens, years, algo, algo_type, filepattern=''):
    '''
    For a given file name format, and directory
    Read the storms from the files
    If vort is in algo_type, then need to read both nh and sh parts
    '''

    if filepattern == '':
        filepattern = "TC_{year}_{resol}_u-{runid}.{track_type}.{hemi}{track_m}.tracks"

    dir_track = dir_in_base
    if not os.path.exists(dir_track):
        os.makedirs(dir_track)

    fnames = []; storms = []
    no_years_with_storms = 0
    
    print('get/read over years ',years)
    fnames = []
    years_exist = []
    for year in years:
        if year > year_end:
            continue
        if year < year_start:
            continue
        print('do year ',year)
        search_fname = os.path.join(dir_track, filepattern.format(runid, str(year), algo_type=algo_type))
        print ('search_fname ',search_fname)
        if not os.path.exists(search_fname):
            print('file does not exist ',search_fname)

        if os.path.exists(search_fname):
            f = search_fname
            storms_list = list(load_TE_nc.load_cmor(f, vort_variable = 'T42'))

            storms.extend(storms_list)
            print('no storms, year ',len(storms_list), year)
            if len(storms_list) > 0:
                no_years_with_storms += 1
                fnames.append(f)
                years_exist.append(year)
                
    # check the metadata to discover which algorithm this is, and hence
    # what feature variable is tracked
    with netCDF4.Dataset(fnames[0], 'r') as nc:
        track_algorithm = nc.getncattr('algorithm')
        if track_algorithm == 'TRACK':
            track_extra = nc.getncattr('algorithm_extra')
            if track_extra == 'T63avg':
                feature_variable = 'vortmean_T63'
            else:
                feature_variable = 'rv850_T42'
        elif 'TempestExtremes' in track_algorithm:
            feature_variable = 'psl'
            track_extra = ''
        else:
            raise Exception('Unrecognised algorithm in netcdf file '+fnames[0])

    return storms, feature_variable, track_extra, no_years_with_storms, years_exist[0], years_exist[-1]

def read_storms_from_highresmip(hemi, run, years, suiteid, runid, algo, experiment = 'highresSST-present'):
    '''
    For a given file name format, and directory
    Read the storms from the netcdf files
    If vort is in algo_type, then need to read both nh and sh parts
    This is the format used for the files uploaded to the CEDA archive
    '''
    rip = 'r1i1p1f1'
    dir_in = '/data/users/hadom/PRIMAVERA/model_derived_data/storm_tracking/{}/MOHC/{}/{}/'+rip+'/tropical/v3/'
    fname_nc = 'TC-?H_{}_{}_{}_'+rip+'_gn_{}0101-{}1231.nc4'

    fnames = []; storms = []
    no_years_with_storms = 0
    for year in years:
        search = fname_nc.format(algo, runid, experiment, str(year), str(year))
        print ('search fname ',os.path.join(dir_in.format(algo, runid, experiment), search))
        #search = fname_nc.format(run['algo_type'])
        fname = glob.glob(os.path.join(dir_in.format(algo, runid, experiment), search))
        if fname:
            fnames.append(fname[0])
            for nf, f in enumerate(fname):
                print ('read f ',f)
                storms_list = list(track_cmor.load_cmor(f, vort_variable = 'rvT63'))
                storms.extend(storms_list)
                print('no storms, year ',len(storms_list), year)
                if nf == 0 and len(storms_list) > 0:
                    no_years_with_storms += 1
                
    print ('fnames ',fnames)
    # check the metadata to discover which algorithm this is, and hence
    # what feature variable is tracked
    with netCDF4.Dataset(fnames[0], 'r') as nc:
        track_algorithm = nc.getncattr('algorithm')
        if track_algorithm == 'TRACK':
            track_extra = nc.getncattr('algorithm_extra')
            if track_extra == 'T63avg':
                feature_variable = 'vortmean_T63'
            else:
                feature_variable = 'rv850_T42'
        elif 'TempestExtremes' in track_algorithm:
            feature_variable = 'slp'
            track_extra = ''
        else:
            raise Exception('Unrecognised algorithm in netcdf file '+path)

    return storms, feature_variable, track_extra, no_years_with_storms

def read_storms_using_track(hemi, run, years, suiteid, runid, algo, track_method = 'T63'):
    '''
    Read TRACK combined file, will have to guess what variables are in columns
    '''

    years_str = str(years[0])+'-'+str(years[-1])
    
    #fname_n = "combined_UVT_%s_%s%s_%s%s.NH%s.tracks"
    #fname_s = "combined_UVT_%s_%s%s_%s%s.SH%s.tracks"
    fname_n = "combined_UVT_*%s.NH%s.tracks"
    fname_s = "combined_UVT_*%s.SH%s.tracks"

    dirs_in = RESULTS_DIR.format('u-'+runid)
    track_type = run['track_type']
    #track_method = 'T42'
    if track_method == 'T42':
        track_m = ''
    elif track_method == 'T63':
        track_m = '.T63'
    fname_nh = fname_n % (track_type, track_m)
    fname_sh = fname_s % (track_type, track_m)

    storms = []
    for hemi in ['nh','sh']:
        if hemi == 'nh':
            fname = fname_nh
        else:
            fname = fname_sh
        search = os.path.join(dirs_in, fname)
        files = sorted(glob.glob(search))
        print ('search ',search)
        if search:
            path = files[-1]
        else:
            print ('no data in search ',search)
        print ('path for track data ',path)
        storms_list = list(ts_model.track.load(path, ex_cols = 3, calendar='360_day'))
        storms.extend(storms_list)

    storms_filtered = filter_time(storms, years, months = range(1,13))        
    no_years_with_storms = len(years)
    print('algo, no storms, no years ',len(storms_filtered), no_years_with_storms)
    if len(storms_filtered) == 0:
        raise Exception('No storms '+ runid)
    return storms_filtered, no_years_with_storms

def read_storms_using_track_ascii_years(runid_info, years, suiteid, runid, algo, track_method = 'T63', resol = 'N216e', fileformat='ascii', filepattern='', variable_indices={}, calendar='gregorian'):
    '''
    Read TRACK ascii individual years
    Will have to guess what variables are in columns
    '''

    years_str = str(years[0])+'-'+str(years[-1])
    
    #fname_pattern = "TC_{year}_{resol}_u-{runid}.{track_type}.{hemi}{track_m}.tracks"

    if filepattern == '':
        fnamepattern = "TC_{year}_{resol}_u-{runid}.{track_type}.{hemi}{track_m}.tracks"
    print('filepattern ',filepattern)

    dirs_in = os.path.join(runid_info['TRACK_results_dir'], 'u-'+runid)

    if track_method == 'T42':
        track_m = ''
    elif track_method == 'T63':
        track_m = '.T63'
    elif track_method == 'T63full':
        track_m = '.T63full'
    elif track_method == 'T42new':
        track_m = ''
    else:
        track_m = ''
    track_type = runid_info['track_type']

    storms = []
    years_with_data = []
    no_years_with_storms = 0
    for year in years:

        for hemi in ['nh','sh']:
            if hemi == 'nh':
                hemipos = 'pos'
            else:
                hemipos = 'neg'
            fname = filepattern.format(runid=runid, resol=resol, year=year, track_type=track_type, hemi=hemi.upper(), track_m=track_m, hemipos=hemipos)
            search = os.path.join(dirs_in, fname)
            files = sorted(glob.glob(search))
            print ('search ',search)
            if len(files) == 1:
                years_with_data.append(year)
                path = files[-1]
                print ('path for track data ',path, fileformat)
                if fileformat == 'ascii':
                    storms_list = list(load_TRACK_ascii.load(path, ex_cols = 3, calendar=calendar, variable_indices=variable_indices))
                else:
                    storms_list = list(load_TE_nc.load_cmor(path, vort_variable = track_method))
                storms.extend(storms_list)
                if hemi == 'nh':
                    no_years_with_storms += 1
                    #print('fname ',fname)
                    #print('storms ',storms_list[:1])
                    #for storm in storms_list[0:2]:
                    #    print('storm ',storm.genesis_date(), storm.obs_at_vmax().mslp, storm.obs_at_vmax().vmax, storm.obs_at_vmax().extras['w10m'])
                    #    for ob in storm.obs:
                    #        print('lat, lon, mslp, vmax ',ob.lat, ob.lon, ob.mslp, ob.vmax, ob.extras)
                    #        stop
            else:
                print ('no data in search ',search)

    year_start = int(min(years_with_data))
    year_end = int(max(years_with_data))
    print('algo, no storms, no years ',len(storms), no_years_with_storms)
    return storms, no_years_with_storms, year_start, year_end

def read_storms_with_appropriate_reader(suiteid, runid, ens, runid_info, algo, algo_type, expt, key, resol_in, track_method, dir_in, years, institute, rip, grid, fileformat='nc', filepattern='', variable_indices={}, calendar='gregorian'):
    '''
    Currently assuming tropical cyclone input dataset
    Identify what kind of input file is going to be specified, and hence
    how to read the data
    Netcdf file - should be able to read TRACK and TE in similar way
    TRACK ascii - need to know/guess variables
    TE ascii - need to know/guess variables
    '''

    if algo == 'TRACK' or algo == 'TempExt':
        print('use read_storms_from_highresmip')
        storms_runid, feature_variable, track_extra, no_years_with_storms_runid, year_start, year_end = read_storms_from_highresmip(dir_in, basin, runid_info, years, suiteid, runid, algo, institute, rip, grid, experiment = expt)

    elif 'track' in algo:
        if fileformat == 'ascii':
            print('use read_storms_using_track_years, ascii')
            storms_runid, no_years_with_storms_runid, year_start, year_end = read_storms_using_track_ascii_years(runid_info, years, suiteid, runid, algo, track_method=track_method, resol = resol_in, fileformat = fileformat, filepattern=filepattern, variable_indices=variable_indices, calendar=calendar)
        else:
            print('use read_storms_using_track_years, nc')
            storms_runid, no_years_with_storms_runid, year_start, year_end = read_storms_using_track_ascii_years(runid_info, years, suiteid, runid, algo, track_method=track_method, resol = resol_in, fileformat = fileformat, filepattern=filepattern, variable_indices=variable_indices, calendar=calendar)

    else:
        print('use read_storms_tempest_nc')
        if fileformat == 'ascii':
            storms_runid, feature_variable, track_extra, no_years_with_storms_runid, year_start, year_end = read_storms_tempest_nc(dir_in, runid_info, suiteid, runid, ens, years, algo, algo_type, filepattern=filepattern)
        else:
            storms_runid, feature_variable, track_extra, no_years_with_storms_runid, year_start, year_end = read_storms_tempest_nc(dir_in, runid_info, suiteid, runid, ens, years, algo, algo_type, filepattern=filepattern)

    return storms_runid, no_years_with_storms_runid, year_start, year_end


