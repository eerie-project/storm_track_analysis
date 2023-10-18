import os

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

def work(suiteid, runid, ens, runid_info, algo, algo_type, expt, key, resol, track_method, dir_in, years, institute, rip, grid, fileformat='nc', filepattern=''):
    '''
    For a given file name format, and directory
    Get tracks from MASS archive
    '''

    if filepattern == '':
        filepattern = "TC_{year}_{resol}_u-{runid}.{track_type}.{hemi}{track_m}.tracks"
    dir_track = dir_in
    if not os.path.exists(dir_track):
        os.makedirs(dir_track)

    algo_search = algo

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

    fname = filepattern.format(runid=runid, resol=resol, year='*', track_type=track_type, hemi='*', track_m=track_m, hemipos='*', algo_type=algo_type)

    if years == []:
        year_start, year_end = get_data_years_from_mass(suiteid, runid, ens, search, dir_track)
        years = np.arange(year_start, year_end+1)
    else:
        year_start = years[0]
        year_end = years[-1]

    print('get years ',years)
    for year in years:
        if year > year_end:
            continue
        if year < year_start:
            continue
        print('do year ',year)
        fname = filepattern.format(runid=runid, resol=resol, year=year, track_type=track_type, hemi='*', track_m=track_m, hemipos='*', algo_type=algo_type)
        search_fname = os.path.join(dir_track, fname)
        print ('search_fname ',search_fname)
        if not os.path.exists(search_fname):
            get_data_from_mass(suiteid, runid, ens, fname, dir_track)

