import sys
import numpy as np
import os.path
import collections

import cftime
import cf_units
import fnmatch

import cartopy.crs as ccrs
import shapely.geometry as sgeom
import calendar

import iris
import iris.coord_systems as icoord_systems
import iris.coords as icoords

tropical_storm_threshold = {'10mwind': 18.0, 'mslp': 1005.0}
knots_to_ms = 1.0 / 1.944

#: Lat/lon locations for each ocean basin for mapping. If set to 
#: None then it returns a global map
MAP_REGION = {'na': (-115, 0, 60, 0),
              'ep': (-180, -70, 40, 0),
              'wp': (-265, -140, 50, 0),
              'cp': (-190, -100, 50, 0),
              'ni': (-310, -260, 30, 0),
              'si': (-340, -260, 0, -40),
              'au': (-270, -195, 0, -40),
              'sp': (-200, -100, 0, -40),
              'sa': ( -90, 0, 0, -40),
              'nh': (-180, 30, 70, 180),
              'sh': (-180, 30, 180, -90),
              None: (-180, 180, 90, -90)
              }

#: Lat/lon locations of tracking regions for each ocean basin. If None
#: then returns a region for the whole globe
TRACKING_REGION = {'na': ([-75, 20, 20, -80, -80, -100, -100, -75, -75], [0, 0, 60, 60, 40, 40, 20, 6, 0]),
                   'ep': ([-140, -75, -75, -100, -100, -140, -140], [0, 0, 6, 20, 60, 60, 0]),
                   #'wp': ([-260, -180, -180, -260, -260], [0, 0, 60, 60, 0]),
#                   'wp': ([-260, -180, -180, -260, -260], [0, 0, 60, 60, 0]),
                   'wp': ([100, 179.9, 179.9, 100, 100], [0, 0, 60, 60, 0]),
                   'cp': ([-180, -140, -140, -180, -180], [0, 0, 50, 50, 0]),
#                   'ni': ([-320, -260, -260, -320, -320], [0, 0, 30, 30, 0]),
#                   'si': ([-330, -270, -270, -330, -330], [-40, -40, 0, 0, -40]),
#                   'au': ([-270, -200, -200, -270, -270], [-40, -40, 0, 0, -40]),
#                   'sp': ([-200, -120, -120, -200, -200], [-40, -40, 0, 0, -40]),
                   'ni': ([40, 100, 100, 40, 40], [0, 0, 30, 30, 0]),
                   'si': ([30, 90, 90, 30, 30], [-40, -40, 0, 0, -40]),
                   'au': ([90, 160, 160, 90, 90], [-40, -40, 0, 0, -40]),
                   'sp': ([160, 240, 240, 160, 160], [-40, -40, 0, 0, -40]),
                   'sa': ([-90, 0, 0, -90, -90], [-40, -40, 0, 0, -40]),
#                   'nh': ([-360, 0, 0, -360, -360],[0, 0, 90, 90 ,0]),
#                   'nh': ([-180, 180, 180, -180, -180],[0, 0, 90, 90 ,0]),
#                   'sh': ([-360, 0, 0, -360, -360],[-90, -90, 0, 0 ,-90]),
#                   'sh': ([-180, 180, 180, -180, -180],[0, 0, -90, -90 ,0]),
                   'sh': ([-180, 179.9, 179.9, -180, -180],[-90, -90, 0, 0 ,-90]),
                   'nh': ([-180, 179.9, 179.9, -180, -180],[0, 0, 90, 90 , 0]),
                   'mdr': ([-80, -20, -20, -80, -80], [10, 10, 20, 20, 10]),
                   'wc': ([-360, 0, 0, -360, -360],[0, 0, 30, 30 ,0]),
                   'ena': ([-60, 30, 30, -60, -60], [5, 5, 40, 40, 5]),
                   'wna': ([-75, -60, -60, -80, -80, -100, -100, -75, -75], [0, 0, 60, 60, 40, 40, 20, 6, 0]),
                   None: ([-360, 0, 0, -360, -360],[-90, -90, 90, 90 ,-90])
                   }    

#: Corresponding full basin names for each abbreviation
BASIN_NAME = {'na': 'North Atlantic',
              'ep': 'Eastern Pacific',
              'wp': 'Western Pacific',
              'cp': 'Central Pacific',
              'ni': 'North Indian Ocean',
              'si': 'Southwest Indian Ocean',
              'au': 'Australian Region',
              'sp': 'South Pacific',
              'nh': 'Northern Hemisphere',
              'sh': 'Southern Hemisphere',
              'ena': 'E NAtlantic',
              'wna': 'W NAtlantic',
              None: 'Global'
              }
#: Corresponding full basin names for each abbreviation
BASIN_NAME_SHORT = {'na': 'NA',
              'ep': 'EP',
              'wp': 'WP',
              'cp': 'CP',
              'ni': 'NI',
              'si': 'SWI',
              'au': 'AU',
              'sp': 'SP',
              'nh': 'NHe',
              'sh': 'SH',
              'ena': 'ENA',
              'wna': 'WNA',
              None: 'Global'
              }

BASINS_LIST = ['nh','sh','na','ep','wp','cp','ni','si','au','sp','sa'] 

BASINS = {'nh': ['nh','N Hemi', 0, 120, 0, 900],'sh': ['sh','S Hemi', 0, 110, 0, 600],'na':['nh','N Atl', 0, 25, 0, 250],'ep': ['nh','E Pac', 0, 30, 0, 350],\
          'wp':['nh','W Pac', 0, 55, 0, 600],\
          'cp':['nh','C Pac', 0, 15, 0, 100],'ni': ['nh','N Ind', 0, 15, 0, 50],'si': ['sh','S Ind', 0, 40, 0, 200],'au': ['sh','Aust', 0, 25, 0, 200],\
          'sp': ['sh','S Pac', 0, 50, 0, 200],'sa': ['sh','S Atl', 0, 10, 0, 50]}

#: Corresponding month name for a given integer value
NUM_TO_MONTH = {1: 'Jan', 2: 'Feb', 3: 'Mar', 4: 'Apr', 5: 'May', 6: 'Jun',
                7: 'Jul', 8: 'Aug', 9: 'Sep',10: 'Oct',11: 'Nov',12: 'Dec'}

def _month_names(months):
    """ Returns list of month names for a given set of integer values """
    names = []
    for month in months:
        names.append(NUM_TO_MONTH.get(month))
    return names
        

def add_time_coord_to_cube(cube, year, months_nh, months_sh, fname):
    '''
    Add time coordinate to track density/genesis cube
    Make this the start of the period i.e. NH year/first month_nh
    Add time period to cube attributes
    '''
    i = cftime.datetime(year, months_nh[0], 1,0,0,0)
    calendar = 'gregorian'
    calendar_unit = 'days since 1950-01-01 00:00:00'
    time_points = cf_units.date2num(i, calendar_unit, calendar=calendar)    
    time_coord = iris.coords.DimCoord(time_points, standard_name='time',
                                     units=cf_units.Unit(calendar_unit,
                                                 calendar=calendar))
    cube.add_aux_coord(time_coord)
    ncube = iris.util.new_axis(cube, 'time')
    ncube.attributes['start_date_nh'] = str(year)+str(months_nh[0]).zfill(2)+'01'
    ncube.attributes['end_date_nh'] = str(year)+str(months_nh[-1]).zfill(2)+'31'
    ncube.attributes['start_date_sh'] = str(year)+str(months_sh[0]).zfill(2)+'01'
    ncube.attributes['end_date_sh'] = str(year+1)+str(months_sh[-1]).zfill(2)+'31'
    ncube.var_name =fname
    
    return ncube

def _storms_in_time_range(storms, year, months):
    """Returns a generator of storms that formed during the desired time period """
    calendar = 'proleptic_gregorian'
    for storm in storms[:1]:
        # derive the calendar from the storm object, and then pass this to ensure that the start/end period has the same calendar for comparison
        cal_type = str(type(storm.genesis_date()))
        cal = cal_type.split('.')[-1][8:-2]
        #print ('calendar ',cal_type, cal)
        if '360' in cal:
            calendar = '360_day'
        elif '365' in cal:
            calendar = '365_day'
        elif 'noleap' in cal or 'NoLeap' in cal:
            calendar = 'noleap'
        elif cal == 'Gregorian':
            calendar = 'gregorian'
        else:
            calendar = 'proleptic_gregorian'

    start_date, end_date = _get_time_range(year, months, calendar=calendar)
    for storm in storms:        
        if (storm.genesis_date() >= start_date) and (storm.genesis_date() < end_date):
            yield storm

def _storm_vmax_in_basin(storm, basin):
    """ Returns True if the maximum intensity of the storm occurred
    in desired ocean basin. """
    rbox = _basin_polygon(basin)  
    if 'obs_at_vmax' in dir(storm):
        xy = ccrs.PlateCarree().transform_point(storm.obs_at_vmax().lon, storm.obs_at_vmax().lat, ccrs.Geodetic())
    else:
        xy = ccrs.PlateCarree().transform_point(storm.obs_at_vmax().lon, storm.obs_at_vmax().lat, ccrs.Geodetic())
    point = sgeom.Point(xy[0], xy[1])
    if point.within(rbox):
        return True
    return False

def _storms_in_strength_range(storms, year, months, storm_types):
    """
    Returns a generator of storms that formed during the desired time 
    period.
    
    """
    for storm in storms:        
        if (storm.max_storm_type() in storm_types):
            yield storm
            
def _basin_polygon(basin, project=True):
    """ 
    Returns a polygon of the tracking region for a particular 
    ocean basin. i.e. a storm must pass through this region 
    in order to be retined. For example, if basin is set to 
    'au' then storms for the Australian region must pass
    through the area defined by -270 to -200W, 0 to -40S.
    
    """
    rbox = sgeom.Polygon(list(zip(*TRACKING_REGION.get(basin))))
    if project: 
        rbox = ccrs.PlateCarree().project_geometry(rbox, ccrs.PlateCarree())
    return rbox
    
def _storm_in_basin(storm, basin):
    """ Returns True if a storm track intersects a defined ocean basin """
    rbox = _basin_polygon(basin)   
    lons, lats = list(zip(*[(ob.lon, ob.lat) for ob in storm.obs]))
    track = sgeom.LineString(list(zip(lons, lats)))       
    projected_track = ccrs.PlateCarree().project_geometry(track, ccrs.Geodetic())
    if rbox.intersects(projected_track):
        return True
    return False

def storm_in_time_range(storm, years, months):
    """Returns True if a storm formed during a given years and months. """
    cal_type = str(type(storm.genesis_date()))
    cal = cal_type.split('.')[-1][8:]
    if '360' in cal:
        calendar = '360_day'
    elif '365' in cal:
        calendar = '365_day'
    elif 'noleap' in cal or 'NoLeap' in cal:
        calendar = 'noleap'
    else:
        calendar = 'proleptic_gregorian'

    for year in years:
        start_date, end_date = _get_time_range(year, months, calendar = calendar)   
        if (storm.genesis_date() >= start_date and 
            storm.genesis_date() < end_date):
            return True

def _get_genesis_months(storms, years, basin):
    """ 
    Returns genesis month of all storms that formed within a 
    given set of years 
    
    """
    genesis_months = []
    for storm in storms:
        if (storm.genesis_date().year in years) and _storm_vmax_in_basin(storm, basin):
            genesis_months.append(storm.genesis_date().month)
    return genesis_months
            
def _get_monthly_storm_count(storms, years, months, basin):
    """ Returns list of storm counts for a desired set of months """
    genesis_months = _get_genesis_months(storms, years, basin)
    monthly_count = []
    for month in months:
        monthly_count.append(genesis_months.count(month))
    return monthly_count

def get_monthly_mean_count(storms, years, months, basin, annual = False, nyears = 0):
    """ 
    Returns list of monthly mean storm counts for a given
    set of years and months 
    """
    monthly_count = _get_monthly_storm_count(storms, years, months, basin)
    if nyears == 0:
        nyears = float(len(years))
    mean_monthly_count = []
    for count in monthly_count:
        mean_monthly_count.append(float(count)/nyears)

    if annual:
        return sum(mean_monthly_count)
    else:
        return mean_monthly_count

def _get_annual_vmax_storm_count(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    count = 0
    for storm in storms:
        if (storm.genesis_date().year in years) and _storm_vmax_in_basin(storm, basin):
            count += 1
    storm_counts.append(count)
    return storm_counts    

def _get_annual_vmax_storm_count_category(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = np.zeros(6)
    count = np.zeros(6)
    for storm in _storms_in_time_range(storms, years[0], months):
        if _storm_vmax_in_basin(storm, basin):
            category = storm_intensity(storm.obs_at_vmax().mslp)
            count[category] += 1
        storm_counts=count
    return storm_counts    

def _get_annual_vmax_storm_ace_category(storms, years, months, basin):
    """ 
    Returns array of storm ACE for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = np.zeros(6); storm_ace = np.zeros(6)
    count = np.zeros(6); ace = np.zeros(6)
    for storm in _storms_in_time_range(storms, years[0], months):
        if _storm_vmax_in_basin(storm, basin):
            category = storm_intensity(storm.obs_at_vmax().mslp)
            count[category] += 1
            #ace[category] += storm.ace_index_no6hrcheck()
            ace[category] += storm.ace_index()
        storm_ace=ace
    return storm_ace   

def get_annual_vmax_mean_count_category(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean storm counts for a given
    set of years and months 
    
    
    """
    annual_category = np.zeros(len(years)*6).reshape(len(years),6)
    mean_annual_category = np.zeros(6)
    for iyr, year in enumerate(years):
        annual_count = _get_annual_vmax_storm_count_category(storms, [year,year], months, basin)
        for cat in range(0,6):
            annual_category[iyr,cat] = float(annual_count[cat])/float(nensemble)
            mean_annual_category[cat] = mean_annual_category[cat] + annual_category[iyr,cat]
    return annual_category

def _get_annual_vmax_storm_ace_category(storms, years, months, basin):
    """ 
    Returns array of storm ACE for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = np.zeros(6); storm_ace = np.zeros(6)
    count = np.zeros(6); ace = np.zeros(6)
    for storm in _storms_in_time_range(storms, years[0], months):
        if _storm_vmax_in_basin(storm, basin):
            category = storm_intensity(storm.obs_at_vmax().mslp)
            count[category] += 1
            #ace[category] += storm.ace_index_no6hrcheck()
            ace[category] += storm.ace_index()
        storm_ace=ace
    return storm_ace   

def get_annual_vmax_mean_ace_category(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean storm ACE for a given
    set of years and months 
    
    
    """
    annual_category = np.zeros(len(years)*6).reshape(len(years),6)
    mean_annual_category = np.zeros(6)
    for iyr, year in enumerate(years):
        annual_ace = _get_annual_vmax_storm_ace_category(storms, [year,year], months, basin)
        for cat in range(0,6):
            annual_category[iyr,cat] = float(annual_ace[cat])/float(nensemble)
            mean_annual_category[cat] = mean_annual_category[cat] + annual_category[iyr,cat]
    return annual_category

def annual_vmax_storm_counts_category(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    annual_category = np.zeros(len(years)*6).reshape(len(years),6)
    for iyr, year in enumerate(years):
        count = 0
        for storm in _storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
#                category = storm_intensity(storm.obs_at_vmax().mslp)
#                category2 = storm_intensity_vmax(storm.obs_at_vmax().vmax)
                category = storm_intensity(storm.obs_at_vmax().mslp)
                category2 = storm_intensity_vmax(storm.obs_at_vmax().vmax)
                annual_category[iyr,category2] += 1
    return annual_category

def annual_vmax_storm_ace_category(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    annual_category = np.zeros(len(years)*6).reshape(len(years),6)
    for iyr, year in enumerate(years):
        count = 0
        for storm in _storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
                category = storm_intensity(storm.obs_at_vmax().mslp)
                category2 = storm_intensity_vmax(storm.obs_at_vmax().vmax)
                annual_category[iyr,category2] += storm.ace_index()
    return annual_category

def _get_annual_vmax_storm_count_hemi(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    count = 0
#    print years, months
    for storm in _storms_in_time_range(storms, years[0], months):
        if _storm_vmax_in_basin(storm, basin):
            count += 1
    storm_counts.append(count)
    return storm_counts    

def _get_annual_vmax_storm_ace_hemi(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_ace = []
    ace = 0
#    print years, months
    for storm in _storms_in_time_range(storms, years[0], months):
        #print 'found storm in ace, anywhere ', storm.ace_index()
        if _storm_vmax_in_basin(storm, basin):
            #ace += storm.ace_index_no6hrcheck()
            ace += storm.ace_index()
    storm_ace.append(ace)
    return storm_ace   

def get_annual_vmax_mean_count(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean storm counts for a given
    set of years and months 
    
    
    """
    mean_annual_count = []
    for year in years:
        annual_count = _get_annual_vmax_storm_count_hemi(storms, [year,year], months, basin)
        for count in annual_count:
            mean_annual_count.append(float(count)/float(nensemble))
    return mean_annual_count

def get_annual_vmax_mean_ace(storms, years, months, basin, nensemble):
    """ 
    Returns list of annual mean ace for a given
    set of years and months 
    """
    mean_annual_ace = []
    for year in years:
        annual_ace = _get_annual_vmax_storm_ace_hemi(storms, [year,year], months, basin)
#        print 'year, annual ace ',annual_ace
        for ace in annual_ace:
            mean_annual_ace.append(float(ace)/float(nensemble))
    return mean_annual_ace

def annual_vmax_storm_counts(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_counts = []
    for year in years:
        count = 0
        for storm in _storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
                count += 1
        storm_counts.append(count)
    return storm_counts   

def annual_vmax_storm_ace(storms, years, months, basin, storm_types=['SS', 'TS', 'HU', 'MH']
):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    storm_ace = []
    for year in years:
        ace = 0
        for storm in _storms_in_time_range(storms, year, months):
            if (storm.max_storm_type() in storm_types) and _storm_vmax_in_basin(storm, basin):
                ace += storm.ace_index()
        storm_ace.append(ace)
    return storm_ace  

def storm_intensity(mslp):
    """ 
    Returns hurricane category based on Saffir-Simpson Hurricane 
    Wind Scale. Non-hurricanes return '--' 
    
    """
#    print 'mslp',mslp
    if (float(mslp) >= 994.):
        category=0 # Category 1
    elif (float(mslp) >= 980. and float(mslp) < 994.):
        category=1 # Category 2
    elif (float(mslp) >= 965. and float(mslp) < 980.):
        category=2 # Category 3  
    elif (float(mslp) >= 945. and float(mslp) < 965.):
        category=3 # Category 4       
    elif (float(mslp) >= 920. and float(mslp) < 945.):
        category=4 # Category 4       
    elif (float(mslp) >= 860. and float(mslp) < 920):
        category=5 # Category 5                   
    else:
        category=0
    return category
    
def storm_intensity_vmax(vmax):
    """ 
    Returns hurricane category based on Saffir-Simpson Hurricane 
    Wind Scale. Non-hurricanes return '--' 
    
    """
    if vmax >= 64 and vmax <= 82:
        category=1 # Category 1
    elif vmax >= 83 and vmax <= 95:
        category=2 # Category 2
    elif vmax >= 96 and vmax <= 112:
        category=3 # Category 3  
    elif vmax >= 113 and vmax <= 136:
        category=4 # Category 4       
    elif vmax >= 137:
        category=5 # Category 5                   
    else:
        category=0
    return category

def smooth(x, window_len=11, window='hanning'):
    """
    smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
   
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
           flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
   
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    #print('x ',x)
    xmax = np.amax(x)
    xmin = np.amin(x)
    #print('smooth, max/min of input ',xmax, xmin)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    
    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
#print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    #y[0:3] = np.ma.masked
    #y[-4:] = np.ma.masked
    #y[miss] = np.ma.masked
    y_to_return = y[(int(window_len/2)-1):-(int(window_len/2)+1)]
    miss = np.where((y_to_return < xmin) | (y_to_return > xmax))[0]
    #print('length ',len(x), len(y_to_return))
    #print('smooth, miss ',miss)
    #print('smoothed data, ',y_to_return)
    y_masked = np.ma.MaskedArray(y_to_return, mask=False)
    y_masked[miss] = np.ma.masked
    #print('masked data, ',y_masked)
    
#    return y
#    return y_to_return
    return y_masked

def storm_density_calculation_obs(storms, years, months, basin):
    """
    Calculate track density for observed storms
    
    """
#    storms = ts_obs.load_data.load_data()   
    lats, lons, count = ts_obs.example_code.storm_track_density.storm_lats_lons(storms, years, months)
    cube = ts_obs.example_code.storm_track_density._binned_cube(lats, lons)
    cube /= len(years)
    cube /= len(months)
    print('Total number of storms in time period:', count)
    return cube,count
    
def get_storms_in_year(storms, year, months):
    storms_in_year = []
    for storm in _storms_in_time_range(storms, year, months):
        storms_in_year.append(storm)
    return storms_in_year

def _get_time_range(year, months, calendar = 'proleptic_gregorian'):
    """ 
    Creates a start and end date (a datetime.date timestamp) for a 
    given year and a list of months. If the list of months overlaps into 
    the following year (for example [11,12,1,2,3,4]) then the end date 
    adds 1 to the original year 
    
    """
    time_unit = 'hours since 1950-01-01 00:00:00'
    start_date = cftime.datetime(year, months[0], 1, calendar = calendar)
    t = cftime.date2num(start_date, time_unit, calendar = calendar)
    t1 = cftime.num2date(t, time_unit, calendar = calendar)
    
    end_year = year
    end_month = months[-1]+1
    if months[-1]+1 < months[0] or months[-1]+1 == 13 or len(months) >= 12:
        end_year = year+1
    if months[-1]+1 == 13:
        end_month = 1
    end_date = cftime.datetime(end_year, end_month, 1, calendar = calendar)
    t = cftime.date2num(end_date, time_unit, calendar = calendar)
    t2 = cftime.num2date(t, time_unit, calendar = calendar)
    return t1, t2

def _storm_genesis_in_basin(storm, basin):
    """ 
    Returns True if the storm formed (i.e. the first lat/lon point) 
    occurred within the desired ocean basin. 
    
    """
    rbox = _basin_polygon(basin)  
    xy = ccrs.PlateCarree().transform_point(storm.obs_at_genesis().lon, 
                                            storm.obs_at_genesis().lat, 
                                            ccrs.Geodetic())
    point = sgeom.Point(xy[0], xy[1])
    if point.within(rbox):
        return True
    return False           

def _cube_data(data):
    """Returns a cube given a list of lat lon information."""
    cube = iris.cube.Cube(data)
    lat_lon_coord_system = icoord_systems.GeogCS(6371229)
    
    step = 4.0
    start = step/2
    count = 90
    pts = start + np.arange(count, dtype=np.float32) * step
    lon_coord = icoords.DimCoord(pts, standard_name='longitude', units='degrees', 
                                 coord_system = lat_lon_coord_system, circular=True)
    lon_coord.guess_bounds()
    
    start = -90
    step = 4.0
    count = 45
    pts = start + np.arange(count, dtype=np.float32) * step
    lat_coord = icoords.DimCoord(pts, standard_name='latitude', units='degrees', 
                                 coord_system = lat_lon_coord_system)
    lat_coord.guess_bounds()
    
    cube.add_dim_coord(lat_coord, 0)
    cube.add_dim_coord(lon_coord, 1)
    return cube 

def _binned_cube(lats, lons):
    """ Returns a cube (or 2D histogram) of lat/lons locations. """   
    data = np.zeros(shape=(45,90))
    binned_cube = _cube_data(data)
    xs, ys = binned_cube.coord('longitude').contiguous_bounds(), binned_cube.coord('latitude').contiguous_bounds()
    binned_data, _, _ = np.histogram2d(lons, lats, bins=[xs, ys])
    binned_cube.data = np.transpose(binned_data)
    binned_cube.attributes.pop('history', None) 
    binned_cube.units = cf_units.Unit(1)
    return binned_cube

def _get_time_period(years, months):
    """ 
    Returns string of time period for a given set of 
    years and months. E.g. months [6,7,8,9] and years 
    numpy.arange(1989,2003) would return a string 
    'June-September 1989-2002'. Note: years and 
    months must be lists or arrays.
    
    
    """    
    start_mon = calendar.month_name[months[0]]
    end_mon = calendar.month_name[months[::-1][0]]
    start_yr = str(years.min())
    end_yr = str(years.max())
    if start_yr == end_yr:
        return '%s-%s %s' % (start_mon, end_mon, start_yr)
    else:
        return '%s-%s %s-%s' % (start_mon, end_mon, start_yr, end_yr)

  
def storm_lats_lons(storms, years, months, basin, genesis=False, 
                 lysis=False, max_intensity=False):
    """ 
    Returns array of latitude and longitude values for all storms that 
    occurred within a desired year, month set and basin. 
    
    To get genesis, lysis or max intensity results set:
    Genesis plot: genesis=True
    Lysis plot: lysis=True
    Maximum intensity (location of max wind): 
    max_intensity=True
    
    """
    lats, lons = [], []
    count = 0
    for storm in storms:
        if (storm_in_time_range(storm, years, months) and 
        _storm_genesis_in_basin(storm, basin)):        
            if genesis:
                lats.extend([storm.obs_at_genesis().lat])
                lons.extend([storm.obs_at_genesis().lon])
            elif lysis:
                lats.extend([storm.obs_at_lysis().lat])
                lons.extend([storm.obs_at_lysis().lon])
            elif max_intensity:
                lats.extend([storm.obs_at_vmax().lat])
                lons.extend([storm.obs_at_vmax().lon])
            else:
                lats.extend([ob.lat for ob in storm.obs])
                lons.extend([ob.lon for ob in storm.obs])
            count += 1
                
    # Normalise lon values into the range 0-360
    norm_lons = []
    for lon in lons:
        norm_lons.append((lon + 720) % 360)
    return lats, norm_lons, count

def smooth(x,window_len=11,window='hanning'):
    """
    smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
   
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
           flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
   
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    #print('x ',x)
    xmax = np.amax(x)
    xmin = np.amin(x)
    #print('smooth, max/min of input ',xmax, xmin)
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")
    
    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")
    
    if window_len<3:
        return x
    
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    
    
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
#print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    #y[0:3] = np.ma.masked
    #y[-4:] = np.ma.masked
    #y[miss] = np.ma.masked
    y_to_return = y[(int(window_len/2)-1):-(int(window_len/2)+1)]
    miss = np.where((y_to_return < xmin) | (y_to_return > xmax))[0]
    #print('length ',len(x), len(y_to_return))
    #print('smooth, miss ',miss)
    #print('smoothed data, ',y_to_return)
    y_masked = np.ma.MaskedArray(y_to_return, mask=False)
    y_masked[miss] = np.ma.masked
    #print('masked data, ',y_masked)
    
#    return y
#    return y_to_return
    return y_masked
