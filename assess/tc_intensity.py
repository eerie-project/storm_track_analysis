""" 
Plots monthly variability of model tropical storms with
corresponding observations. Produces one plot with 
results for *all* individual ocean basins.

Note: this routine requires the 'ts_obs' module in order 
to work. This is a seperate Python module that needs to
be imported from fcm.


"""
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import numpy
import numpy as np
#import multipolyfit.multipolyfit as mpf
import os, glob, sys
import pylab as P
import math, copy
import random
import pickle
#import colormaps as cmaps
#plt.register_cmap(name='viridis', cmap=cmaps.viridis)

TC_SCRIPT_PATH='/data/users/hadom/branches/git/eerie-project/storm_track_analysis/assess'
sys.path.append(TC_SCRIPT_PATH)
import ts_obs.example_code
import ts_obs.load_data
import tc_assess_code

#st_assess = '/data/users/hadom/tracking'
#sys.path.append(st_assess)
#sys.path.append('/home/h06/hadom/workspace/PRIMAVERA-H2020_git_repository/HighResMIP-storm_reader')
#import storm_assess.track as st_track
#import storm_assess.functions as st_func

#sys.path.append('/data/users/hadom/tracking')
#import ts_obs.example_code
#import ts_obs.load_data
#import ts_model.track as ts_track
#import ts_model.example_code
#import ts_autoassess

symbols = {'glm':'v', 'nd_3dsm':'<', 'eg_sabl':'>', 'eg_ga6_entx2_cape15':'^', \
           'eg_ga6':'d', 'eg_ga6_entx1.5':'D'}
basins = ['na','ep','wp','ni'] #,'si','au','sp']
#    basins = ['na','ep']
basin_name=['N Atl','E Pac','W Pac','N Ind']
basin_ymax = [16.0, 30.0, 30.0, 10.0, 10.0, 20.0, 20.0]
SYM_SIZE = 4

colournames = {}
colournames['notpaired'] = ['red','blue','green','plum','aquamarine','coral','blueviolet','indianred','lime','lightgoldenrodyellow','pink','lightblue','plum']*3+['black']
colournames['paired'] = ['red','red','blue','blue','green','green','plum','plum','aquamarine','aquamarine','coral','coral','lime','lime','pink','pink','grey','darkgrey','slategrey','lightgrey']*3+['black']
colournames['paired'] = ['red','red','blue','blue','green','green','plum','plum','aquamarine','aquamarine','coral','coral','pink', 'pink', 'lime','lime','pink','grey','darkgrey','slategrey','lightgrey']*3+['black']
alpha = [0.4]*20
styles = {}
styles['notpaired'] = ['-','-','-','-','-','-','-','-','-','-', '-']
styles['paired'] = [':','-']*8 + ['-','-','-', '-']
dashList = [(2,1), (5,0)]*7 + [(5,0), (5,0), (5,0), (5,0)]
markers = ['v','x','*','s','^','o','1','v','x','*','s','o','v','x','*','s','o','v','x','*','s']

hurr_cat = ['1','2','3','4','5']
knots_to_ms = 1.0 / 1.944
mslp_limits = [994,980,965,945,920,880]
wind_limits = [33,43,49,58,70,80]
wind_max = 80
mslp_range = np.arange(890, 1025, 5)
life_range = np.arange(0, 26, 2)
lat_range = np.arange(2, 60, 2)
vmax_range = np.arange(0, 100, 5)
core_range = np.arange(-20, 40, 3)
vort_range = np.arange(2, 30, 2)

convert_925_10mwind_factor = 0.76

def hires_storms():
    '''
    Some example storms from the 4km model
    '''
    ernesto = {}
    ernesto['glm'] = (994, 53)
    ernesto['nd_3dsm'] = (996, 56)
    ernesto['eg_sabl'] = (1000, 44)
    ernesto['eg_ga6_entx2_cape15'] = (954, 88)
    ernesto['eg_ga6'] = (1007, 34)
    ernesto['eg_ga6_entx1.5'] = (1003, 56)
    
    leslie = {}
    leslie['glm'] = (950, 77)
    leslie['nd_3dsm'] = (948, 92)
    leslie['eg_sabl'] = (919, 110)
    leslie['eg_ga6_entx2_cape15'] = (918, 97)
    leslie['eg_ga6'] = (937, 82)
    leslie['eg_ga6_entx1.5'] = (951, 74)
    
    bopha = {}
    bopha['glm'] = (976, 62)
    bopha['nd_3dsm'] = (956, 82)
    bopha['eg_sabl'] = (936, 90)
    bopha['eg_ga6_entx2_cape15'] = (911, 109)
    bopha['eg_ga6'] = (945, 90)
    bopha['eg_ga6_entx1.5'] = (916, 109)
    
    bolaven = {}
    bolaven['glm'] = (940, 73)
    bolaven['nd_3dsm'] = (949, 88)
    bolaven['eg_sabl'] = (920, 109)
    bolaven['eg_ga6_entx2_cape15'] = (896, 112)
    bolaven['eg_ga6'] = (939, 82)
    bolaven['eg_ga6_entx1.5'] = (893, 109)
    
    storms_hires = [ernesto, leslie, bopha, bolaven]
    
    return storms_hires
    

def _get_annual_storms(storms, years, months, basin):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    for storm in storms:
        if (storm.genesis_date().year in years) and tc_assess_code._storm_vmax_in_basin(storm, basin):
#            print storm.genesis_date().year
            yield storm

def _get_annual_storms_obs(storms, years, months, basin, storm_type):
    """ 
    Returns array of storm counts for each year for a given set of months 
    and storm types. Default storm type is to count all tropical cyclones 
    (1 min winds > 33 kts) 
    
    """
    for storm in storms:        
#        if (storm.genesis_date().year in years) and _storm_in_basin(storm, basin) and _storm_in_strength_range(storm, storm_type):
        if (storm.genesis_date().year in years) and tc_assess_code._storm_vmax_in_basin(storm, basin) and \
                ts_obs.example_code._storm_in_strength_range(storm, storm_type):
#            print storm.genesis_date().year
            yield storm

def load_example_storms(fnames):
    """ 
    Returns list of example model tropical storms
    for the northern and southern hemisphere. 
    
    Note: Northern Hemsiphere file covers the period
    May-November 2000-2011; the Southern Hemisphere
    file the period October-April 2000-2011.  
    
    """
    if 'list' not in str(type(fnames)):
        fnames = [fnames]
    storms = []
    for fname in fnames:
#        storms.extend(list(ts_model.track.load(fname, ex_cols=3, 
#                                               calendar='netcdftime'))
#                      )    
        storms.extend(list(ts_track.load(fname, ex_cols=3, 
                                               calendar='netcdftime'))
                      )    
    return storms
    
    
def lineplot_setup(xaxis_range,yaxis_range,title,xtitle,ytitle):
    """ 
    Set up a lineplot 
    
    """
    fig.subplots_adjust(right=0.83, left=0.17, bottom=0.12)
    plt.grid()
    plt.title('%s' % title)  
    ax = plt.gca()
    ax.set_xlabel(xtitle, fontsize='large')
    ax.set_ylabel(ytitle, ha='center', color='black', fontsize='large')
    ax.set_xlim(xaxis_range[0]-1, xaxis_range[len(xaxis_range)-1]+1)
    ax.set_ylim(yaxis_range[0],yaxis_range[1])
    

def setup_Reanalysis_barchart_data(years,months,basins):
    """ 
    Set up all the standard annual counts for reanalyses for given set of years 
    
    """
    base_dir = '/data/local/hadom/TRACK/TCTRACK_DATA_vn1.40/'
    base_filename = 'combined_ff_trs.vor_10m_fullgrid_N512_xgxqi_L5.new.wc_19852011.date'
    base_name = 'combined_ff_trs.vor_10m_fullgrid_'
    directory = '/data/local/hadom/TRACK/TCTRACK_DATA_vn1.40/xgxqi/TOTAL/'

# standard filename years
    year_start = 1985
    year_end = 2011

    runid_erai = 'erai'; resol_erai = 'ERAI'
    runid_jra25 = 'jra25'; resol_jra25 = 'JRA25'
    runid_merra = 'merra'; resol_merra = 'MERRA'
    print(years)

    fnames_erai = []
    fnames_erai.append(base_dir + runid_erai+'/TOTAL/'+base_name + resol_erai + '_' + runid_erai + '_L5.new.wc_' + str(year_start) + str(year_end)+ '.date')
    fnames_jra25 = []
    fnames_jra25.append(base_dir + runid_jra25+'/TOTAL/'+base_name + resol_jra25 + '_' + runid_jra25 + '_L5.new.wc_' + str(year_start) + str(year_end)+ '.date')
    fnames_merra = []
    fnames_merra.append(base_dir + runid_merra+'/TOTAL/'+base_name + resol_merra + '_' + runid_merra + '_L5.new.wc_' + str(year_start) + str(year_end)+ '.date')
    print(fnames_erai,fnames_jra25,fnames_merra)
# read in reanalysis st
    storms_erai = load_example_storms(fnames_erai)
    storms_jra25 = load_example_storms(fnames_jra25)
    storms_merra = load_example_storms(fnames_merra)
   
    obs_storms = ts_obs.load_data.load_data()

# set up arrays for mean counts
    annual_count_erai=numpy.zeros(len(basins)*len(years)).reshape(len(basins),len(years))
    annual_count_jra25=numpy.zeros(len(basins)*len(years)).reshape(len(basins),len(years))
    annual_count_merra=numpy.zeros(len(basins)*len(years)).reshape(len(basins),len(years))
    annual_count_hurdath=numpy.zeros(len(basins)*len(years)).reshape(len(basins),len(years))
    annual_count_hurdata=numpy.zeros(len(basins)*len(years)).reshape(len(basins),len(years))

    print(len(years))
    iyr=0
# calculate the number of storms in each basin for each model
    for i, basin in enumerate(basins):        
        nensemble=1.0
        annual_count_erai[i]  = tc_assess_code.get_annual_vmax_mean_count(storms_erai, years, months, basin, nensemble)
        annual_count_jra25[i] = tc_assess_code.get_annual_vmax_mean_count(storms_jra25, years, months, basin, nensemble)
        annual_count_merra[i] = tc_assess_code.get_annual_vmax_mean_count(storms_merra, years, months, basin, nensemble)
        annual_count_hurdath[i] = ts_obs.example_code.annual_storm_counts(obs_storms, years, months, basin, storm_types=['HU', 'MH'])
        annual_count_hurdata[i] = ts_obs.example_code.annual_storm_counts(obs_storms, years, months, basin, storm_types=['TS', 'HU', 'MH'])

    return annual_count_erai,annual_count_jra25,annual_count_merra,annual_count_hurdath,annual_count_hurdata

def sample_data():
    subsample = 10
# do a random sampling of the data
    data_size = len(x)
    idx = random.sample(list(range(data_size)), min(2000, data_size))
#        ax.scatter(x[::subsample], y[::subsample],s=10,color='grey',alpha=alpha[i],marker=markers[i])
#        ax.scatter(x[::subsample], y[::subsample],s=10,color='grey',alpha=alpha[i],marker=markers[i])
    print(idx)
    float_x = np.zeros(len(x))
    float_y = np.zeros(len(x))
    for ind in np.arange(len(x)):
        float_x[ind] = float(x[ind])
        float_y[ind] = float(y[ind])
# random subsampling
#        ax.scatter(float_x[idx], float_y[idx],s=10,color='grey',alpha=alpha[i],\
#                   marker=markers[i], rasterized = True)
    pass

def calculate_fit_obs(storms, basin = 'nh'):
    x, y = [], []
    knapp_vmax = []; knapp_mslp = []
    vmax_10min = []
    for storm in storms:
        if (storm.obs_at_vmax().mslp != -999 and tc_assess_code._storm_vmax_in_basin(storm, basin) ):
            x.append(storm.obs_at_vmax().vmax * knots_to_ms)
            y.append(storm.obs_at_vmax().mslp)
            
            # knapp model fit
            delta_p = storm.obs_at_vmax().mslp - 1014.3
            kvmax = (18.633 - 14.96*1. - 0.755*storm.obs_at_vmax().lat - 0.518*delta_p \
                              + 9.738 * (math.fabs(delta_p))**0.5 + 1.5*(9.6**0.63))/1.944
            knapp_vmax.append(kvmax)
            knapp_mslp.append(storm.obs_at_vmax().mslp)
            vmax_10min.append(0.88*storm.obs_at_vmax().vmax * knots_to_ms)

    #print('hurdat x, ',x)
    #print('hurdat y, ',y)
    coeffs = numpy.polyfit(x,y,2)
    coeffs_10min = numpy.polyfit(vmax_10min,y,2)
    x2 = numpy.arange(min(x)-1, max(x)+1, .01)
    y2 = numpy.polyval(coeffs, x2)
    y_10min = numpy.polyval(coeffs_10min, x2)
    #y_10min = x2

    return x2, y_10min

def calculate_fit_model(storms, use_925 = False, basin = 'nh'):
    '''
    calculate 10mwind max, corresponding mslp min, and calculate fitcurve
    if use_925 is true, use 925 winds instead (scaled by 0.76)
    '''
    x, x1, y = [], [], []
    for storm in storms:
        if tc_assess_code._storm_vmax_in_basin(storm, basin):
            if storm.obs_at_vmax().extras['w10m'] < 1.0e10 and storm.obs_at_vmax().mslp > 600 and np.abs(storm.obs_at_vmax().lat) < 40.0 and np.abs(storm.obs_at_vmax().vmax) < 1.0e10:
                w10m = storm.obs_at_vmax().extras['w10m']
                w925 = storm.obs_at_vmax().vmax
                mslp = storm.obs_at_vmax().mslp
                if mslp > 1500. or mslp < 600:
                    continue
                try:
                    if np.isfinite(w10m) and np.isfinite(w925) and np.isfinite(mslp):
                        x.append(storm.obs_at_vmax().extras['w10m'])
                        x1.append(storm.obs_at_vmax().vmax)
                        y.append(storm.obs_at_vmax().mslp)
                except:
                    pass

    if len(x) > 0:
        #print('max 10m wspeed ',np.max(x))
        #print('max 925 wspeed ',np.max(x1))
        #print('min mslp ',np.min(y))
        x = numpy.array(x)
        x1 = numpy.array(x1)
        x1 = x1 * convert_925_10mwind_factor
        coeffs = numpy.polyfit(x, y, 2)
        coeffs1 = numpy.polyfit(x1, y, 2)
        if not use_925:
            x2 = numpy.arange(min(x)-1, max(x)+1, .01)
            y_fit = numpy.polyval(coeffs, x2)
        else:
            x2 = numpy.arange(min(x1)-1, max(x1)+1, .01)
            y_fit = numpy.polyval(coeffs1, x2)
            
    else:
        x2 = numpy.arange(0,50,5)
        y_fit = np.zeros((len(x2)))
    
    return x2, y_fit

def max_storm_per_year(storms, years):
    '''
    Find the strongest storm (min lat) per year
    '''
    storm_min = {}
    storm_min['mslp'] = []
    storm_min['year'] = []
    for year in years:
        storm_min['year'].append(year)
        min_mslp = 1200.
        for storm in storms:
            if storm.obs_at_max().year == year:
                if storm.obs_at_vmax().mslp < min_mslp:
                    min_mslp = storm.obs_at_vmax().mslp
        storm_min['mslp'] = min_mslp

    return storm_min

def calculate_histogram_latitudes(storms):
    lat_storm = []
    for storm in storms:
#        if (storm.obs_at_vmax().mslp != -999):
        lat_storm.append(storm.obs_at_vmax().lat)
    # histogram for storm latitudes
    hist, bin_edges = numpy.histogram(lat_storm, bins = lat_range, density = True)
    return hist, bin_edges
    

def calculate_histogram_lifetime(storms):
    lifetime_storm = []
    for storm in storms:
        length = len(storm.obs)/4.0
        lifetime_storm.append(length)
    # histogram for storm latitudes
    hist, bin_edges = numpy.histogram(lifetime_storm, bins = numpy.arange(14), density = True)
    return hist, bin_edges
    
def calculate_histogram_maxmslp(storms):
    maximum_storm = []
    for storm in storms:
        if storm.obs_at_vmax().mslp > 500. and storm.obs_at_vmax().mslp < 2000.:
            max_strength = storm.obs_at_vmax().mslp
            maximum_storm.append(max_strength)
    # histogram for storm latitudes
    hist, bin_edges = numpy.histogram(maximum_storm, bins = mslp_range, range = (mslp_range[0], mslp_range[-1]), density = True)
    return hist, bin_edges
    
def calculate_histogram_maxvmax(storms, scaling = 1.0, wind10m=False, basin = 'nh'):
    maximum_storm = []
    for storm in storms:
        if tc_assess_code._storm_vmax_in_basin(storm, basin):
            if not wind10m:
                max_strength = storm.obs_at_vmax().vmax*scaling
            else:
                max_strength = storm.obs_at_vmax().extras['w10m']
            maximum_storm.append(max_strength)
    # histogram for storm latitudes
    hist, bin_edges = numpy.histogram(maximum_storm, bins = vmax_range, range = (vmax_range[0], vmax_range[-1]), density = True)
    return hist, bin_edges
    
def calculate_histogram_maxlife(storms):
    maximum_storm = []
    lifetime_storm = []
    for storm in storms:
        if storm.obs_at_vmax().lat <= 60.0 and storm.obs_at_vmax().mslp < 2000:
            max_strength = storm.obs_at_vmax().mslp
            if max_strength > 1500:
                print(storm.genesis_date())
                print(storm.obs_at_vmax().lat)
            maximum_storm.append(max_strength)
            length = len(storm.obs)/4.0
            lifetime_storm.append(length)
    print('min life ',np.min(lifetime_storm))
    print('max mslp ',np.max(maximum_storm))
    # histogram for storm latitudes
    xedges = mslp_range
    yedges = life_range
    hist, xedges, yedges = numpy.histogram2d(maximum_storm, lifetime_storm, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def calculate_histogram_maxlat(storms, basin = 'nh'):
    maximum_storm = []
    lat_storm = []
    for storm in storms:
        if tc_assess_code._storm_vmax_in_basin(storm, basin):
            if storm.obs_at_vmax().lat <= 60.0 and storm.obs_at_vmax().mslp < 2000:
                max_strength = storm.obs_at_vmax().mslp
                maximum_storm.append(max_strength)
                lat_storm.append(storm.obs_at_vmax().lat)
    # histogram for storm latitudes
    xedges = mslp_range
    yedges = lat_range
    hist, xedges, yedges = numpy.histogram2d(maximum_storm, lat_storm, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    

def calculate_histogram_warmcorevmax(storms, runid):
    maximum_core = []
    vmax_storm = []
    for storm in storms:
        if np.absolute(storm.obs_at_vmax().lat) <= 60.0:
            max_strength = storm.obs_at_vmax().vmax
            vmax_storm.append(max_strength)
            if runid == 'nicam':
                max_core = (storm.obs_at_vmax().t63_2 - storm.obs_at_vmax().t63_5)
            else:
                max_core = (storm.obs_at_vmax().t63_diff)
            maximum_core.append(max_core)
    # histogram for storm latitudes
    xedges = vmax_range
    yedges = core_range
    hist, xedges, yedges = numpy.histogram2d(vmax_storm, maximum_core, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def calculate_histogram_v850_vgrad(storms, runid):
    print('calc_vortmax, runid ',runid)
    maximum_vgrad = []
    vort_storm = []
    for storm in storms:
        if np.absolute(storm.obs_at_vmax().lat) <= 60.0:
            #max_strength = storm.obs_at_vmax().vort_max
            #print 'storm str ',storm.obs_at_vmax()
            max_strength = storm.obs_at_vmax().vort
            vort_storm.append(max_strength)
#            max_strength = storm.obs_at_vmax().vort
#            vort_storm.append(max_strength)
            if runid == 'nicam':
                    #max_core = (storm.obs_at_vmax().t63_2 - storm.obs_at_vmax().t63_5)
                max_core = storm.obs_at_vmax().t63_diff
            else:
#max_core = (storm.obs_at_vmax().t63_diff)
                max_core = storm.obs_at_vmax().t63_diff
            maximum_vgrad.append(max_core)
    # histogram for storm latitudes
    xedges = vort_range
    yedges = core_range
    hist, xedges, yedges = numpy.histogram2d(vort_storm, maximum_vgrad, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def calculate_histogram_vmaxmslp(storms, runid):
    min_mslp = []
    vmax_storm = []
    for storm in storms:
        if np.absolute(storm.obs_at_vmax().lat) <= 60.0 and storm.obs_at_vmax().mslp > 600 and storm.obs_at_vmax().mslp < 2000 :
            max_strength = storm.obs_at_vmax().vmax
            if runid == 'hurdat':
                max_strength*= 0.88 * knots_to_ms

            vmax_storm.append(max_strength)
            min_mslp_storm = storm.obs_at_vmax().mslp
            min_mslp.append(min_mslp_storm)
    # histogram for storm latitudes
    xedges = vmax_range
    yedges = mslp_range
    hist, xedges, yedges = numpy.histogram2d(vmax_storm, min_mslp, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def calculate_histogram_lifelat(storms):
    lifetime_storm = []
    lat_storm = []
    for storm in storms:
        length = len(storm.obs)/4.0
        lifetime_storm.append(length)
        lat_storm.append(storm.obs_at_vmax().lat)
    # histogram for storm latitudes
    xedges = life_range
    yedges = lat_range
    hist, xedges, yedges = numpy.histogram2d(lifetime_storm,lat_storm, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def calculate_histogram_lifegenlat(storms):
    lifetime_storm = []
    genlat_storm = []
    for storm in storms:
        length = len(storm.obs)/4.0
        lifetime_storm.append(length)
        genlat_storm.append(storm.obs_at_genesis().lat)
    # histogram for storm latitudes
    xedges = life_range
    yedges = lat_range
    hist, xedges, yedges = numpy.histogram2d(lifetime_storm,genlat_storm, bins = (xedges, yedges), density = True)
    hist = hist.T
    return hist, xedges, yedges
    
def plot_hart_lifecycle(storms):
    fig = plt.figure(figsize=(10,6),dpi=100)
    ax = fig.add_subplot(111)
    for storm in storms:
        xphase = []; yphase = []
        for ob in storm.obs:
            v850_m_v600 = -(ob.t63_1 - ob.t63_3)
            v600_m_v250 = -(ob.t63_3 - ob.t63_5)
            if v850_m_v600 > 0. and v600_m_v250 > 0.:
                xphase.append(v850_m_v600)
                yphase.append(v600_m_v250)
        plt.plot([xp for xp in xphase], [yp for yp in yphase],
                     linewidth=1.2)
    print('no storms ',len(xphase))
    ax.set_xlim(-25,25)
    ax.set_ylim(-25,25)
    #plt.show()
        

def plot_obs_mslp_wspeed(storms, ax, x_fit=[], y_fit = [], basin = 'nh', linewidth = 2.5, plot_scatter = True):
    for storm in storms:
#        ax1.scatter(storm.obs_at_vmax().vmax*knots_to_ms, storm.obs_at_vmax().lat,s=20, \
#                    color='black', rasterized = True)

        #if ts_obs.example_code._storm_in_basin(storm, basin) and (storm.obs_at_vmax().mslp != -999):
        if tc_assess_code._storm_vmax_in_basin(storm, basin) and (storm.obs_at_vmax().mslp != -999):
            ax.scatter(0.88*storm.obs_at_vmax().vmax*knots_to_ms, storm.obs_at_vmax().mslp, s=SYM_SIZE, \
                       color='black', rasterized = True)
    
    if len(x_fit) != 0:
        ax.plot(x_fit, y_fit, '-' , color='black', linewidth=linewidth, label='OBS (10min)', \
                rasterized = True)

def plot_model_mslp_wspeed(storms, i, ax, x_fit=[], y_fit =[], label = '', use_925 = False, paired = False, basin = 'nh', linewidth = 2, plot_scatter = True):
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if plot_scatter:
        for storm in storms:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                if not use_925:
                    wspeed = storm.obs_at_vmax().extras['w10m']
                    if wspeed > 20:
                        ax.scatter(storm.obs_at_vmax().extras['w10m'],storm.obs_at_vmax().mslp, s=SYM_SIZE/2,\
                                   color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
                else:
                    ax.scatter(storm.obs_at_vmax().vmax*convert_925_10mwind_factor,storm.obs_at_vmax().mslp, s=SYM_SIZE/2,\
                                   color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
    
    if len(x_fit) != 0:
        if label != '':
            if paired:
                ax.plot(x_fit, y_fit, '-' , color=colourname[i], dashes = dashList[i], linewidth=linewidth, label=label, rasterized = True)
            else:
                ax.plot(x_fit, y_fit, '-' , color=colourname[i], linewidth=linewidth, label=label, rasterized = True)
        else:
            if paired:
                ax.plot(x_fit, y_fit, '-' , color=colourname[i], dashes = dashList[i], linewidth=linewidth, rasterized = True)
            else:
                ax.plot(x_fit, y_fit, '-' , color=colourname[i], linewidth=linewidth, rasterized = True)

def plot_model_mslp_wspeed_noscatter(storms, i, ax, x_fit=[], y_fit =[], label = '', use_925 = False, paired = False, basin = 'nh', plot_scatter = False, linewidth = 2):
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if plot_scatter:
        for storm in storms:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                if not use_925:
                    ax.scatter(storm.obs_at_vmax().extras['w10m'],storm.obs_at_vmax().mslp, s=SYM_SIZE,\
                color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
                else:
                    ax.scatter(storm.obs_at_vmax().vmax*convert_925_10mwind_factor,storm.obs_at_vmax().mslp, s=SYM_SIZE,\
                color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
    
    if len(x_fit) != 0:
        if label != '':
            ax.plot(x_fit, y_fit, '-' , color=colourname[i], dashes = dashList[i], linewidth=linewidth, label=label, rasterized = True)
        else:
            ax.plot(x_fit, y_fit, '-' , color=colourname[i], dashes = dashList[i], linewidth=linewidth, rasterized = True)

def plot_model_mslp_size(storms, i, ax, x_fit=[], y_fit =[], label = '', use_925 = False, paired = False, basin = 'nh', linewidth = 2, plot_scatter = True):
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if plot_scatter:
        for storm in storms:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                wspeed = storm.obs_at_vmax().extras['w10m']
                mslp = storm.obs_at_vmax().mslp
                try:
                    size = storm.obs_at_vmax().radius_8
                    ax.scatter(mslp, size, s=SYM_SIZE/2,\
                                   color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
                    if size > 6:
                        print('wspeed, mlsp, size ', wspeed, mslp, size, storm.genesis_date(), storm.obs_at_vmax().lat, storm.obs_at_vmax().lon, storm.obs_at_genesis().lat, storm.obs_at_genesis().lon, storm.warm_core)
                except:
                    continue

def plot_model_warm_core(storms, i, ax, x_fit=[], y_fit =[], label = '', use_925 = False, paired = False, basin = 'nh', linewidth = 2, plot_scatter = True):
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if plot_scatter:
        for storm in storms:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                for ob in storm.obs:
                    delta_z_250 = ob.zg_max_250 - ob.zg_min_250
                    delta_z_500 = ob.zg_max_500 - ob.zg_min_500
                    z_diff.append(delta_z_500 - delta_z_250)
                    z_max = np.amax(z_diff)
                try:
                    size = storm.obs_at_vmax().radius_8
                    ax.scatter(mslp, size, s=SYM_SIZE/2,\
                                   color=colourname[i], alpha=alpha[i], marker=markers[i], rasterized = True)
                    if size > 6:
                        print('wspeed, mlsp, size ', wspeed, mslp, size, storm.genesis_date(), storm.obs_at_vmax().lat, storm.obs_at_vmax().lon, storm.obs_at_genesis().lat, storm.obs_at_genesis().lon, storm.warm_core)
                except:
                    continue

def plot_latitude_storm(storms):
    fig = plt.figure(figsize=(10,6),dpi=100)
    ax = fig.add_subplot(111)
    for ist, storm in enumerate(storms):
        ax.scatter(ist, storm.obs_at_vmax().lat)
    #plt.show()

def plot_latitude_histogram(histogram, runid, resols, model_desc, plot_dir, edges = [], x_label = 'Latitude', title='', paired = False, plot_cat = True, savefname = ''):
    fig = plt.figure(figsize=(10,6),dpi=100)
    ax = fig.add_subplot(111)
    if edges == []:
        bins = numpy.arange(49)
    else:
        bins = edges[:-1]
    #print 'bins ',bins
    offset = (bins[1] - bins[0])/2.0
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if 'hurdat' in histogram:
        plt.plot(bins+offset,histogram['hurdat'][:], color='black', label='OBS', lw =3)

    for i, run in enumerate(runid):
        resol = resols[i]
        key = run+resol
        if paired:
            if i%2 == 1 or i > 13:
                label = model_desc[i]
            else:
                label = ''
        else:
            label = model_desc[i]
            
        plt.plot(bins+offset,histogram[key][:], color=colourname[i], dashes = dashList[i], label=label, lw =3)

    ax.set_ylabel('Normalised frequency')
    ax.legend(loc="center right", fontsize='medium', fancybox = True, framealpha = 1.0)
    plt.title(title, loc = 'left')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0.0, ymax)
    if 'windspeed' in x_label:
        ax.set_xlim(0.0, 85.0)
        ax.set_xlabel(x_label+' (ms$^{-1}$)')
    elif 'Mslp' in title:
        ax.set_xlim(880.0, 1020.0)
        ax.set_xlabel(x_label+' (hPa)')

    if 'windspeed' in x_label and plot_cat:
        # plot the wind speed category marks
        ax.plot([wind_limits[0],wind_limits[0]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[1],wind_limits[1]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[2],wind_limits[2]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[3],wind_limits[3]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[4],wind_limits[4]],[0.0, ymax], color='black', lw=1)
        #ax.plot([wind_limits[5],wind_limits[5]],[0.0, ymax], color='black', lw=1)
        for il in np.arange(0,5):
            if (wind_limits[il]+wind_limits[il+1])/2 <= wind_max:
                ax.text((wind_limits[il]+wind_limits[il+1])/2,ymax-ymax*0.05,hurr_cat[il])

    if savefname == '':
        savefname = title
    plot_filename = plot_dir+'/intensity_scatter_latitude_'+savefname+'.pdf'
    #plt.savefig(plot_filename)
    plt.savefig(plot_filename[:-3]+'png')
    #plt.show()    

def plot_latitude_histogram_nofig(histogram, runid, resols, model_desc, plot_dir, ax, edges = [], x_label = 'Latitude', title='', paired = False, plot_cat = True, savefname = '', fig_no = 1):
    if edges == []:
        bins = numpy.arange(49)
    else:
        bins = edges[:-1]
    #print 'bins ',bins
    offset = (bins[1] - bins[0])/2.0
    if paired:
        colourname = colournames['paired']
        style = styles['paired']
    else:
        colourname = colournames['notpaired']
        style = styles['notpaired']

    if 'hurdat' in histogram:
        plt.plot(bins+offset,histogram['hurdat'][:], color='black', label='OBS', lw =3)

    for i, run in enumerate(runid):
        resol = resols[i]
        key = run+resol
        if i%2 == 1 or i > 11:
            plt.plot(bins+offset,histogram[key][:], color=colourname[i], dashes = dashList[i], label=model_desc[i], lw =3)
        else:
            plt.plot(bins+offset,histogram[key][:], color=colourname[i], dashes = dashList[i], lw =3)
    ax.set_ylabel('Normalised frequency')
    ax.legend(loc="center right", fontsize='medium', fancybox = True, framealpha = 1.0)
    plt.title(title, loc = 'left')
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(0.0, ymax)
    if 'windspeed' in x_label:
        ax.set_xlim(0.0, 85.0)
        ax.set_xlabel(x_label+' (ms$^{-1}$)')

    if plot_cat:
        # plot the wind speed category marks
        ax.plot([wind_limits[0],wind_limits[0]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[1],wind_limits[1]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[2],wind_limits[2]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[3],wind_limits[3]],[0.0, ymax], color='black', lw=1)
        ax.plot([wind_limits[4],wind_limits[4]],[0.0, ymax], color='black', lw=1)
        #ax.plot([wind_limits[5],wind_limits[5]],[0.0, ymax], color='black', lw=1)
        for il in np.arange(0,5):
            if (wind_limits[il]+wind_limits[il+1])/2 <= wind_max:
                ax.text((wind_limits[il]+wind_limits[il+1])/2,ymax-ymax*0.05,hurr_cat[il])


def plot_2d_histogram(histogram, runid, resol, desc, plot_dir, fig, ax, nx, ny, xedges = [], yedges= [], x_label = 'Latitude', title='', xlabel = 'Min. MSLP', ylabel = 'Lifetime (days)', max_val = 0.005, min_val = 0.0005, show_title = True):

    x, y = np.meshgrid(xedges, yedges)
    print(runid)
    #print histogram

    ioff = 1
    iplt = 0
    cs = ax.pcolormesh(x, y, histogram, vmin=min_val, vmax=max_val)
    ax.set_title(desc)
    #ax.locator_params(axis='x', numticks=4)
    ax.set_xticks(np.arange(900, 1050, 50))

    #if not np.mod(i, nx) == 0:
    #    #ax.set_yticks([])
    #    ax.set_yticklabels([])
    #    ax.set_ylabel('')
    #else:
    #    ax.set_ylabel(ylabel)

    ax.set_xlabel(xlabel)

    ax.set_xlim(880, 1020)

    cs.cmap.set_under(color='white')
    #ax.set_xticks(np.arange(900, 1050, 50))
    #ax.set_xlim(880, 1020)
    #ax.set_yticks([])
    #ax.set_ylabel('')
    #ax.set_xlabel(xlabel)

    if show_title:
        plt.suptitle(title)
    
def setup_mslp_windspeed_plot(ax, basin = 'nh', fig_no = 1):
    ax.plot([wind_limits[0],wind_limits[0]],[880,1020], color='black')
    ax.plot([wind_limits[1],wind_limits[1]],[880,1020], color='black')
    ax.plot([wind_limits[2],wind_limits[2]],[880,1020], color='black')
    ax.plot([wind_limits[3],wind_limits[3]],[880,1020], color='black')
    ax.plot([wind_limits[4],wind_limits[4]],[880,1020], color='black')
    ax.plot([wind_limits[5],wind_limits[5]],[880,1020], color='black')

    for il in np.arange(0,5):
        if (wind_limits[il]+wind_limits[il+1])/2 <= wind_max:
            ax.text((wind_limits[il]+wind_limits[il+1])/2,1015,hurr_cat[il])
    for iran in range(len(mslp_limits)-1):
        ax.plot([0,wind_max],[mslp_limits[iran],mslp_limits[iran]], color='black', rasterized = True)
        y1 = (mslp_limits[iran]+mslp_limits[iran+1])/2.
#        ax.text(wind_max-2,y1,hurr_cat[iran])
        ax.text(wind_max,y1,hurr_cat[iran])
    ax.set_xlim(0., wind_max)
    ax.set_ylim(880.,1020)
#    ax.set_xlabel('Maximum lifetime 10 m wind speed '+r'$(m s^{-1})$')
    if fig_no == 1:
        ax.set_ylabel('Minimum MSLP (hPa)')
        ax.set_xlabel('')
        #ax.set_xticklabels([])
    elif fig_no == 2:
        ax.set_ylabel('')
        #ax.set_yticklabels([])
        ax.set_xlabel('Max. lifetime 10m wind speed (ms$^{-1}$)')
    elif fig_no == 3:
        ax.set_xlabel('Max. lifetime 10m wind speed (ms$^{-1}$)')
        ax.set_ylabel('Minimum MSLP (hPa)')
        pass

    ax2 = ax.twiny()
    ax3 = ax.twinx()
    mstoknots = (1.0 / knots_to_ms)
    ax2.set_xlim(0, wind_max*mstoknots)
    ax2.set_xticklabels([0, 20, 40, 60, 80, 100, 120, 140, 160])
    ax2.set_ylabel('MSLP category')
    #ax3.tick_params(axis='both', which=u'both',length=0)
    ax3.yaxis.set_tick_params(length = 0)
    if fig_no == 1:
        ax2.set_xlabel('Max. lifetime 10m wind speed (knots)')
        ax3.set_ylabel('')
        #ax3.set_xticklabels([])
        ax3.set_yticklabels([])
    elif fig_no == 2:
        ax2.set_xlabel('Max. lifetime 10m wind speed (knots)')
        ax3.set_ylabel('MSLP category')
        ax3.set_yticklabels([])
    elif fig_no == 3:
        ax2.set_xlabel('')
        ax3.set_ylabel('MSLP category')
        #ax3.set_xticklabels([])
        ax3.set_yticklabels([])

def get_model_storm_data(runid, resol, datadir, i, storm_type, years, months, basin, method, algo = ''):
    fname_nh, start_year, end_year = find_last_year(datadir, runid, resol, storm_type, method[i], algo)
    if runid == 'merra2':
        storms = list(ts_track.load(fname_nh, calendar='netcdftime', ex_cols=3))
    elif runid == 'nicam':
        storms = list(ts_track.load(fname_nh, calendar='netcdftime', ex_cols=0))

    else:
        storms = list(ts_track.load(fname_nh, calendar='netcdftime'))
#    storms = list(ts_model.track.load(fname_nh, calendar='netcdftime'))
    storms_yr = list(_get_annual_storms(storms, years, months, basin))
    if len(storms_yr) < 1:
        raise Exception('No storms found in '+fname_nh)
    return storms_yr

def find_last_year(datadir, runid, resol, track_type, method, algo = ''):
    if algo == '':
        files = glob.glob(os.path.join(datadir, runid+resol, 'combined_UVT_*'+track_type+'*.NH'+method+'*tracks'))
        print('search ',os.path.join(datadir, runid+resol, 'combined_UVT_*'+track_type+'*.NH'+method+'*tracks'))
    else:
        files = glob.glob(os.path.join(datadir, runid+resol, 'TempExt.combined_UVT_*.NH'+method+'*tracks'))
        print('search ',os.path.join(datadir, runid+resol, 'TempExt.combined_UVT_*.NH'+method+'*tracks'))
    files = sorted(files)
    print('load file ',files[-1])
    
    f_cmp = os.path.basename(files[-1]).split('_')
    start_year = f_cmp[2]
    end_year = f_cmp[3][:4]

    print(start_year, end_year)
    return files[-1], start_year, end_year

def work(directory, years, months, datadir, basin, track_type, run_dict, plot_dir, paired=False, algo = '', do_new = False, do_new_obs = False, T63_tracking = False):
    year_period = str(years[0])+'-'+str(years[-1])
    title = '10m wind vs MSLP '+str(years[0])+'-'+str(years[len(years)-1])
    title1 = '925hPa wind vs MSLP '+str(years[0])+'-'+str(years[len(years)-1])
    
    runid_model = run_dict['runid_model']
    method = run_dict['method']
    model_desc = run_dict['model_desc']
    model_resol = run_dict['resol_model']
    resols = run_dict['resol']

    storms_yr = {}

    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        filename = 'scatter_plot_highresmip_assessment_'+year_period+runid+resol+track_type+basin+algo+'.pkl'
        print('filename',os.path.join(directory,filename))
        if not os.path.isfile(os.path.join(directory,filename)) or do_new:

            print('making new file')
            storms = get_model_storm_data(runid, resol, datadir, i, track_type, years, months, basin, method, algo = algo)
            fh = open(os.path.join(directory, filename), 'wb') # write binary file mode
            pickle.dump(storms, fh)
            fh.close()
            storms_yr[key] = storms
        else:
            fh = open(os.path.join(directory, filename), 'rb') #read binary mode
            storms = pickle.load(fh)
            fh.close()
            storms_yr[key] = storms
    
    filename = 'scatter_plot_highresmip_assessment_'+year_period+basin+'hurdat2'+'.pkl'
    years_obs = np.arange(1979, 2015)
    print ('obs file ', os.path.join(directory,filename))
    if not os.path.exists(os.path.join(directory,filename) or do_new_obs):
        obs_storms = ts_obs.load_data.load_data()
        print(len(obs_storms))
        storm_type = ['TD','TS','HU','MH']
        obs_storms_hurr = list(_get_annual_storms(obs_storms, years_obs, months, basin))

        fh = open(os.path.join(directory,filename), 'wb') # write binary file mode
        pickle.dump(obs_storms_hurr, fh)
        fh.close()
    else:
#reading
        fh = open(os.path.join(directory,filename), 'rb') #read binary mode
        obs_storms_hurr = pickle.load(fh)
        fh.close()
    
#    for i, runid in enumerate(runid_model):
#        print 'hart runid ',runid
#        storms = storms_yr[runid]
#        plot_hart_lifecycle(storms)

    use_925 = False

    basins = ['na', 'wp', 'ep']
    fig_title = ['(a)','(b)','(c)']
    
    fig = plt.figure(figsize=[10, 10], dpi=150)
    plt.rcParams.update({'font.size': 8})    
    obs_storms_hurr_use = copy.copy(obs_storms_hurr)
    nmodels = 7
    for ib, basin in enumerate(basins):
        ax = fig.add_subplot(2,2,ib+1)
        plt.title('')

        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            print('mslp_wind runid ',runid, basin)
            if 'jra' in key:
                continue
            storms = storms_yr[key]
            if runid == 'jra55':
                plot_latitude_storm(storms)
            x_fit, y_fit = calculate_fit_model(storms, use_925 = use_925, basin = basin)
            if paired:
                if i%2 == 1 or i >= nmodels*2:
                    label = run_dict['model_desc'][i]
                else:
                    label = ''
            else:
                label = run_dict['model_desc'][i]
            plot_model_mslp_wspeed(storms, i, ax, x_fit=x_fit, y_fit=y_fit, label=label, use_925 = use_925, paired = paired, basin = basin, linewidth = 1.5, plot_scatter = True)
    
        x_fit_obs, y_fit_obs = calculate_fit_obs(obs_storms_hurr_use, basin = basin)
        plot_obs_mslp_wspeed(obs_storms_hurr_use, ax, x_fit=x_fit_obs, y_fit = y_fit_obs, basin = basin, linewidth = 1)
        ax.legend(loc=3, fontsize='x-small')
        setup_mslp_windspeed_plot(ax, basin = basin, fig_no = ib+1)    
        #plt.title(fig_title[ib]+' '+basin.upper(), loc = 'left', pad = 1.0)
        ax.text(-4, 1029, fig_title[ib]+' '+basin.upper(), fontsize=12)

    #plt.figure(fig.number)
    fig.subplots_adjust(bottom=0.05, left=0.07, right=0.93, top=0.94, wspace=0.15, hspace=0.15)
    plot_filename = plot_dir+'/intensity_mslp_wspeed_highresmip_resols_basins.pdf'
    plt.savefig(plot_filename)
    plt.savefig(plot_filename[:-3]+'png')
    plt.show()
    
    i = 0
    fig = plt.figure(figsize=[10,10], dpi=150)
    for ir in range(0, len(runid_model), 2):
        i += 1
        obs_storms_hurr_use = copy.copy(obs_storms_hurr)
        ax = fig.add_subplot(3, 3, i)
        plt.title('')
        for ib, basin in enumerate(basins):
            print('ir, ',ir)
            for im in range(2):
                if ir + im >= len(runid_model):
                    continue
                runid_s = runid_model[ir+im]
                resol = resols[ir+im]
                key = runid_model[ir+im]+resol
                print('mslp_wind runid ',key)
                storms = storms_yr[key]
                if runid == 'jra55':
                    plot_latitude_storm(storms)
                x_fit, y_fit = calculate_fit_model(storms, use_925 = use_925, basin = basin)
                #label = run_dict['model_desc'][i]
                label = resol+'_'+basin

                plot_model_mslp_wspeed_noscatter(storms, ir+im, ax, x_fit=x_fit, y_fit=y_fit, label=label, use_925 = use_925, paired = paired, basin = basin, linewidth = 1)
    
        x_fit_obs, y_fit_obs = calculate_fit_obs(obs_storms_hurr_use, basin = basin)
        plot_obs_mslp_wspeed(obs_storms_hurr_use, ax, x_fit=x_fit_obs, y_fit = y_fit_obs, basin = basin, linewidth = 1)
        ax.legend(loc=3, fontsize='x-small')
        setup_mslp_windspeed_plot(ax, basin = basin)    
        plt.suptitle(fig_title[ib]+' '+basin.upper(), x=0.1)
    
        plt.figure(fig.number)
    plot_filename = plot_dir+'/intensity_mslp_wspeed_highresmip_resols_subplots.pdf'
    #plt.savefig(plot_filename)
    plt.savefig(plot_filename[:-3]+'png')
        #plt.show()


    # calc and plot intensity-lifetime
    max_life_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        print(runid, resol)
        max_life_histogram[key], xedges, yedges = calculate_histogram_maxlife(storms_yr[key])
    obs_storms_hurr_use = copy.copy(obs_storms_hurr)
    max_life_histogram['hurdat'], xedges, yedges = calculate_histogram_maxlife(obs_storms_hurr_use)
    #print 'edges ',xedges
    plot_2d_histogram(max_life_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Min. MSLP', title = 'Max_intensity-lifetime', min_val = 0.0001)
    
    if T63_tracking:
    # calc and plot max wspeed-core strength
        maxv_core_histogram = {}
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            print(runid)
            maxv_core_histogram[key], xedges, yedges = calculate_histogram_warmcorevmax(storms_yr[key], runid)
        plot_2d_histogram(maxv_core_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Max vwind (925hPa)', title = 'Max_vwind-core', ylabel = 'Core strength', max_val = 0.003, min_val = 0.00002)
    
        # calc and plot 850 vorticity vs vorticity gradient (top - bottom)
        v850_vgrad_histogram = {}
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            print(runid)
            v850_vgrad_histogram[key], xedges, yedges = calculate_histogram_v850_vgrad(storms_yr[key], runid)
        plot_2d_histogram(v850_vgrad_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Max vort (850hPa)', title = 'Max_vort-core', ylabel = 'Core strength', max_val = 0.003, min_val = 0.00002)
    
    # calc and plot max wspeed-mslp min
    maxv_mslp_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        print(runid)
        maxv_mslp_histogram[key], xedges, yedges = calculate_histogram_vmaxmslp(storms_yr[key], runid)
    obs_storms_hurr_use = copy.copy(obs_storms_hurr)
    maxv_mslp_histogram['hurdat'], xedges, yedges = calculate_histogram_vmaxmslp(obs_storms_hurr_use, 'hurdat')
    plot_2d_histogram(maxv_mslp_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Max 925wind', title = 'Max_925wind-MSLP_min', ylabel = 'MSLP min', max_val = 0.003, min_val = 0.00002)
    
    # calc and plot intensity-latitude
    max_lat_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        max_lat_histogram[key], xedges, yedges = calculate_histogram_maxlat(storms_yr[key])
    max_lat_histogram['hurdat'], xedges, yedges = calculate_histogram_maxlat(obs_storms_hurr)
    print('edges ',xedges)
    plot_2d_histogram(max_lat_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'MSLP (hPa)', title = 'Max_intensity-latitude', ylabel = 'Latitude ($^{o}$N)', max_val = 0.003, min_val = 0.00002, show_title = False)
    
    # calc and plot lifetime-latitude
    life_lat_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        life_lat_histogram[key], xedges, yedges = calculate_histogram_lifelat(storms_yr[key])
    obs_storms_hurr_use = copy.copy(obs_storms_hurr)
    life_lat_histogram['hurdat'], xedges, yedges = calculate_histogram_lifelat(obs_storms_hurr)
    print('edges ',xedges)
    plot_2d_histogram(life_lat_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Lifetime', title = 'Lifetime-latitude_max_intensity', ylabel = 'Latitude (deg)', max_val = 0.003, min_val = 0.00002)
    
    # calc and plot lifetime-genesis-latitude
    life_genlat_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        life_genlat_histogram[key], xedges, yedges = calculate_histogram_lifegenlat(storms_yr[key])
    life_genlat_histogram['hurdat'], xedges, yedges = calculate_histogram_lifegenlat(obs_storms_hurr)
    print('edges ',xedges)
    plot_2d_histogram(life_genlat_histogram, runid_model, resols, model_resol, plot_dir, xedges = xedges, yedges = yedges, xlabel = 'Lifetime', title = 'Lifetime-genesis-latitude', ylabel = 'Latitude (deg)', max_val = 0.003, min_val = 0.00002)
    
    # plot min MSLP at peak
    maximum_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        maximum_histogram[key], edges = calculate_histogram_maxmslp(storms_yr[key])
    maximum_histogram['hurdat'], edges = calculate_histogram_maxmslp(obs_storms_hurr)
    print('edges ',edges)
    plot_latitude_histogram(maximum_histogram, runid_model, resols, model_desc, plot_dir, edges = edges, x_label = 'Min. MSLP', title = 'Mslp_min_at_peak_intensity', paired = paired)
    
    # plot max 925 windspeed at peak
    maximum_histogram = {}
    for basin in ['nh', 'wp', 'na']:
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            maximum_histogram[key], edges = calculate_histogram_maxvmax(storms_yr[key], basin = basin)
    #maximum_histogram['hurdat'], edges = calculate_histogram_maxvmax(obs_storms_hurr, scaling = knots_to_ms*(1.0/0.88))
        title = '(a) Max. 925 hPa windspeed in '+basin.upper()
        savefname = 'Max_925_windspeed_at_peak_vmax_'+basin
        plot_latitude_histogram(maximum_histogram, runid_model, resols, model_desc, plot_dir, edges = edges, x_label = 'Max. 925 hPa windspeed', title = title, paired = paired, savefname = savefname)
    
    # plot max 10m windspeed at peak
    for basin in ['nh', 'wp', 'na']:
        maximum_histogram = {}
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            storms = storms_yr[key]
            obs_storms_hurr_use = copy.copy(obs_storms_hurr)
            maximum_histogram[key], edges = calculate_histogram_maxvmax(storms, wind10m = True, basin = basin)
        maximum_histogram['hurdat'], edges = calculate_histogram_maxvmax(obs_storms_hurr_use, scaling = knots_to_ms, basin = basin)
        title = '(b) Max. 10m windspeed in '+basin.upper()
        savefname = 'Max_10m_windspeed_at_peak_vmax_'+basin
        plot_latitude_histogram(maximum_histogram, runid_model, resols, model_desc, plot_dir, edges = edges, x_label = 'Max. 10m windspeed', title = title, paired = paired, plot_cat = True, savefname = savefname)
    
    plt.rcParams.update({'font.size': 10})    
    # plot both max 925 and 10m windspeed at peak
    maximum_histogram = {}
    for basin in ['nh', 'wp', 'na']:
        fig = plt.figure(figsize=(8,10),dpi=100)
        ax = fig.add_subplot(2, 1, 1)
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            maximum_histogram[key], edges = calculate_histogram_maxvmax(storms_yr[key], basin = basin)
    #maximum_histogram['hurdat'], edges = calculate_histogram_maxvmax(obs_storms_hurr, scaling = knots_to_ms*(1.0/0.88))
        title = '(a) 925 hPa wind in '+basin.upper()
        savefname = 'Max_925_10m_windspeed_at_peak_vmax_'+basin
        plot_latitude_histogram_nofig(maximum_histogram, runid_model, resols, model_desc, plot_dir, ax, edges = edges, x_label = 'Lifetime max. 925 hPa windspeed', title = title, paired = paired, savefname = savefname, fig_no = 1)
    
        maximum_histogram = {}
        ax = fig.add_subplot(2, 1, 2)
        for i, runid in enumerate(runid_model):
            resol = resols[i]
            key = runid+resol
            storms = storms_yr[key]
            obs_storms_hurr_use = copy.copy(obs_storms_hurr)
            maximum_histogram[key], edges = calculate_histogram_maxvmax(storms, wind10m = True, basin = basin)
        maximum_histogram['hurdat'], edges = calculate_histogram_maxvmax(obs_storms_hurr_use, scaling = knots_to_ms, basin = basin)
        title = '(b) 10m wind in '+basin.upper()
        savefname = 'Max_925_10m_windspeed_at_peak_vmax_'+basin
        plot_latitude_histogram_nofig(maximum_histogram, runid_model, resols, model_desc, plot_dir, ax, edges = edges, x_label = 'Lifetime max. 10m windspeed', title = title, paired = paired, plot_cat = True, savefname = savefname, fig_no = 2)
    
        fig.subplots_adjust(bottom=0.1, left=0.08, right=0.93, top=0.94, wspace=0.1, hspace=0.2)
        plot_filename = plot_dir+'/intensity_scatter_latitude_'+savefname+'.pdf'
        plt.savefig(plot_filename)
        plt.savefig(plot_filename[:-3]+'png')
    #plt.show()    

    # plot latitude at peak
    latitude_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        latitude_histogram[key], edges = calculate_histogram_latitudes(storms_yr[key])
    latitude_histogram['hurdat'], edges = calculate_histogram_latitudes(obs_storms_hurr)
    plot_latitude_histogram(latitude_histogram, runid_model, resols, model_desc, plot_dir, edges = edges, title='Latitude at peak intensity', paired = paired)
    
    lifetime_histogram = {}
    for i, runid in enumerate(runid_model):
        resol = resols[i]
        key = runid+resol
        lifetime_histogram[key], edges = calculate_histogram_lifetime(storms_yr[key])
    lifetime_histogram['hurdat'], edges = calculate_histogram_lifetime(obs_storms_hurr)
    
    print('edges ',edges)
    plot_latitude_histogram(lifetime_histogram, runid_model, resols,  model_desc, plot_dir, edges=edges, x_label = 'Lifetime (6hr intervals)', title='Lifetime', paired = paired)
    

if __name__ == '__main__':
    years_amip = numpy.arange(1950,2015)
    months = [5,6,7,8,9,10,11]

    runid_erai = 'erai'; resol_erai = 'ERAI'
    runid_jra25 = 'jra25'; resol_jra25 = 'JRA25'
    runid_merra = 'merra'; resol_merra = 'MERRA'
    
    years = years_amip

    base_dir = '/home/users/mjrobert/track_148b/TRACK/results/'
    base_name = 'combined_ff_trs.vor_10m_fullgrid_'

    directory = '/home/users/mjrobert/track_148b/TRACK/results/intensity/'

    work(directory, years, months)
