import time
import numpy as np
import os
import pickle
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker 
import iris
import iris.plot as iplt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LATITUDE_FORMATTER,LONGITUDE_FORMATTER

import ts_obs
import tc_assess_code
import tc_interannual
import tc_intensity
import tc_calculate

TC_SCRIPT_PATH='/data/users/hadom/branches/git/eerie-project/storm_track_analysis/assess'
IBTRACS_FILE = os.path.join(TC_SCRIPT_PATH, 'ts_obs/obs_data', 'IBTrACS.since1980.v04r00.nc')

BASINS_LIST = tc_assess_code.BASINS_LIST
BASINS = tc_assess_code.BASINS
BASIN_NAME = tc_assess_code.BASIN_NAME
tropical_storm_threshold = tc_assess_code.tropical_storm_threshold
knots_to_ms = tc_assess_code.knots_to_ms

mslp_limits = [994,980,965,945,920,880]
wind_limits = [33,43,49,58,70,80]
mslp_limits_new = [880,920,945,965,980,995,1020]
wind_limits_new = [0,33,43,49,58,70,90]
wind_max = 80
hurr_cat = ['1','2','3','4','5']

colournames = ['red','blue','green','plum','aquamarine','coral','blueviolet','indianred','lime','lightgoldenrodyellow','pink','lightblue','plum']*3+['black']

def monthly_variability_with_obs(storms, years, months, runid_info, fontsize = 10, keys = [], years_obs=[]):
    """ 
    Plots monthly mean tropical storm counts for a set of
    model results and corresponding observations. Results
    for all ocean basins are plotted. Note: model counts
    assume only one ensemble member; further work would be 
    needed to normalise with respect to ensemble count.
    
    
    """
    basins = ['na','ep','wp','ni','si','au','sp']
    
    # Load observations and get monthly storm count

    fig = plt.figure(figsize=[7,9])
    
    for i, basin in enumerate(basins): 
        obs_storms = ts_obs.load_data.load_data(basin=basin)       
        obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
        obs_monthly_count = tc_assess_code.get_monthly_mean_count(obs_storms_tsandup, years_obs, months, basin)
        #Set up plot with multiple panels
        plt.subplot(4,2,i+1)

        # Plot bar plot for observations and line for model results
        plt.bar(np.arange(1, len(months)+1), obs_monthly_count, 
                align='center', color='black', label='Observations')

        for ir, runid in enumerate(runid_info['model']):
            for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
                key = keys[ir]
                # Get monthly storm count for the model  
                nyears = storms[(runid, key, 'nyears')]
                model_monthly_count = tc_assess_code.get_monthly_mean_count(storms[(runid, key)], years, months, basin, nyears = nyears)
                #print('monthly ',runid, basin, model_monthly_count)
        
                plt.plot(months, model_monthly_count, linewidth=3,
                         label=runid+'-'+algo, c=colournames[ir])        
        
                # Add plot details and titles etc.
                plt.title('%s' % tc_assess_code.BASIN_NAME.get(basin), fontsize = fontsize)  
                xlabels = tc_assess_code._month_names(months)  
                xlab_letter = []
                [xlab_letter.append(c[0:1]) for c in xlabels]

                plt.xticks(np.arange(1,len(months)+1), xlab_letter, fontsize = fontsize)
                plt.margins(0.05)  
                plt.ylim(0, 8)
                plt.ylabel('No. TCs', fontsize = fontsize - 2)

    # Add extra space around plots, leaving space for the title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Add main title for plot
    plt.suptitle('Monthly Mean Tropical Storm Frequency \n %s' % 
                 tc_assess_code._get_time_period(years, months), 
                 fontsize=13
                 )  

    # Add legend in bottom right corner        
    plt.legend(loc=8, ncol=1, numpoints=1, bbox_to_anchor = (0.8, 0.0), 
               bbox_transform=plt.gcf().transFigure, 
               prop={'size': fontsize}, 
               frameon=False
               )

    fig.subplots_adjust(bottom=0.07, left=0.06, right=0.95, top=0.90, wspace=0.14, hspace=0.4)


def storm_tracks(storms, years, months, basin, title, fig, ax, algorithm, hemi, genesis=False, lysis=False, max_intensity=False, warmcore = False, yoff=0., fontsize = 12):
    """ 
    Plots storm tracks, genesis locations, lysis locations 
    or location of maximum intensity for all model tropical
    storms that form in a desired set of years, months and 
    ocean basin. Default plot is of storm tracks. 
    
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
    1996 and 30 April 1997 inclusive, assuming those 
    data are available in the input file.
    
    Setting the basin will return any storms which
    *passed through* (not just formed) in the 
    designated basin region. 
    
    
    """   
    count = 0
    for year in years:
        for storm in tc_assess_code._storms_in_time_range(storms, year, months):
            if genesis:
                variable = 'Genesis'
                ax.plot(storm.obs_at_genesis().lon, storm.obs_at_genesis().lat,
                         'bo', markersize=3, transform=ccrs.Geodetic())

            elif lysis:
                variable = 'Lysis'
                ax.plot(storm.obs_at_lysis().lon, storm.obs_at_lysis().lat,
                         'go', markersize=3, transform=ccrs.Geodetic())

            elif max_intensity:
                variable = 'Maximum Intensity'
                ax.plot(storm.obs_at_vmax().lon, storm.obs_at_vmax().lat,
                         'ro', markersize=3, transform=ccrs.Geodetic())

            elif warmcore:
                variable = 'Tracks'  
                if storm.warm_core:
                    ax.plot([ob.lon for ob in storm.obs], [ob.lat for ob in storm.obs],
                     linewidth=1.2, transform=ccrs.Geodetic())
            else:
                variable = 'Tracks'  
                ax.plot([ob.lon for ob in storm.obs], [ob.lat for ob in storm.obs],
                     linewidth=1.2, transform=ccrs.Geodetic())
            count += 1

    if count != 0:
        fig.gca().coastlines() 
        title1 = 'Model Tropical Storm %s\n (%s - %s) using %s \n %s' % \
                  (variable, str(years[0]), str(years[-1]), algorithm, title)
        #print(title1)
        ax.set_title(title1, fontsize=fontsize)
        s = ('Total tropical storms for %s: %s' % (hemi, count))
        plt.text(0.02, -0.08-yoff, s, ha='left', va='center', transform=plt.gca().transAxes)
        #print s   
    else:
        print('No storms found')

def set_colour(vmax):
    """ 
    Returns colour based on wind speed (kts). Names from
    xkcd https://xkcd.com/color/rgb/. 
    
    """
    if vmax < 34:
        return 'lightsteelblue'# dark sky blue   
    elif vmax >= 34 and vmax < 64:
        return 'lightsteelblue' #bright cyan
    elif vmax >= 64 and vmax < 83:
        return 'skyblue' #cream
    elif vmax >= 83 and vmax < 96:
        return 'palegreen' #light gold
    elif vmax >= 96 and vmax < 113:
        return 'peachpuff' #goldenrod
    elif vmax >= 113 and vmax < 137:
        return 'darkorange' #pumpkin orange
    elif vmax >= 137:
        return 'red' #fire engine red      
    else:
        raise UserWarning('Wind speed has no associated colour', 
                         vmax)

def set_colour_mslp(vmax):
    """ 
    Returns colour based on mslp. Names from
    xkcd https://xkcd.com/color/rgb/. 
    
    """
    if vmax >= 994.0:
        return 'lightsteelblue'# dark sky blue   
    elif vmax >= 980.0 and vmax < 994.0:
        return 'skyblue' #cream
    elif vmax >= 946.0 and vmax < 980.0:
        return 'palegreen' #light gold
    elif vmax >= 945.0 and vmax < 965.0:
        return 'peachpuff' #goldenrod
    elif vmax >= 920.0 and vmax < 945.0:
        return 'darkorange' #pumpkin orange
    elif vmax < 920.0 and vmax > 0:
        return 'red' #fire engine red     
    elif vmax == -999.0:
        return 'lightsteelblue'  
    else:
        raise UserWarning('Wind speed has no associated colour', 
                         vmax)
 
def storm_tracks_intensity(storms, years, months, basin, title, fig, ax, algorithm, hemi, genesis=False, lysis=False, max_intensity=False, warmcore = False, yoff=0., fontsize = 12, obs = False):

    if not obs:
        factor = 1.944
    else:
        factor = 1.0

    count = 0
    for year in years:
        for storm in tc_assess_code._storms_in_time_range(storms, year, months):
            variable = 'Tracks'  
            #ax.plot([ob.lon for ob in storm.obs], [ob.lat for ob in storm.obs],
            #         linewidth=1.2, transform=ccrs.Geodetic(), color = 'lightblue')
            lon = [ob.lon for ob in storm.obs]
            lat = [ob.lat for ob in storm.obs]
            mslp = [ob.mslp for ob in storm.obs]
            if min(mslp) < 920:
                print('mslp < 920 ',min(mslp), lon[0], lat[0], storm.genesis_date())
            if not obs:
                for i in range(len(storm.obs)-1):
                    ax.plot([lon[i], lon[i+1]], [lat[i], lat[i+1]], color = set_colour_mslp(mslp[i]), 
                     linewidth=1.2, transform=ccrs.Geodetic())
            else:
                vmax= [ob.vmax for ob in storm.obs]
                for i in range(len(storm.obs)-1):
                    ax.plot([lon[i], lon[i+1]], [lat[i], lat[i+1]], color = set_colour(vmax[i]), 
                     linewidth=1.2, transform=ccrs.Geodetic())
            #plt.scatter([ob.lon for ob in storm.obs], [ob.lat for ob in storm.obs],
            #        color=[set_colour_mslp(ob.mslp) for ob in storm.obs], 
            #        transform=ccrs.Geodetic(), zorder=10, 
            #        edgecolor='black', s=10)
            count += 1

    if count != 0:
        fig.gca().coastlines() 
        title1 = 'Model Tropical Storm %s\n (%s - %s) using %s \n %s' % \
                  (variable, str(years[0]), str(years[-1]), algorithm, title)
        #print(title1)
        ax.set_title(title1, fontsize=fontsize)
        s = ('Total tropical storms for %s: %s' % (hemi, count))
        plt.text(0.02, -0.08-yoff, s, ha='left', va='center', transform=plt.gca().transAxes)
        #print s   
    else:
        print('No storms found')


def storm_density_plot(cube, years, fig, ax, basin=None, title = '', fontsize = 12, plot_diff = False, genesis = False, bias = False, factor = 1.0):
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
    # Plot map and contour plot of track density
    #fig = plt.figure(figsize=(9,6), dpi=100)
    #ts_model.example_code.load_map(basin=basin)
    #if cube.data.max() > 1:
    #    maxval = np.max([(cube.data.max() // 1), 1])
    #else:
    #    maxval = cube.data.max()
    #maxval = 2.5
    #interval = maxval / 10
    #minval = interval
    #print ('maxval, interval ',maxval, interval)
    #levels = np.arange(minval, maxval, interval)
    
    if not plot_diff:
        cmap = matplotlib.colors.ListedColormap(['#ffffff','#d8cecb','#b39d99','#8c6c66','#977b61','#beab75','#e3dc8a','#e8fa94','#9ceb85','#51dc76',\
                                             '#05cc66','#04a8d0','#0f7ae0','#333399'])
        if genesis:
            levels = [0.001, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.075, 0.1]
        else:
            levels = [0.05,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5,1.75,2.0,3.,4.]
    else:
        name = 'bwr'
        cmap = plt.get_cmap(name)
        if genesis:
            levels = np.arange(-0.05, 0.052, 0.002)
            levels = [-0.05, -0.04, -0.03, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]
        else:
            levels = np.arange(-1.2, 1.3, 0.1)
            levels = [-1.25, -1.0, -0.75, -0.5, -0.25, -0.1, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25]
    if factor != 1.0:
        levels = [lev*factor for lev in levels]

    norml = matplotlib.colors.BoundaryNorm(levels, cmap.N)
    cs = iplt.pcolormesh(cube, cmap=cmap, norm=norml)
#    cs = iplt.contourf(cube, 
#                       levels=levels, 
#                       extend='max', 
#                       cmap='Reds')
    plt.gca().coastlines() 
    gl = ax.gridlines(draw_labels=True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    #gl.top_labels = False
    #gl.right_labels = False
    gl.yformatter = LATITUDE_FORMATTER
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlocator = mticker.FixedLocator([0, 90, -180, -90, 270])
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}

    # Add colorbar
    cb = plt.colorbar(cs, orientation='horizontal')
    cb.set_label('Storm transits per month', fontsize=fontsize-2)
    
    # Add plot title and storm count 
    if bias:
        bias_title = ' bias'
    else:
        bias_title = ''
    if genesis:
        title_new = 'Storm Track Genesis%s, %s' % (bias_title, tc_assess_code.BASIN_NAME.get(basin))
    else:
        title_new = 'Storm Track Density%s, %s' % (bias_title, tc_assess_code.BASIN_NAME.get(basin))
    title_new += '\n '+title
    plt.title(title_new, fontsize=fontsize)
    #s = ('Total storms: %s' % (count))
    #plt.text(0.02, 0.05, s, ha='left', va='center', 
    #         transform=plt.gca().transAxes)
    #plt.show()

def plot_piechart(dict_in, basins, ax, fig, obs_count, title = '', plot_dir = './', plot_key = False, fontsize = 12):
    '''
    See
    https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.pyplot.pie.html
    '''
    colors = ['#ff6666', '#ffcc99', '#cc9966', 'cc6666', '#99ff99', '#66b3ff', '#c2c2f0', '#ff6666', '#ffcc99', '#99ff99', '#66b3ff', '#c2c2f0']
    colors = ['#ff6666', '#ffcc99', '#cc9966', '#cc6666', '#66b3ff', '#c2c2f0', '#6666ff', '#33cccc']
    colors = ['#ff6666', '#ffcc99', '#cc9966', '#cc6666', '#66b3ff']
    #colors_gender = ['#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6', '#c2c2f0','#ffb3e6']
    basins = ['na','wp','ep','ni','other']
    labels = basins
    labelnames = ['NA','WP','EP','NI','Other']
 
# Plot
    sizes = []; labels = []
    for ib, basin in enumerate(basins):
        if basin != 'nh' and basin != 'sh' and basin != 'other':
            sizes.append(dict_in[basin])
            labels.append(labelnames[ib])
    #print('nh total, sum(sizes) ',np.mean(dict_in['nh']), np.sum(sizes))
    rest = np.mean(dict_in['nh']) - np.sum(sizes)
    if rest < 0:
        rest = 0.
    sizes.append(rest)
    labels.append('Other')
    #total_storms = dict_in['nh'][0] + dict_in['sh'][0]
    total_storms = dict_in['nh']
    total_storms_sh = dict_in['sh']
    total_storms_obs = obs_count['nh']

    radius_ref = 1.5
    radius = (total_storms / total_storms_obs) * radius_ref
    radius = np.max([radius, 0.53])

    #plt.rc('font', size=fontsize) 
    #print ('sizes, total_storms ',sizes, total_storms)
    #print (len(labels), len(basins), len(colors))
    wedges, texts, autotexts = ax.pie(sizes, colors=colors, startangle=90, frame=True, autopct='%1.0f', radius = radius, labeldistance = 1.0, pctdistance = 0.7, textprops = dict(color='black'))
    centre_circle = plt.Circle((0,0), 0.5, color='black', fc='white', linewidth=0)
    #fig = plt.gcf()
    fig.gca().add_artist(centre_circle)
    #model = runid.split('_')[0]
    #plt.title(resol+ ' '+model+', '+str(total_storms)[:5])
    title = title
    fontdict = {'fontsize': fontsize, 'fontweight' : 4}
    plt.title(title, loc = 'left', fontdict = fontdict)
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    plt.axis('equal')
    plt.axis('off')
    #plt.text(-0.1, 0.1, str(total_storms)[:5])
    yoffset = 0.13 * (radius / radius_ref) 
    xoffset = -0.28 * radius / radius_ref
    plt.text(xoffset, yoffset - 0.09, str(round(total_storms, 1)), fontsize = fontsize, fontweight = 70)
    plt.text(xoffset, -yoffset - 0.09, str(round(total_storms_sh, 1)), fontsize = fontsize, fontweight = 70)
    #plt.tight_layout()

    if plot_key:
        ax.legend(wedges, labels,
          title="NH Basin",
          loc="center left",
                  bbox_to_anchor=(0.9, 0, 0.25, 1.0), fontsize = fontsize-2, borderaxespad=0)

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
    hist, xedges, yedges = numpy.histogram2d(vmax_storm, min_mslp, bins = (xedges, yedges), normed = True)
    hist = hist.T
    return hist, xedges, yedges
            
def calculate_histogram_maxvmax(storms, vmax_range, scaling = 1.0, variable = 'wind10m', basin = 'nh', obs = False, data_freq = 6):
    '''
    Calculate a histogram for the input variable from the tracks
    '''

    maximum_storm = []
    for storm in storms:
        max_strength = 0.
        if np.absolute(storm.obs_at_vmax().lat) <= 40.0:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                if variable == '10mwind':
                    if obs:
                        max_strength = storm.obs_at_vmax().vmax*scaling
                    else:
                        max_strength = storm.obs_at_vmax().extras['w10m']
                elif variable == 'slp' or variable == 'mslp':
                    if storm.obs_at_vmax().mslp > 600 and storm.obs_at_vmax().mslp < 2000 :
                        max_strength = storm.obs_at_vmax().mslp
                elif variable == 'vmax':
                    max_strength = storm.obs_at_vmax().vmax
                elif variable == 'rvT63':
                    max_strength = storm.obs_at_vmax().vort
                    if max_strength == 1.0:
                        if storm.obs_at_genesis().lat > 0:
                            max_strength = storm.obs_at_vmax().vort
                        else:
                            max_strength = storm.obs_at_vmax().vort_min*-1.0
                    if np.abs(max_strength) < 1e-1:
                        max_strength *= 1.0e5
                    #if basin == 'si':
                        #print ('basin, vort min ', basin, max_strength, storm.obs_at_vmax().vort, storm.obs_at_vmax().lat, storm.obs_at_vmax().lon)
                        #vorts = [ob.vort for ob in storm.obs]
                        #print (vorts)
                elif variable == 'latitude':
                    max_strength = storm.obs_at_vmax().lat
                elif variable == 'lifetime':
                    freq_convert_to_days = 24.0 / data_freq
                    max_strength = len(storm.obs)/freq_convert_to_days

#                if max_strength > 0.:
                maximum_storm.append(max_strength)

    hist, bin_edges = np.histogram(maximum_storm, bins = vmax_range, range = (vmax_range[0], vmax_range[-1]), density = True)
    return hist, bin_edges

def calculate_histogram_maxintens(storms, vmax_range, scaling = 1.0, variable = 'wind10m', basin = 'nh', obs = False, data_freq = 6):
    '''
    Calculate a histogram for the maximum intensification for each storm
    '''
    min_storm_len = int(3 * 24 / data_freq)
    day_len = int(24 / data_freq)
    maximum_storm = []
    for storm in storms:
        max_strength = 0.
        delta_max = 0.
        if np.absolute(storm.obs_at_vmax().lat) <= 40.0:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                if variable == '10mwind':
                    values = []
                    if obs:
                        for ob in storm.obs:
                            values.append(ob.vmax * scaling)
                        #print('obs ',values)
                    else:
                        for ob in storm.obs:
                            values.append(ob.extras['w10m'])
                    if len(values) > min_storm_len and np.min(values) > 0:
                        delta = np.roll(np.asarray(values), -day_len) - np.asarray(values) 
                        delta_max = np.amax(delta[day_len:-day_len])
                        delta_max = np.amax([delta_max, 0.0])
                        #if delta_max > 40:
                        #    print('big change ',values, delta_max)
                elif variable == 'slp' or variable == 'mslp':
                    #print('storm max ',storm.obs_at_vmax().mslp)
                    if storm.obs_at_vmax().mslp > 600 and storm.obs_at_vmax().mslp < 2000 :
                        values = []
                        for ob in storm.obs:
                            values.append(ob.mslp)
                        if len(values) > min_storm_len and np.min(values) > 800:
                            delta = np.roll(np.asarray(values), -day_len) - np.asarray(values) 
                            delta_max = np.amin(delta[day_len:-day_len])
                            delta_max = np.amin([delta_max, 0.0])
                            #if delta_max < -80:
                            #    print('big change ',values, delta_max)
                elif variable == 'vmax':
                    values = []
                    for ob in storm.obs:
                        values.append(ob.vmax)
                    delta = np.roll(np.asarray(values), -day_len) - np.asarray(values) 
                    delta_max = np.amax(delta[day_len:-day_len])

                if delta_max != 0.:
                    maximum_storm.append(delta_max)

    hist, bin_edges = np.histogram(maximum_storm, bins = vmax_range, range = (vmax_range[0], vmax_range[-1]), density = True)
    return hist, bin_edges

def calculate_histogram_maxintens_24h(storms, vmax_range, scaling = 1.0, variable = 'wind10m', basin = 'nh', obs = False, data_freq = 6):
    '''
    Calculate a histogram for the maximum intensification for each storm for each 24h period
    '''

    min_storm_len = int(2 * 24 / data_freq)
    day_len = int(24 / data_freq)
    maximum_storm = []
    inten_rate = []
    for storm in storms:
        max_strength = 0.
        delta_max = 0.
        if np.absolute(storm.obs_at_vmax().lat) <= 40.0:
            if tc_assess_code._storm_vmax_in_basin(storm, basin):
                if variable == '10mwind':
                    values = []
                    if obs:
                        for ob in storm.obs:
                            values.append(ob.vmax * scaling)
                        #print('obs ',values)
                    else:
                        for ob in storm.obs:
                            values.append(ob.extras['w10m'])
                    if len(values) > min_storm_len:
                        for iv in np.arange(day_len, len(values), day_len):
                            if values[iv] > 0. and values[iv-day_len] > 0:
                                delta = values[iv] - values[iv-day_len]
                                inten_rate.append(delta)
                        #if delta_max > 40:
                        #    print('big change ',values, delta_max)
                elif variable == 'slp' or variable == 'mslp':
                    #print('storm max ',storm.obs_at_vmax().mslp)
                    if storm.obs_at_vmax().mslp > 600 and storm.obs_at_vmax().mslp < 2000 :
                        values = []
                        for ob in storm.obs:
                            values.append(ob.mslp)
                        if len(values) > min_storm_len and np.min(values) > 800:
                            for iv in np.arange(day_len, len(values), day_len):
                                delta = values[iv] - values[iv-day_len]
                                inten_rate.append(delta)
                                                     
                elif variable == 'vmax':
                    values = []
                    for ob in storm.obs:
                        values.append(ob.vmax)
                    for iv in np.arange(day_len, len(values), day_len):
                        if values[iv] > 0. and values[iv-day_len] > 0:
                            delta = values[iv] - values[iv-day_len]
                            inten_rate.append(delta)

    hist, bin_edges = np.histogram(inten_rate, bins = vmax_range, range = (vmax_range[0], vmax_range[-1]), density = True)
    return hist, bin_edges

def plot_latitude_histogram_nofig(histogram, runid, algo, ax, edges = [], x_label = 'Latitude', title='', plot_cat = True, ymax = 1.0, color = ''):

    bins = edges[:-1]
    #print 'bins ',bins
    offset = (bins[1] - bins[0])/2.0

    #if 'hurdat' in histogram:
    #    plt.plot(bins+offset,histogram['hurdat'][:], color='black', label='OBS', lw =3)

    if color == '':
        plt.plot(bins+offset,histogram[:], label=runid+algo, lw =3)
    else:
        plt.plot(bins+offset,histogram[:], label=runid+algo, lw =3, color = color)
    ax.set_ylabel('Normalised frequency')
    ax.set_xlabel(x_label)
    ax.legend(loc="upper right", fontsize='medium', fancybox = True, framealpha = 0.8)
    plt.title(title, loc = 'left')
    #ymin, ymax = ax.get_ylim()
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

def do_storm_tracks(storms, years, runid_info, plot_dir, months_nh, months_sh, plot_matched = False, title_extra = '', plot_name_extra = '', warmcore = False, keys = [], regionx = [], regiony = [], years_obs=[]):

    if plot_matched:
        algo_types_match = ['tc_vort_T63', 'tc_vort_T63wc']
        runid_info_match = {'model': [runid_info['model'][0]], 'resol': resol, 'grid': model_grid, 'algorithm': algorithm, 'algo_type': algo_types_match}
        storms_match = track_matching.work(runid_info_match, CMOR_DIR, years, years)

    if regionx != [] or regiony != []:
        central_longitude = 0.
        plot_name_region = '_region'
        hspace = 0.4
        wspace = 0.3
        bottom = 0.1
        top = 0.90
    else:
        central_longitude = -160.
        plot_name_region = ''
        hspace = 0.05
        wspace = 0.08
        bottom = 0.02
        top = 0.95

    fig = plt.figure(figsize=(11,7), dpi=100)
    ny = (len(runid_info['model']) // 2)+1
    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            key = keys[ir]
            #print ('storm tracks ',ir, ia)
            ax = fig.add_subplot(ny,2,ir+ia+1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
            ax.set_global()
            yoff = 0
            for hemi in ['nh', 'sh']:
                storms_hemi = []
                #print ('algos ',runid, algo, hemi)
                for storm in storms[(runid, key)]:
                    if warmcore:
                        if not storm.warm_core:
                            continue
                    if storm.obs_at_genesis().lat >= 0.0 and hemi == 'nh':
                        storms_hemi.append(storm)
                        months = months_nh
                    elif storm.obs_at_genesis().lat < 0.0 and hemi == 'sh':
                        storms_hemi.append(storm)
                        months = months_sh

                    if hemi == 'sh': yoff = 0.08

                title = runid+'-'+runid_info['resol'][ir]+', '+algo+' '+title_extra
                if len(storms_hemi) > 0:
                    #print('years to plot tracks ',years)
                    storm_tracks_intensity(storms_hemi, years, months, hemi, title, fig, ax, algo, hemi, yoff=yoff, max_intensity = False, genesis = False, warmcore = warmcore)
            if regionx != [] or regiony != []:
                if regionx != []:
                    ax.set_xlim([regionx[0], regionx[1]])
                if regiony != []:
                    ax.set_ylim([regiony[0], regiony[1]])

    if plot_matched:
        ax = fig.add_subplot(2,2,4,projection=ccrs.PlateCarree(central_longitude=central_longitude))
        ax.set_global()
        yoff = 0
        for hemi in ['nh', 'sh']:
            storms_hemi = []
            for storm in storms_match:
                if storm.obs_at_genesis().lat > 0.0 and hemi == 'nh':
                    storms_hemi.append(storm)
                    months = months_nh
                elif storm.obs_at_genesis().lat < 0.0 and hemi == 'sh':
                    storms_hemi.append(storm)
                    months = months_sh

                if hemi == 'sh': yoff = 0.08

            title = runid+'-'+runid_info['resol'][ir]+', '+algo
            storm_tracks_intensity(storms_hemi, years, months, hemi, title, fig, ax, algorithm, hemi, yoff=yoff, max_intensity = False, genesis = False)

    # do observed tracks
    if not plot_matched:
        obs_storms = ts_obs.load_data.load_data()       
        obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
        yoff = 0
        ax = fig.add_subplot(2,2,4,projection=ccrs.PlateCarree(central_longitude=central_longitude))
        for hemi in ['nh', 'sh']:
            storms_hemi = []
            for storm in obs_storms_tsandup:
                if storm.obs_at_genesis().lat >= 0.0 and hemi == 'nh':
                    storms_hemi.append(storm)
                    months = months_nh
                elif storm.obs_at_genesis().lat < 0.0 and hemi == 'sh':
                    storms_hemi.append(storm)
                    months = months_sh

            if hemi == 'sh': yoff = 0.08

            title = 'OBS  storms '
            storm_tracks_intensity(storms_hemi, years_obs, months, hemi, title, fig, ax, '', hemi, yoff=yoff, max_intensity = False, genesis = False, warmcore = warmcore, obs = True)
        if regionx != [] or regiony != []:
            if regionx != []:
                ax.set_xlim([regionx[0], regionx[1]])
            if regiony != []:
                ax.set_ylim([regiony[0], regiony[1]])

    fig.subplots_adjust(bottom=bottom, left=0.05, right=0.95, top=top, wspace=wspace, hspace=hspace)
    #current_dir = os.getcwd()
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_tracks'+plot_name_extra+plot_name_region+'.png')
    if warmcore:
        figname = figname[:-4]+'_warmcore.png'
    
    plt.savefig(figname)
    plt.close()
    print(years)

def callback_del_attributes(cube, field, filename):
    '''
    Remove metadata that iris is not happy about
    '''
    attributes = ['start_date_nh', 'end_date_nh', 'start_date_sh', 'end_date_sh']
    for attrib in attributes:
        try:
            del(cube.attributes[attrib])
        except:
            pass

def do_track_density_plot(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis=False, keys = [], regionx = [], regiony = [], years_obs = [], storedir='./', basin=None):
    '''
    Firstly calculate the track density for each year and save
    Then read and mean these files for plotting
    '''
    cube_density = {}
    if len(runid_info['model']) > 1:
        nx = len(runid_info['model']); ny = 2
    else:
        nx = len(runid_info['algo_type']); ny = 2

    if nx > 2:
        fig = plt.figure(figsize=(12,8), dpi=100)
    else:
        fig = plt.figure(figsize=(8,8), dpi=100)

    if regionx != [] or regiony != []:
        central_longitude = 0.
        plot_name_region = '_region'
        hspace = 0.4
        wspace = 0.3
        bottom = 0.1
        top = 0.90
    else:
        central_longitude = -160.
        plot_name_region = ''
        hspace = 0.26
        wspace = 0.0
        bottom = 0.08
        top = 0.92

    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            ens = runid_info['ens'][ir]
            if ens == '':
                store_dir = storedir.format(runid, algo)
            else:
                store_dir = storedir.format(runid+'/'+ens, algo)
            key = keys[ir]
            ax = fig.add_subplot(ny,nx,ir+ia+1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
            storms_runid = storms[(runid, key)]
            nyears = storms[(runid, key, 'nyears')]
            years_model = storms[(runid, key, 'years')]
            #cube_nh, count_nh = storm_density_calculation(storms_runid, years, months_nh, 'nh', genesis=genesis, nyears = nyears)
            #cube_sh, count_sh = storm_density_calculation(storms_runid, years, months_sh, 'sh', genesis=genesis, nyears = nyears)
            #cube = cube_nh + cube_sh
            cube_files = []
            print('years model for density ',years_model)
            for year in years_model:
                if genesis:
                    fname = 'genesis'
                else:
                    fname = 'density'
                cube_file = os.path.join(store_dir, fname+'_'+str(year)+'.nc')
                if os.path.exists(cube_file):
                    cube_files.append(cube_file)
                else:
                    print('no file found for density ',cube_file)
            #print('cube list in density plot ',runid, algo, cube_files)
            if len(cube_files) < 1:
                raise Exception('No cube list files')
            cube_list = iris.cube.CubeList()
            for ic, f in enumerate(cube_files):
                c = iris.load_cube(f, callback=callback_del_attributes)
                #index_coord = iris.coords.AuxCoord([ic], long_name = 'Index', units = '1')
                #c.add_aux_coord(index_coord)
                cube_list.append(c)
            #cube_l = cube_list.merge_cube()
            cube_l = cube_list.concatenate_cube()
            #cube = cube_l.collapsed('Index', iris.analysis.MEAN)
            cube = cube_l.collapsed('time', iris.analysis.MEAN)
            #cube = iris.load_cube(cube_file)
            cube_density[(runid, key)] = cube
            title = runid+'-'+runid_info['resol'][ir]+', '+str(years_model[0])+'-'+str(years_model[-1])+', '+algo
            # set multiplicative factor for density plot
            factor = 1.0
            if '_ew' in algo:
                factor = 3.0
            if regionx != [] or regiony != []:
                if regionx != []:
                    ax.set_xlim([regionx[0], regionx[1]])
                if regiony != []:
                    ax.set_ylim([regiony[0], regiony[1]])
            storm_density_plot(cube, years_model, fig, ax, basin=None, title = title, genesis=genesis, factor = factor)

    if len(runid_info['model']) > 1:
        for ir, runid in enumerate(runid_info['model'][1:]):
            model_exp = runid
            model_ctl = runid_info['model'][0]
            #algo_ctl = runid_info['algo_type'][0]
            #algo_exp = runid_info['algo_type'][ir+1]
            algo_ctl = keys[0]
            algo_exp = keys[ir+1]
            cube_diff = cube_density[(model_exp, algo_exp)] - cube_density[(model_ctl, algo_ctl)]
            ax = fig.add_subplot(ny,nx,nx+ir+1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
            title = model_exp + '-' + model_ctl
            if regionx != [] or regiony != []:
                if regionx != []:
                    ax.set_xlim([regionx[0], regionx[1]])
                if regiony != []:
                    ax.set_ylim([regiony[0], regiony[1]])
            storm_density_plot(cube_diff, years, fig, ax, basin=None, title = title, plot_diff = True, genesis=genesis)
    elif len(runid_info['algo_type']) > 1 and len(runid_info['model']) == 1:
        for ia, algo in enumerate(runid_info['algo_type'][1:]):
            model_ctl = runid_info['model'][0]
            #algo_ctl = runid_info['algo_type'][0]
            #algo = runid_info['algo_type'][ia]
            algo_ctl = keys[0]
            algo = keys[ia]
            cube_diff = cube_density[(model_ctl, algo)] - cube_density[(model_ctl, algo_ctl)]
            ax = fig.add_subplot(ny,nx,nx+ia+1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
            title = model_ctl + ', '+algo+'-' + algo_ctl
            # set multiplicative factor for density plot
            factor = 1.0
            if 'ew' in algo:
                factor = 3.0
            if regionx != [] or regiony != []:
                if regionx != []:
                    ax.set_xlim([regionx[0], regionx[1]])
                if regiony != []:
                    ax.set_ylim([regiony[0], regiony[1]])
            storm_density_plot(cube_diff, years, fig, ax, basin=None, title = title, plot_diff = True, genesis=genesis, factor = factor)

    # do observed density
    obs_storms = ts_obs.load_data.load_data(basin=basin)       
    obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
    track_density_nh,count = tc_calculate.storm_density_calculation(obs_storms_tsandup, years_obs, months_nh, 'nh', genesis=genesis)
    track_density_sh,count = tc_calculate.storm_density_calculation(obs_storms_tsandup, years_obs, months_sh, 'sh', genesis=genesis)
    cube = track_density_nh + track_density_sh
    ax = fig.add_subplot(ny,nx,ny*nx,projection=ccrs.PlateCarree(central_longitude=central_longitude))
    title = 'OBS, '+str(years_obs[0])+'-'+str(years_obs[-1])
    if regionx != [] or regiony != []:
        if regionx != []:
            ax.set_xlim([regionx[0], regionx[1]])
        if regiony != []:
            ax.set_ylim([regiony[0], regiony[1]])
    storm_density_plot(cube, years_obs, fig, ax, basin=None, title = title, genesis=genesis)
    
    fig.subplots_adjust(bottom=0.08, left=0.05, right=0.95, top=0.92, wspace=0.0, hspace=0.26)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    if genesis:
        fout = models+'-'+algos+plot_name_region+'_genesis.png'
    else:
        fout = models+'-'+algos+plot_name_region+'_density.png'
    figname = os.path.join(plot_dir, fout)

    plt.savefig(figname)
    plt.close()

def do_track_density_bias(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis=False, keys = [], regionx = [], regiony = [], years_obs=[], basin=None):
    cube_density = {}
    if len(runid_info['model']) > 1:
        nx = len(runid_info['model']); ny = 2
    else:
        nx = len(runid_info['algo_type']); ny = 2

    if nx > 2:
        fig = plt.figure(figsize=(12,8), dpi=100)
    else:
        fig = plt.figure(figsize=(8,8), dpi=100)

    if regionx != [] or regiony != []:
        central_longitude = 0.
        plot_name_region = '_region'
        hspace = 0.4
        wspace = 0.3
        bottom = 0.1
        top = 0.90
    else:
        central_longitude = -160.
        plot_name_region = ''
        hspace = 0.26
        wspace = 0.0
        bottom = 0.08
        top = 0.92

    # do observed density
    obs_storms = ts_obs.load_data.load_data(basin=basin)       
    obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
    track_density_nh,count = tc_calculate.storm_density_calculation(obs_storms_tsandup, years_obs, months_nh, 'nh', genesis=genesis)
    track_density_sh,count = tc_calculate.storm_density_calculation(obs_storms_tsandup, years_obs, months_sh, 'sh', genesis=genesis)
    cube_obs = track_density_nh + track_density_sh

    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            key = keys[ir]
            ax = fig.add_subplot(ny,nx,ir+ia+1,projection=ccrs.PlateCarree(central_longitude=central_longitude))
            storms_runid = storms[(runid, key)]
            nyears = storms[(runid, key, 'nyears')]
            years_model = storms[(runid, key, 'years')]
            data_freq = storms[(runid, key, 'freq')]
            #print('density bias ',years)
            cube_nh, count_nh = tc_calculate.storm_density_calculation(storms_runid, years, months_nh, 'nh', genesis=genesis, nyears = nyears, data_freq=data_freq)
            cube_sh, count_sh = tc_calculate.storm_density_calculation(storms_runid, years, months_sh, 'sh', genesis=genesis, nyears = nyears, data_freq=data_freq)
            cube_model = (cube_nh + cube_sh)
            cube = cube_model - cube_obs
            cube_density[(runid, key)] = cube
            title = runid+'-'+runid_info['resol'][ir]+', '+str(years_model[0])+'-'+str(years_model[-1])+', '+algo
            if regionx != [] or regiony != []:
                if regionx != []:
                    ax.set_xlim([regionx[0], regionx[1]])
                if regiony != []:
                    ax.set_ylim([regiony[0], regiony[1]])
            storm_density_plot(cube, years, fig, ax, basin=None, title = title, plot_diff = True, genesis = genesis, bias = True)

    ax = fig.add_subplot(ny,nx,ny*nx,projection=ccrs.PlateCarree(central_longitude=central_longitude))
    title = 'OBS, '+str(years_obs[0])+'-'+str(years_obs[-1])
    if regionx != [] or regiony != []:
        if regionx != []:
            ax.set_xlim([regionx[0], regionx[1]])
        if regiony != []:
            ax.set_ylim([regiony[0], regiony[1]])
    storm_density_plot(cube_obs, years_obs, fig, ax, basin=None, title = title, genesis = genesis)
    
    fig.subplots_adjust(bottom=0.08, left=0.05, right=0.95, top=0.92, wspace=0.0, hspace=0.26)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    if genesis:
        fout = models+'-'+algos+plot_name_region+'_genesis_bias.png'
    else:
        fout = models+'-'+algos+plot_name_region+'_density_bias.png'
    figname = os.path.join(plot_dir, fout)

    plt.savefig(figname)
    plt.close()

def do_mslp_wspeed(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, keys = [], years_obs=[]):

    fig = plt.figure(figsize=[9,10], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = 2; ny = 2

    for ib, basin in enumerate(basins):
        ax = fig.add_subplot(ny,nx,ib+1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            x_fit_obs, y_fit_obs = tc_intensity.calculate_fit_obs(obs_storms, basin = basin)
            tc_intensity.plot_obs_mslp_wspeed(obs_storms, ax, x_fit=x_fit_obs, y_fit = y_fit_obs, basin = basin, linewidth = 1)

        for ir, runid in enumerate(runid_info['model']):
            for algo in runid_info['algo_type'][ir:ir+1]:
                key = keys[ir]
                #print ('mslp_wspeed ',runid, algo)
                storms_runid = storms[(runid, key)]
                x_fit, y_fit = tc_intensity.calculate_fit_model(storms_runid, use_925 = False, basin = basin)
                #print ('x_fit, y_fit ',x_fit, y_fit)
                plt.title('')
                tc_intensity.plot_model_mslp_wspeed(storms_runid, i, ax, x_fit=x_fit, y_fit=y_fit, label=runid+'-'+algo, use_925 = False, paired = False, basin = basin)
                i+=1

                ax.legend(loc=3, fontsize='x-small')

        tc_intensity.setup_mslp_windspeed_plot(ax, basin = basin, fig_no = ib+1)    

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_mslp_windspeed.png')
    plt.savefig(figname)

def do_size(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False):

    fig = plt.figure(figsize=[9,10], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = 2; ny = 2

    for ib, basin in enumerate(basins):
        i = 0
        ax = fig.add_subplot(ny,nx,ib+1)
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            x_fit_obs, y_fit_obs = tc_intensity.calculate_fit_obs(obs_storms, basin = basin)
            tc_intensity.plot_obs_mslp_wspeed(obs_storms, ax, x_fit=x_fit_obs, y_fit = y_fit_obs, basin = basin, linewidth = 1)

        for runid in runid_info['model']:
            for algo in runid_info['algo_type']:
                storms_runid = storms[(runid, algo)]
                #x_fit, y_fit = tc_intensity.calculate_fit_model(storms_runid, use_925 = False, basin = basin)
                #print ('x_fit, y_fit ',x_fit, y_fit)
                plt.title('')
                tc_intensity.plot_model_mslp_size(storms_runid, i, ax, label=runid+'-'+algo, use_925 = False, paired = False, basin = basin)
                ax.legend(loc=3, fontsize='x-small')
                i+=1


        #tc_intensity.setup_mslp_windspeed_plot(ax, basin = basin, fig_no = ib+1)    

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_mslp_size.png')
    plt.savefig(figname)

def plot_radius_max_wind_vmax(storms, runid_info, plot_dir, do_radmax = True, basins = ['nh'], keys = []):
    '''
    Keywords
    do_radmax - True - use the radius of maximum wind
              - False - use the outer 8m/s contour
    '''
    fig = plt.figure(figsize=[10,10])
    cols = ['blue', 'orange', 'red', 'purple', 'cyan']
    fontsize = 18
    fontdict = {'fontsize': fontsize, 'fontweight' : 8}
    fontdict_label = {'fontsize': fontsize, 'fontweight' : 6}
    for ib, basin in enumerate(basins):
        ax = fig.add_subplot(2, 2, ib+1)
        for ir, runid in enumerate(runid_info['model']):
            print('runid ',runid)
            for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
                key = keys[ir]
                title = runid+', '+algo
                storms_runid = storms[(runid, key)]
                vm = []
                rad = []
                for storm in storms_runid:
                    if tc_assess_code._storm_vmax_in_basin(storm, basin):
                        try:
                            vmax = storm.obs_at_vmax().vmax
                            if do_radmax:
                                radius = storm.radius_max_wind()
                            else:
                                radius = storm.obs_at_vmax().radius_8
                            rad.append(radius)
                            vm.append(vmax)
                        except:
                            #rad.append(0.0)
                            pass
                resol = runid+', '+runid_info['resol'][ir]
                plt.scatter(rad, vm, label = runid+' '+resol.upper(), color = cols[ir])
        ax.set_xlim([0, 12])
        ax.set_ylim([0, 80])
        if do_radmax:
            plt.suptitle('Relation of radius_max and Vmax (925hPa)', fontdict = fontdict)
            plt.title(basin.upper(), fontdict = fontdict)
            plt.xlabel('Radius of max wind ($^{o}$)', fontdict = fontdict_label)
        else:
            plt.suptitle('Relation of radius_8 and Vmax (925hPa)', fontdict = fontdict)
            plt.title(basin.upper(), fontdict = fontdict)
            plt.xlabel('Radius of 8$ms^{-1}$ wind ($^{o}$)', fontdict = fontdict_label)
        plt.ylabel('Max. wind at 925 hPa', fontdict = fontdict_label)
        plt.legend(loc = 'upper left', fontsize = 8)

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    if do_radmax:
        figname = os.path.join(plot_dir, models+'-'+algos+'_size_radmax_vmax.png')
    else:
        figname = os.path.join(plot_dir, models+'-'+algos+'_size_rad8_vmax.png')

    plt.savefig(figname)

def plot_interannual_variability(years, runid_info, plot_dir, storedir, months_nh, months_sh, basins = ['na'], keys = [], coupled = True):
    '''
    Calculate interannual variability
    '''
    years_hurdat = np.arange(1980, 2015)
    colour_runid = ['red', 'blue', 'green', 'cyan', 'orange', 'purple']
    colour_obs = ['black', 'gray']
    nensemble = 1
    annual_count_obs = {}
    annual_count_models = {}
    annual_corr_models = {}
    mean_count = {}
    std_count = {}
    mean_count_cat = {}
    std_count_cat = {}

    obs_ds = 'hurdat2'
    obs_datasets = [obs_ds]
    plot_hurdat = True
    annual_count_obs[obs_ds] = {}
    annual_count_obs['years'] = years_hurdat
    obs_tc_pickle_dir = storedir.format(obs_ds, '')
    if not os.path.exists(obs_tc_pickle_dir):
        os.makedirs(obs_tc_pickle_dir)
    obs_tc_pickle_file = os.path.join(obs_tc_pickle_dir, obs_ds+'_'+str(years_hurdat[0])+'-'+str(years_hurdat[-1])+'.pkl')
    if not os.path.exists(obs_tc_pickle_file):
        print('calculating obs data')
        obs_storms = ts_obs.load_data.load_data()
        tc_interannual.calculate_hurdat_frequency(obs_storms, annual_count_obs, obs_ds, BASINS, months_nh, months_sh)
        with open(obs_tc_pickle_file, 'wb') as fh:
            pickle.dump(annual_count_obs[obs_ds], fh)
    else:
        with open(obs_tc_pickle_file, 'rb') as fh:
            annual_count_obs[obs_ds] = pickle.load(fh)

    #print('obs var ',annual_count_obs['years'], annual_count_obs[obs_ds][('tc','na')])

    runids = []; runkeys = []; resols = []
    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            key = keys[ir]
            ens = runid_info['ens'][ir]
            runids.append(runid)
            runkey = runid+'_'+key
            for basin in BASINS:
                std_count[(runkey, basin)] = 0.0
            runkeys.append(runkey)
            resols.append(runid_info['resol'][ir])
            annual_count_models[runkey] = {}
            for basin in BASINS:
                annual_count_models[runkey][basin] = []
                annual_count_models[runkey]['ace', basin] = []
                for icat in range(6):
                    annual_count_models[runkey][basin, icat] = []
                    annual_count_models[runkey]['ace', basin, icat] = []
            if ens == '':
                outdir = storedir.format(runid, algo)
            else:
                outdir = storedir.format(runid+'/'+ens, algo)
            fname_format = os.path.join(outdir, 'tc_count_ace_{}_{}.nc')
            tc_assess_code.calc_interann.read_model_count_ace(fname_format, years, runkey, annual_count_models, outdir)
            #print('model annual ',runid, algo, annual_count_models[runkey]['nh'], annual_count_models[runkey][('years', 'nh')])

            tc_assess_code.calc_interann.calc_tc_frequency_mean_std(runid, runkey, BASINS, annual_count_obs, annual_count_models, mean_count, std_count)

    if not coupled:
    # calculate correlations
        correlations = tc_interannual._basin_correlations(runids, runkeys, BASINS, plot_dir, annual_count_models, annual_count_obs, annual_corr_models)

        tc_interannual._create_correlation_plot(runids, runkeys, annual_corr_models, BASINS, BASIN_NAME, BASINS_LIST, plot_dir, colour_runid, resols)

        tc_interannual._create_correlation_plot(runids, runkeys, annual_corr_models, BASINS, BASIN_NAME, BASINS_LIST, plot_dir, colour_runid, resols, ace=True)

    #if not do_metric:
    for do_ace in [False, True]:
        tc_interannual._create_interann_plot(runids, runkeys, resols, annual_count_obs, \
                                                     annual_count_models, plot_dir, BASINS, obs_datasets, colour_runid, colour_obs, years, std_count, plot_hurdat, plot_all_years=True, show_plots=False, use_model_years=True, ace=do_ace)
    
        tc_interannual._create_basin_plot(runids, runkeys, mean_count, annual_count_obs, BASINS, BASIN_NAME, BASINS_LIST, plot_dir, colour_runid, resols)

def do_tc_intensity(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, variable = '10mwind', keys = [], years_obs=[]):
    '''
    Plot the proportion of TCs in each intensity category
    '''

    colors = ['#ff6666', '#ffcc99', '#cc9966', '#cc6666', '#66b3ff', '#c2c2f0']
    labels = ['0','1','2','3','4','5']
    fontsize = 10
    fontdict = {'fontsize': fontsize, 'fontweight' : 4}
    fig = plt.figure(figsize=[9,10], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = len(basins)+1
        ny = np.amax([len(runid_info['model']), len(runid_info['algo_type'])])+1

    radius_ref = 0.5
    storm_cats = {}
    if variable == '10mwind':
        vmax_range = np.arange(0, 100, 1)
        plot_cat = True
        limits_new = wind_limits_new
    elif variable == 'slp' or variable == 'mslp':
        vmax_range = np.arange(890, 1025, 1)
        plot_cat = False
        limits_new = mslp_limits_new

    for ib, basin in enumerate(basins):
        ax = fig.add_subplot(ny,nx,nx*(ny-1)+ib+1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            obs_histogram, edges = calculate_histogram_maxvmax(obs_storms_tsandup, vmax_range, scaling = knots_to_ms, basin = basin, variable = variable, obs = True)
            storm_cats[('obs',basin)] = []
            for icat in range(len(limits_new)-1):
                limits = np.where((limits_new[icat] <= edges)&(edges < limits_new[icat+1]))[0]
                cat = np.sum(obs_histogram[limits])
                storm_cats[('obs', basin)].append(cat)
            if variable == 'mslp' or variable == 'slp':
                storm_cats[('obs',basin)].reverse()
            #print ('storm_cats obs ',variable, basin, storm_cats[('obs', basin)])
            wedges, texts, autotexts = ax.pie(storm_cats[('obs', basin)], startangle=90, frame=True, autopct='%1.0f', labeldistance = 1.0, pctdistance = 0.7, textprops = dict(color='black'), radius = radius_ref, colors=colors)
            ax.axes.get_yaxis().set_visible(False)
            ax.axes.get_xaxis().set_visible(False)
            plt.axis('equal')
            plt.axis('off')
            title = 'OBS, '+basin.upper()
            plt.title(title, loc = 'left', fontdict = fontdict)

        for ir, runid in enumerate(runid_info['model']):
            for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
                key = keys[ir]
                title = runid+', '+algo+', '+basin.upper()
                storm_cats[(runid, algo, basin)] = []
                storms_runid = storms[(runid, key)]
                data_freq = storms[(runid, key, 'freq')]
                hist, edges = calculate_histogram_maxvmax(storms_runid, vmax_range, basin = basin, variable = variable, data_freq=data_freq)
                for icat in range(len(limits_new)-1):
                    lower_limit = limits_new[icat]
                    limits = np.where((lower_limit <= edges)&(edges < limits_new[icat+1]))[0]
                    cat = np.sum(hist[limits])
                    storm_cats[(runid, algo, basin)].append(cat)
                if variable == 'mslp' or variable == 'slp':
                    storm_cats[(runid, algo, basin)].reverse()
                #print('model cats ',runid, algo, variable, basin, storm_cats[(runid, algo, basin)])
                #print('hist ',hist)
                ax = fig.add_subplot(ny,nx,ir*nx+ia*nx+ib+1)
                wedges, texts, autotexts = ax.pie(storm_cats[(runid, algo, basin)], startangle=90, frame=True, autopct='%1.0f', labeldistance = 1.0, pctdistance = 0.7, textprops = dict(color='black'), radius = radius_ref, colors=colors)
                ax.axes.get_yaxis().set_visible(False)
                ax.axes.get_xaxis().set_visible(False)
                plt.axis('equal')
                plt.axis('off')
                plt.title(title, loc = 'left', fontdict = fontdict)
    
    plt.suptitle('Intensity based on '+variable)

    ax = fig.add_subplot(ny,nx,nx*ny)
    plot_key = True
    if plot_key:
        ax.legend(wedges, labels,
          title="TC Cat %",
          loc="center left",
                  bbox_to_anchor=(0.7, 0, 0.25, 1.0), fontsize = fontsize-2, borderaxespad=0)
        ax.axes.get_yaxis().set_visible(False)
        ax.axes.get_xaxis().set_visible(False)
        plt.axis('equal')
        plt.axis('off')

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_pie_intensity_'+variable+'.png')
    plt.savefig(figname)

def do_seasonal_cycle(storms, years, runid_info, plot_dir, keys = [], years_obs=[]):
    months = range(1,13)
    monthly_variability_with_obs(storms, years, months, runid_info, keys = keys, years_obs=years_obs)

    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_seasonal_cycle.png')
    plt.savefig(figname)

def filter_obs(storms, storm_types = ['SS', 'TS', 'HU', 'MH']):
    storms_filtered = []
    for storm in storms:
        if (storm.type() in storm_types):
            storms_filtered.append(storm)
    return storms_filtered

def filter_model(storms, tc_thres_mslp = False, tc_thres_10mwind = False):
    if not tc_thres_mslp and not tc_thres_10mwind:
        return []

    storms_filtered = []
    if tc_thres_mslp:
        for storm in storms:
            if (storm.obs_at_vmax().mslp <= tropical_storm_threshold['mslp']):
                storms_filtered.append(storm)
    elif tc_thres_10mwind:
        for storm in storms:
            if (storm.obs_at_vmax().extras['w10m'] >= tropical_storm_threshold['10mwind']):
                storms_filtered.append(storm)

    return storms_filtered

def filter_storms_threshold(storms, runid, algo, vort_threshold = 6.0):
    '''
    Filter the storms based on a threshold of a variable
    '''
    storms_filtered = []
    storms_to_filter = storms[(runid, algo)]
    for storm in storms_to_filter:
        vort_max = storm.obs_at_vmax().t63_850
        if vort_max == 1.0:
            if storm.obs_at_genesis().lat > 0:
                vort_max = storm.obs_at_vmax().vort
            else:
                vort_max = storm.obs_at_vmax().vort_min*-1.0
            if np.abs(vort_max) < 1e-1:
                vort_max *= 1.0e5
        if  np.abs(vort_max) >= vort_threshold:
            storms_filtered.append(storm)

    storms[(runid, algo)] = storms_filtered

def filter_time(storms, years, months):
    storms_filtered = []
    for year in years:
        for storm in _storms_in_time_range(storms, year, months):
            storms_filtered.append(storm)
    return storms_filtered

def filter_space(storms, region):
    storms_filtered = []
    for storm in storms:
        if _storm_in_basin(storm, region):
            if storm.obs_at_genesis().lat > -90. and storm.obs_at_genesis().lat < 90.:
                storms_filtered.append(storm)
    return storms_filtered

def do_piechart(storms, years, runid_info, plot_dir, tc_thres_mslp = False, tc_thres_10mwind = False, keys = [], years_obs=[]):
    basins = ['nh', 'na', 'wp', 'ep', 'ni', 'sh']
    months = {}
    months['nh'] = [5,6,7,8,9,10,11]
    months['sh'] = [10,11,12,1,2,3,4,5]

    fig = plt.figure(figsize=[8,6], dpi=150)
    #nx = len(runid_info['algo_type'])
    nn = len(runid_info['model'])+1
    nx = ((nn+1) // 2)
    ny = ((nn+1) // nx)

    fontsize = 12
    if ny > 2:
        fontsize = 10

    for ir, runid in enumerate(runid_info['model']):
        for ia, algo in enumerate(runid_info['algo_type'][ir:ir+1]):
            #print ('ir, ia, plt ',ir, ia, nx,ny,ir+1)
            key = keys[ir]
            ax = fig.add_subplot(ny,nx,ir+1)

            if tc_thres_mslp:
                storms_filter = filter_model(storms[(runid, key)], tc_thres_mslp)
                file_ending = '_mslpfilter'
            elif tc_thres_10mwind:
                storms_filter = filter_model(storms[(runid, key)], tc_thres_10mwind)
                file_ending = '_10mwindfilter'
            else:
                storms_filter = storms[(runid, key)]
                file_ending = ''

            mean_count = {}; obs_count = {}
            for basin in basins:
                hemi = BASINS[basin][0]
                months_basin = months[hemi]
                nyears = storms[(runid, key, 'nyears')]
                mean_count[basin] = tc_assess_code.get_monthly_mean_count(storms_filter, years, months_basin, basin, annual = True, nyears = nyears)
                obs_storms = ts_obs.load_data.load_data(basin=basin)
                obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
                obs_count[basin] = tc_assess_code.get_monthly_mean_count(obs_storms_tsandup, years_obs, months_basin, basin, annual = True)
                #print('mean count obs ',basin, obs_count[basin])

            ct = ''
            plot_piechart(mean_count, basins, ax, fig, obs_count, title = runid+'\n'+key, plot_dir = './', plot_key = False, fontsize = fontsize)

    ax = fig.add_subplot(ny,nx,ny*nx)
    ct = ''
    plot_piechart(obs_count, basins, ax, fig, obs_count, title = 'Obs', plot_dir = './', plot_key = True)

    title = 'TC NH piechart '
    if file_ending != '':
        title += file_ending[1:]
    plt.suptitle(title)

    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_piechart_basinfreq'+file_ending+'.png')
    plt.savefig(figname)

def do_intensification_pdf(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, variable = '10mwind', keys = [], years_obs = []):

    fig = plt.figure(figsize=[10,9], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = 2; ny = 2

    ymax = 0.2
    if variable == '10mwind':
        vmax_range = np.arange(0, 39, 3)
        plot_cat = False
    elif variable == 'slp' or variable == 'mslp':
        vmax_range = np.arange(-80, 0, 4)
        plot_cat = False
    elif variable == 'vmax':
        vmax_range = np.arange(0, 50, 2)
        plot_cat = False
        do_obs = False
    elif variable == 'rvT63':
        vmax_range = np.arange(0, 1, 0.01)
        plot_cat = False
        do_obs = False
        ymax = 0.2

    for ib, basin in enumerate(basins):
        title = basin.upper()+', '+variable
        if variable == 'vmax':
            title += ' (925hPa wind)'
        ax = fig.add_subplot(ny,nx,ib+1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            obs_histogram, edges = calculate_histogram_maxintens(obs_storms_tsandup, vmax_range, scaling = knots_to_ms, basin = basin, variable = variable, obs = True)
            plot_latitude_histogram_nofig(obs_histogram, 'OBS', '', ax, edges = edges, x_label = 'Max intensification rate (/day)', title=title, plot_cat = plot_cat, ymax = ymax, color = 'black')
            #stop

        for ir, runid in enumerate(runid_info['model']):
            for algo in runid_info['algo_type'][ir:ir+1]:
                key = keys[ir]
                storms_runid = storms[(runid, key)]
                data_freq = storms[(runid, key, 'freq')]
                hist, edges = calculate_histogram_maxintens(storms_runid, vmax_range, basin = basin, variable = variable, data_freq=data_freq)
                plot_latitude_histogram_nofig(hist, runid, algo, ax, edges = edges, x_label = 'Max intensification rate (/day)', title=title, plot_cat = plot_cat, ymax = ymax, color=colournames[ir])

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_'+variable+'_intensification_pdf.png')
    plt.savefig(figname)

def do_intensification_24h_pdf(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, variable = '10mwind', keys = [], years_obs = []):

    fig = plt.figure(figsize=[10,9], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = 2; ny = 2

    ymax = 0.2
    if variable == '10mwind':
        vmax_range = np.arange(-30, 32, 3)
        plot_cat = False
    elif variable == 'slp' or variable == 'mslp':
        vmax_range = np.arange(-40, 40, 4)
        plot_cat = False
    elif variable == 'vmax':
        vmax_range = np.arange(-40, 40, 4)
        plot_cat = False
        do_obs = False
    elif variable == 'rvT63':
        vmax_range = np.arange(-1, 1, 0.01)
        plot_cat = False
        do_obs = False
        ymax = 0.2

    for ib, basin in enumerate(basins):
        title = basin.upper()+', '+variable
        if variable == 'vmax':
            title += ' (925hPa wind)'
        ax = fig.add_subplot(ny,nx,ib+1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin) 
            obs_histogram, edges = calculate_histogram_maxintens_24h(obs_storms, vmax_range, scaling = knots_to_ms, basin = basin, variable = variable, obs = True)
            plot_latitude_histogram_nofig(obs_histogram, 'OBS', '', ax, edges = edges, x_label = 'Max intensification rate (/24h)', title=title, plot_cat = plot_cat, ymax = ymax, color = 'black')
            #print('obs, var, >15 ',basin, variable, np.sum(obs_histogram[15:]))
            obs_storms_ibtracs = ts_obs.track_ibtracs_netcdf.load_cmor(IBTRACS_FILE)
            #obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            obs_histogram_ib, edges_ib = calculate_histogram_maxintens_24h(obs_storms_ibtracs, vmax_range, scaling = 1.0, basin = basin, variable = variable, obs = True)
            plot_latitude_histogram_nofig(obs_histogram_ib, 'OBS_IB', '', ax, edges = edges_ib, x_label = 'Max intensification rate (/24h)', title=title, plot_cat = plot_cat, ymax = ymax, color = 'gray')
            #print('ibtracs, var, >15 ',basin, variable, np.sum(obs_histogram_ib[15:]))

        for ir, runid in enumerate(runid_info['model']):
            for algo in runid_info['algo_type'][ir:ir+1]:
                key = keys[ir]
                storms_runid = storms[(runid, key)]
                hist, edges = calculate_histogram_maxintens_24h(storms_runid, vmax_range, basin = basin, variable = variable)
                plot_latitude_histogram_nofig(hist, runid, algo, ax, edges = edges, x_label = 'Max intensification rate (/24h)', title=title, plot_cat = plot_cat, ymax = ymax, color = colournames[ir])
                #print('model, var, >15 ',runid, basin, variable, np.sum(hist[15:]))

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_'+variable+'_intensification_24h_pdf.png')
    plt.savefig(figname)

def do_1D_pdf(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, variable = '10mwind', keys = []):

    fig = plt.figure(figsize=[10,9], dpi=100)
    if len(basins) == 1:
        nx = 1; ny = 2
    elif len(basins) <= 4:
        nx = 2; ny = 2

    ymax = 0.1
    if variable == '10mwind':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = True
    elif variable == 'slp' or variable == 'mslp':
        vmax_range = np.arange(890, 1025, 5)
        plot_cat = False
    elif variable == 'vmax':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = False
        do_obs = False
    elif variable == 'latitude':
        vmax_range = np.arange(2, 45, 2)
        plot_cat = False
    elif variable == 'lifetime':
        vmax_range = np.arange(0, 26, 2)
        plot_cat = False
        ymax = 0.2
    elif variable == 'rvT63':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = False
        do_obs = False
        ymax = 0.2

    for ib, basin in enumerate(basins):
        title = basin.upper()+', '+variable
        if variable == 'vmax':
            title += ' (925hPa wind)'
        ax = fig.add_subplot(ny,nx,ib+1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            obs_histogram, edges = calculate_histogram_maxvmax(obs_storms_tsandup, vmax_range, scaling = knots_to_ms, basin = basin, variable = variable, obs = True)
            plot_latitude_histogram_nofig(obs_histogram, 'OBS', '', ax, edges = edges, x_label = variable, title=title, plot_cat = plot_cat, ymax = ymax, color = 'black')

        for ir, runid in enumerate(runid_info['model']):
            for algo in runid_info['algo_type'][ir:ir+1]:
                key = keys[ir]
                storms_runid = storms[(runid, key)]
                data_freq = storms[(runid, key, 'freq')]
                hist, edges = calculate_histogram_maxvmax(storms_runid, vmax_range, basin = basin, variable = variable, data_freq=data_freq)
                plot_latitude_histogram_nofig(hist, runid, algo, ax, edges = edges, x_label = variable, title=title, plot_cat = plot_cat, ymax = ymax, color=colournames[ir])

    fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
    models = '_'.join(runid_info['model'])
    algos = '_'.join(runid_info['algo_type'])
    figname = os.path.join(plot_dir, models+'-'+algos+'_'+variable+'_pdf.png')
    plt.savefig(figname)

def do_2D_pdf(storms, runid_info, plot_dir, basins = ['nh'], do_obs = False, variable_x = 'mslp', variable_y = 'latitude', keys = []):

    nplots = len(runid_info['model'])+1
    #nplots = 4
    nx = int((nplots+1)/2)+1
    ny = int((nplots+1)/nx)+1
    #nx = 2
    #ny = 2

    ymax = 0.1
    if variable_x == '10mwind':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = True
    elif variable_x == 'slp' or variable_x == 'mslp':
        vmax_range = np.arange(890, 1025, 5)
        plot_cat = False
    elif variable_x == 'vmax':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = False
        do_obs = False
    elif variable_x == 'latitude':
        vmax_range = np.arange(2, 45, 2)
        plot_cat = False
    elif variable == 'lifetime':
        vmax_range_x = np.arange(0, 26, 2)
        plot_cat = False
        ymax = 0.2
    elif variable_x == 'rvT63':
        vmax_range = np.arange(0, 100, 5)
        plot_cat = False
        do_obs = False
        ymax = 0.2

    for ib, basin in enumerate(basins):
        fig = plt.figure(figsize=(12,8),dpi=100)

        title = basin.upper()+', '+variable_x
        if variable_x == 'vmax':
            title += ' (925hPa wind)'
        ax = fig.add_subplot(ny,nx,1)
        i = 0
        if do_obs:
            obs_storms = ts_obs.load_data.load_data(basin=basin)       
            obs_storms_tsandup = filter_obs(obs_storms, storm_types = ['TS', 'HU', 'MH'])
            max_lat_histogram_obs, xedges, yedges = tc_intensity.calculate_histogram_maxlat(obs_storms_tsandup, basin = basin)
            #print('xedges ',xedges)
            tc_intensity.plot_2d_histogram(max_lat_histogram_obs, 'OBS', 'OBS', 'OBS', plot_dir, fig, ax, nx, ny, xedges = xedges, yedges = yedges, xlabel = 'Min. MSLP', title = 'Max_intensity-lifetime', min_val = 0.0001)

        for ir, runid in enumerate(runid_info['model']):
            for algo in runid_info['algo_type'][ir:ir+1]:
                key = keys[ir]
                ax = fig.add_subplot(ny,nx,ir+2)
                resol = runid
                model_resol = runid
                storms_runid = storms[(runid, key)]
                max_lat_histogram, xedges, yedges = tc_intensity.calculate_histogram_maxlat(storms_runid, basin = basin)
                tc_intensity.plot_2d_histogram(max_lat_histogram, runid, resol, model_resol, plot_dir, fig, ax, nx, ny, xedges = xedges, yedges = yedges, xlabel = 'Min. MSLP', title = 'Max_intensity-lifetime, '+basin, min_val = 0.0001)

        fig.subplots_adjust(bottom=0.08, left=0.08, right=0.95, top=0.92, wspace=0.27, hspace=0.27)
        models = '_'.join(runid_info['model'])
        algos = '_'.join(runid_info['algo_type'])
        figname = os.path.join(plot_dir, models+'-'+algos+'_'+basin+'_2d_pdf.png')
        plt.savefig(figname)

def plot_etc_properties(runid_info, storms, years, storedir, do_meanstate, do_variability, years_obs, plot_dir, keys, months_nh, months_sh, coupled=False):
    '''
    Plot all TC metrics/diagnostics
    '''
    #do_size(storms, runid_info, plot_dir, basins = ['na', 'wp', 'ep', 'si'])
    months_nh = range(1,13); months_sh = months_nh
    years_plot = years[:2]
    years_plot_obs = years_obs[:2]

    time_var_start = time.perf_counter()
    if do_variability:
        plot_interannual_variability(years, runid_info, plot_dir, storedir, months_nh, months_sh, keys = keys, coupled=coupled)
        plt.close('all')
    time_var = time.perf_counter()
    print('Variability took ',(time_var - time_var_start),' seconds')

    if do_meanstate:
        #do_seasonal_cycle(storms, years, runid_info, plot_dir, keys = keys, years_obs=years_obs)
        #do_storm_tracks(storms, years_plot, runid_info, plot_dir, months_nh, months_sh, warmcore = False, keys = keys, regionx = [100, 180], regiony = [0, 50.])
        do_storm_tracks(storms, years_plot, runid_info, plot_dir, months_nh, months_sh, warmcore = False, keys = keys, years_obs=years_plot_obs)

        time_tracks = time.perf_counter()
        print('Tracks took ',(time_tracks - time_var),' seconds')

        do_track_density_plot(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = False, keys = keys, regionx = [100, 180], regiony = [0, 50.], years_obs=years_obs, storedir=storedir)
        do_track_density_plot(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = False, keys = keys, years_obs=years_obs, storedir=storedir)
        do_track_density_plot(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = True, keys = keys, regionx = [100, 180], regiony = [0, 50.], years_obs=years_obs, storedir=storedir)
        do_track_density_plot(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = True, keys = keys, years_obs=years_obs, storedir=storedir)

        time_density = time.perf_counter()
        print('Density took ',(time_density - time_tracks),' seconds')

        do_piechart(storms, years, runid_info, plot_dir, keys = keys, years_obs=years_obs)
        do_piechart(storms, years, runid_info, plot_dir, tc_thres_mslp = True, keys = keys, years_obs=years_obs)
        do_piechart(storms, years, runid_info, plot_dir, tc_thres_10mwind = True, keys = keys, years_obs=years_obs)
        time_pie = time.perf_counter()
        print('Pie took ',(time_pie - time_density),' seconds')

        for variable in ['10mwind', 'mslp']:
            do_intensification_24h_pdf(storms, runid_info, plot_dir, basins = ['na', 'wp', 'ep'], do_obs = True, variable = variable, keys = keys, years_obs=years_obs)
            do_intensification_pdf(storms, runid_info, plot_dir, basins = ['na', 'wp', 'ep'], do_obs = True, variable = variable, keys = keys, years_obs=years_obs)
        time_intense = time.perf_counter()
        print('Intensify took ',(time_intense - time_pie),' seconds')

        do_mslp_wspeed(storms, runid_info, plot_dir, basins = ['na', 'wp', 'ep'], do_obs = True, keys = keys, years_obs=years_obs)

        do_tc_intensity(storms, runid_info, plot_dir, basins = ['na','wp','ep'], do_obs = True, variable = '10mwind', keys = keys, years_obs=years_obs)

        do_tc_intensity(storms, runid_info, plot_dir, basins = ['na','wp','ep'], do_obs = True, variable = 'mslp', keys = keys, years_obs=years_obs)

        time_intense2 = time.perf_counter()
        print('Intense took ',(time_intense2 - time_intense),' seconds')
        for variable in ['10mwind', 'mslp', 'vmax', 'lifetime', 'latitude', 'rvT63']:
            do_1D_pdf(storms, runid_info, plot_dir, basins = ['na', 'wp', 'ep', 'si'], do_obs = True, variable = variable, keys = keys)
        time_pdf = time.perf_counter()
        print('Pdf took ',(time_intense2 - time_intense),' seconds')

        do_2D_pdf(storms, runid_info, plot_dir, basins = ['wp', 'na', 'ep'], do_obs = True, variable_x = 'mslp', variable_y = 'latitude', keys = keys)

        #plot_radius_max_wind_vmax(storms, runid_info, plot_dir, do_radmax = False, basins = ['wp', 'na', 'ep'], keys = keys)
        #plot_radius_max_wind_vmax(storms, runid_info, plot_dir, do_radmax = True, basins = ['wp', 'na', 'ep'], keys = keys)

        #do_storm_tracks(storms, years_plot, runid_info, plot_dir, months_nh, months_sh, warmcore = True)

        do_track_density_bias(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = False, keys = keys, years_obs=years_obs)
        do_track_density_bias(runid_info, storms, years, plot_dir, months_nh, months_sh, genesis = True, keys = keys, years_obs=years_obs)

        time_mean1 = time.perf_counter()
        print('Rest took ',(time_mean1 - time_pdf),' seconds')
    plt.show()

    time_mean = time.perf_counter()
    print('Mean state took ',(time_mean - time_var),' seconds')

    months_nh = range(1,13); months_sh = months_nh
    years_plot = years[:6]
    print ('do storm tracks')
    #try:
    #    do_storm_tracks(storms, years_plot, runid_info, plot_dir, months_nh, months_sh)
    #except:
    #    print('plot storm tracks failed')
    plt.show()

