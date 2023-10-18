""" 
Load function to read in TempestExtremes output from inline tracking. 
Written for the CMORised netcdf files
Should work for both an individual month and all months 

Probably need to replace netcdftime with cftime for forwards compatibility with iris2.1

"""
import os.path, sys
import collections
import datetime
import netCDF4
import cftime
import numpy as np
from itertools import groupby
from operator import itemgetter

""" 
Store model storm observations as a named tuple (this enables access to data by its name instead of position index) 
"""

Observation = collections.namedtuple('Observation', ['date', 'lat', 'lon', 'vort', 'vort_min', 'vmax', 
                                                     'mslp', 't63_850', 't63_700', 't63_600', 't63_500', 't63_250', 'zg_avg_925', 'zg_max_925', 'zg_min_925', 'zg_avg_600', 'zg_max_600', 'zg_min_600', 'zg_avg_500', 'zg_max_500', 'zg_min_500', 'zg_avg_300', 'zg_max_300', 'zg_min_300', 'zg_avg_250', 'zg_max_250', 'zg_min_250', 'radius_8', 'ace', 'ace_psl', 'ike', 'pdi', 'rprof', 'rprof_slp', 'extras'])

class Observation(Observation):
    """  Represents a single observation of a model tropical storm. """
    
    def six_hourly_timestep(self):
        """ Returns True if a storm record is taken at 00, 06, 12 or 18Z only """
        return self.date.hour in (0, 6, 12, 18) and self.date.minute == 0 and self.date.second == 0
    
    def add_to_axes(self, ax):
        """ Instructions on how to plot a model tropical storm observation """
        ax.plot(self.lon, self.lat)
        

class Storm(object):
    def __init__(self, snbr, obs, extras=None):
        """ Stores information about the model storm (such as its storm number) and corresponding 
        observations. Any additional information for a storm should be stored in extras as a dictionary. """
        self.snbr = snbr
        self.obs = obs
        
        # Create an empty dictionary if extras is undefined
        if extras is None:
            extras = {}
        self.extras = extras 

    @property
    def vmax(self):
        """ The maximum wind speed attained by the storm during its lifetime """
        return max(ob.vmax for ob in self.obs)
    
    @property
    def mslp_min(self):
        """ The minimum central pressure reached by the storm during its lifetime (set to 
        -999 if no records are available) """
        mslps = [ob.mslp for ob in self.obs if ob.mslp != 1e12]  
        if not mslps:
            mslps = [-999]
        return min(mslps)
    
    @property
    def vort_max(self):
        """ The maximum 850 hPa relative vorticity attained by the storm during its lifetime """
        if self.obs_at_genesis().lat > 0:
            return max(ob.vort for ob in self.obs)
        else:
            return min(ob.vort_min for ob in self.obs)
    
    @property
    def t63_850_max(self):
        """ The maximum T63 850 vorticity attained by the storm during its lifetime """
        if self.obs_at_genesis().lat > 0:
            return max(ob.t63_850 for ob in self.obs)
        else:
            return min(ob.t63_850 for ob in self.obs)
    
    @property
    def t63_250_max(self):
        """ The maximum T63 250 vorticity attained by the storm during its lifetime """
        return max(ob.t63_250 for ob in self.obs)
    
    @property
    def t63_850_250_diff_max(self):
        """ The maximum T63 850 - 250 vorticity attained by the storm during its lifetime """
        return max(ob.t63_850_250_diff for ob in self.obs)
    
    @property
    def warm_core(self):
        """ Use the zg values to determine if this storm has a upper warm core """
        delta_z_250 = []; delta_z_500 = []; delta_z_925 = []; z_diff_upper = []; z_diff_lower = []
        for ob in self.obs:
            #delta_z_250.append(ob.zg_max_250 - ob.zg_min_250)
            #delta_z_500.append(ob.zg_max_500 - ob.zg_min_500)
            #delta_z_925.append(ob.zg_max_925 - ob.zg_min_925)
            z_diff_upper.append((ob.zg_max_600 - ob.zg_min_600) - (ob.zg_max_300 - ob.zg_min_300))
            z_diff_lower.append((ob.zg_max_925 - ob.zg_min_925) - (ob.zg_max_600 - ob.zg_min_600))
        #print('z_diff_upper ',z_diff_upper)
        pos_upper = np.where(np.asarray(z_diff_upper) > 0.)[0]
        pos_lower = np.where(np.asarray(z_diff_lower) > 0.)[0]
        # find maximum consecutive points
        max_consecutive = 0; consec_values = []
        for k, g in groupby(enumerate(pos_upper), lambda ix : ix[0] - ix[1]):
            consec = list(map(itemgetter(1),g))
            if len(consec) > max_consecutive: 
                max_consecutive = len(consec)
                consec_values = consec
        common_warmcore = sorted(list(set(consec_values).intersection(pos_lower)))
        max_wc = 0
        for ind in consec_values:
            max_wc = np.amax([max_wc, z_diff_upper[ind]])
        if len(common_warmcore) > 4 and max_wc > 30:
            print('warmcore ',max_consecutive, len(common_warmcore), z_diff_upper[consec_values[0]], max_wc)
            return True
        else:
            print('no warmcore ',max_consecutive, len(common_warmcore), max_wc)
            return False
    
    def __len__(self):
        """ The total number of observations for the storm """
        return len(self.obs)
    
    def nrecords(self):
        """ The total number of records/observations for the storm """
        return len(self)
    
    def number_in_season(self):
        """ Returns storm number of the storm (the number of the storm for that year and ensemble
        member number) """
        return self.snbr    
    
    def lifetime(self):
        """ The total length of time that the storm was active. This uses all observation
        points, no maximum wind speed threshold has been set """
        return max(ob.date for ob in self.obs)-min(ob.date for ob in self.obs)
        
    def genesis_date(self):
        """ The first observation date that a storm becomes active """
        #return min(ob.date for ob in self.obs)
        return self.obs_at_genesis().date
    
    def lysis_date(self):
        """ The final date that a storm was active """
        #return max(ob.date for ob in self.obs)
        return self.obs_at_lysis().date
    
    def ace_index(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_index = 0
        for ob in self.obs:
            if ob.six_hourly_timestep():
                ace_index += np.square(ob.extras['vmax_kts'])/10000.
        return round(ace_index, 2)
    
    def ace_index_no6hrcheck(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_index = 0
        for ob in self.obs:
            ace_index += np.square(ob.extras['vmax_kts'])/10000.
        return round(ace_index, 2)
    
    def ace_storm(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_storm = 0
        for ob in self.obs:
            if ob.six_hourly_timestep():
                ace_storm += ob.ace
        return round(ace_storm, 2)
    
    def ace_psl_storm(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_psl_storm = 0
        for ob in self.obs:
            if ob.six_hourly_timestep():
                ace_psl_storm += ob_psl.ace
        return round(ace_psl_storm, 2)
    
    def ike_storm(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ike_storm = 0
        for ob in self.obs:
            if ob.six_hourly_timestep():
                ike_storm += ob.ike
        return round(ike_storm, 2)
    
    def obs_at_vmax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.vmax)
    
    def obs_at_vortmax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.vort)
    
    def obs_at_ikemax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.ike)
    
#    def obs_at_coremax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
#	t63_diff = [ob.t63_diff for ob in self.obs]
#	peak = find_local_maximum(t63_diff)
        #return max(self.obs, key=lambda ob: ob.t63_diff)
#        return [ob for ob in self.obs][peak]
    
    def obs_at_mslpmin(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return min(self.obs, key=lambda ob: ob.mslp)
    
    def obs_at_genesis(self):
        """Returns the Observation instance for the first date that a storm becomes active """       
        for ob in self.obs:
            return ob
        else:
            raise ValueError('model storm was never born :-(')

    def obs_at_lysis(self):
        """Returns the Observation instance for the last date that a storm was active """    
        return [ob for ob in self.obs][-1]
    
    def radius_max_wind(self):
        """ Look through profile and find maximum """
        if self.obs_at_vmax().lat < 40.0:
            max_value = max(self.obs_at_vmax().rprof)
            max_pt = np.where(max_value == self.obs_at_vmax().rprof)[-1]
            #print('max_pt ',max_pt, max_value)
            radius_max_wind = 0.125 * float(max_pt+1)
            #print('radius_max_wind ',radius_max_wind)
            return radius_max_wind
        else:
            return 0.0

def load_cmor(fh, vort_variable = ''):
    '''
    Load cmor netcdf format of track files
    The input file should contain (at least):
       Dimensions:
          ntracks: number of storms in file
          record: total number of time points
          plev: number of pressure levels (with associated data)
       Variables:
          TRACK_ID(tracks): Storm number
          FIRST_PT(tracks): Index to first point in each storm
          NUM_PTS(tracks): The number of points in this storm
          index(record): storm track sequence number (index to this storm in the whole record diemsion)
          vortmean_T63(record): feature tracked variable
          lon(record): longitude of feature tracked variable
          lat(record): latitude of feature tracked variable
    '''
    # final variables needed
    # storm number snbr
    # date, lat, long, vort, vmax, mslp, T63[nlev], vmax_kts, w10m

    scaling_ms_knots = 1.944
    print('fh type', type(fh))
    fh_type = str(type(fh))
    if 'str' in fh_type:
        fh = [fh]
    
    # for each file in the file handle            
    for fname in fh:
        if not os.path.exists(fname):
            raise Exception('Input file does not exist '+fname)
        else:
            print('fname, track_inline ',fname)
            with netCDF4.Dataset(fname, 'r') as nc:
                try:
                    track_algorithm = nc.getncattr('algorithm')
                    if vort_variable == '':
                        if 'TRACK' in track_algorithm:
                            if nc.getncattr('algorithm_extra') == 'T63avg':
                                vort_variable = 'vortmean_T63'
                            else:
                                vort_variable = 'rvT42_850'
                        elif 'TempestExtremes' in track_algorithm:
                            vort_variable = 'psl'
                except:
                    vort_variable = 'rv850'

                # number of storms in the file
                ntracks = int(nc.dimensions['tracks'].size)
                try:
                    plev = int(nc.dimensions['plev'].size)
                except:
                    plev = 0

                variables = nc.variables
                # Loop through each storm, and create a class object containing the storm properties
                psl_var_options = ['psl', 'slp', 'slp_min', 'psl_min', 'air_pressure_at_mean_sea_level']
                for ps in psl_var_options:
                    if ps in variables:
                        psl = nc.variables[ps]
                        psl_var = ps

                psls = nc.variables[psl_var][:]
                divisor = 1.0
                try:
                    if psl.units == 'Pa':
                        divisor = 100.0
                except:
                    if np.median(psls) > 2000:
                        divisor = 100.0
                psl_scaling = 1.0 / divisor

                # read the time variable and convert to a more useful format
                time_var = nc.variables['time']
                try:
                    calendar = time_var.calendar
                except:
                    calendar = time_var.time_calendar
                dtime = cftime.num2date(time_var[:], time_var.units, calendar = calendar)

                # index to the first point of each track
                first_pts = nc.variables['FIRST_PT'][:]
                storm_lengths = nc.variables['NUM_PTS'][:]
                indices = nc.variables['index'][:]

                lat_var = ['lat', 'latitude']
                lon_var = ['lon', 'longitude']

                for lat in lat_var:
                    if lat in variables:
                        lats = nc.variables[lat][:]

                for lon in lon_var:
                    if lon in variables:
                        lons = nc.variables[lon][:]

                vort_var_options = [vort_variable, 'rvT63_max']
                vorts = np.ones(len(lons))
                for v in vort_var_options:
                    if v in variables:
                        vorts = nc.variables[vort_variable][:]

                vort_var_min_options = [vort_variable+'_1', 'rvT63_min']
                vorts_min = np.ones(len(lons))
                for v in vort_var_min_options:
                    if v in variables:
                        vorts_min = nc.variables[vort_variable][:]

                sfcw_options = ['sfcWind', 'sfcWind_max', 'wind_speed']
                sfcWinds = np.ones(len(lons))
                sfcWinds_in_data = False
                for sf in sfcw_options:
                    if sf in variables:
                        sfcWinds = nc.variables[sf][:]
                        sfcWinds_in_data = True

                ws925_options = ['ws925', 'wind925', 'wind_max_925']
                vmaxs = np.ones(len(lons))
                vmax_in_data = False
                for ws in ws925_options:
                    if ws in variables:
                        vmaxs = nc.variables[ws][:]
                        vmax_in_data = True
                if not vmax_in_data and sfcWinds_in_data:
                    vmaxs = sfcWinds
                
                # number of pressure levels for variables
                if plev >= 5:
                    rvT63_850 = nc.variables['rvT63_850'][:]
                    rvT63_700 = nc.variables['rvT63_700'][:]
                    rvT63_600 = nc.variables['rvT63_600'][:]
                    rvT63_500 = nc.variables['rvT63_500'][:]
                    rvT63_250 = nc.variables['rvT63_250'][:]
                elif plev > 0 and plev < 5:
                    rvT63_850 = nc.variables['rvT63_850'][:]
                    rvT63_500 = nc.variables['rvT63_500'][:]
                    rvT63_250 = nc.variables['rvT63_250'][:]
                elif plev == 0:
                    rvT63_850 = np.ones(len(lons))
                    rvT63_500 = np.ones(len(lons))
                    rvT63_250 = np.ones(len(lons))

                zg_avg_925 = np.ones(len(lons))
                zg_max_925 = np.ones(len(lons))
                zg_min_925 = np.ones(len(lons))
                zg_avg_500 = np.ones(len(lons))
                zg_max_500 = np.ones(len(lons))
                zg_min_500 = np.ones(len(lons))
                zg_avg_600 = np.ones(len(lons))
                zg_max_600 = np.ones(len(lons))
                zg_min_600 = np.ones(len(lons))
                zg_avg_250 = np.ones(len(lons))
                zg_max_250 = np.ones(len(lons))
                zg_min_250 = np.ones(len(lons))
                zg_avg_300 = np.ones(len(lons))
                zg_max_300 = np.ones(len(lons))
                zg_min_300 = np.ones(len(lons))
                if 'zg_avg_250' in variables:
                    zg_avg_250 = variables['zg_avg_250'][:]
                    zg_min_250 = variables['zg_min_250'][:]
                    zg_max_250 = variables['zg_max_250'][:]
                if 'zg_avg_300' in variables:
                    zg_avg_300 = variables['zg_avg_300'][:]
                    zg_max_300 = variables['zg_max_300'][:]
                    zg_min_300 = variables['zg_min_300'][:]
                if 'zg_avg_500' in variables:
                    zg_avg_500 = variables['zg_avg_500'][:]
                    zg_min_500 = variables['zg_min_500'][:]
                    zg_max_500 = variables['zg_max_500'][:]
                if 'zg_avg_600' in variables:
                    zg_avg_600 = variables['zg_avg_600'][:]
                    zg_min_600 = variables['zg_min_600'][:]
                    zg_max_600 = variables['zg_max_600'][:]
                if 'zg_avg_925' in variables:
                    zg_avg_925 = variables['zg_avg_925'][:]
                    zg_min_925 = variables['zg_min_925'][:]
                    zg_max_925 = variables['zg_max_925'][:]

                radius_8 = np.ones(len(lons))
                if 'radius_8' in variables:
                    radius_8 = variables['radius_8'][:]

                ace = np.ones(len(lons))
                if 'ace' in variables:
                    ace = nc.variables['ace'][:]
                ace_psl = np.ones(len(lons))
                if 'ace_psl' in variables:
                    ace_psl = nc.variables['ace_psl'][:]
                ike = np.ones(len(lons))
                if 'ike' in variables:
                    ike = nc.variables['ike'][:]
                pdi = np.ones(len(lons))
                if 'pdi' in nc.variables:
                    pdi = variables['pdi'][:]

                profile_len = 79
                rprof = np.ones(len(lons)*profile_len)
                if 'rprof' in variables:
                    rprof = nc.variables['rprof'][:]
                    print ('read rprof ',rprof)
                rprof_slp = np.ones(len(lons)*profile_len)
                if 'rprof_slp' in variables:
                    rprof_slp = nc.variables['rprof_slp'][:]

                print('ntracks ',ntracks)
                #print('mslp, w10m, vmax, vort ',psl[0:2], '; ', sfcWinds[0:2], '; ', vmaxs[0:2], '; ', vorts[0:2])

                for storm_no in range(ntracks):
                #for storm_no in range(2):
                    storm_obs = []
                    tcid = storm_no
                    #print 'storm_no ',storm_no
                    first_pt = first_pts[storm_no]
                    storm_length = storm_lengths[storm_no]
                    record_no = storm_length
                    index = indices[first_pt:first_pt+storm_length]
                    #print storm_no, index
                    for ip in index:
                        i = ip+first_pt
                        date = dtime[i]
                        #print storm_no, ip, i, first_pt, date

                        # lat, lon are the latitude and longitude associated with the tracked feature
                        lat = lats[i]
                        lon = lons[i]
                        vort = vorts[i]
                        vort_min = vorts_min[i]
                        psl = psls[i] * psl_scaling
                        sfcWind = sfcWinds[i]
                        vmax = vmaxs[i]
                        vmax_kts = vmax * scaling_ms_knots

                        rv850_t63_this = rvT63_850[i]
                        rv500_t63_this = rvT63_500[i]
                        rv250_t63_this = rvT63_250[i]
                        if plev >= 5:
                            rv700_t63_this = rvT63_700[i]
                            rv600_t63_this = rvT63_600[i]
                        else:
                            rv700_t63_this = -99.9
                            rv600_t63_this = -99.9

                        zg_avg_925_this = zg_avg_925[i]
                        zg_max_925_this = zg_max_925[i]
                        zg_min_925_this = zg_min_925[i]
                        zg_avg_500_this = zg_avg_500[i]
                        zg_max_500_this = zg_max_500[i]
                        zg_min_500_this = zg_min_500[i]
                        zg_avg_600_this = zg_avg_500[i]
                        zg_max_600_this = zg_max_500[i]
                        zg_min_600_this = zg_min_500[i]
                        zg_avg_250_this = zg_avg_250[i]
                        zg_max_250_this = zg_max_250[i]
                        zg_min_250_this = zg_min_250[i]
                        zg_avg_300_this = zg_avg_300[i]
                        zg_max_300_this = zg_max_300[i]
                        zg_min_300_this = zg_min_300[i]
                        radius_8_this = radius_8[i]
                        ace_this = ace[i]
                        ace_psl_this = ace_psl[i]
                        ike_this = ike[i]
                        pdi_this = pdi[i]
                        rprof_this = rprof[i:i+profile_len]
                        rprof_slp_this = rprof_slp[i:i+profile_len]
            
                        storm_obs.append(Observation(date, lat, lon, vort, vort_min, vmax, psl, rv850_t63_this, rv700_t63_this, rv600_t63_this, rv500_t63_this, rv250_t63_this, zg_avg_925_this, zg_max_925_this, zg_min_925_this, zg_avg_600_this, zg_max_600_this, zg_min_600_this, zg_avg_500_this, zg_max_500_this, zg_min_500_this, zg_avg_300_this, zg_max_300_this, zg_min_300_this, zg_avg_250_this, zg_max_250_this, zg_min_250_this, radius_8_this, ace_this, ace_psl_this, ike_this, pdi_this, rprof_this, rprof_slp_this, extras={'vmax_kts':vmax_kts, 'w10m':sfcWind}))

            # Yield storm
                    yield Storm(tcid, storm_obs, extras={})

if __name__ == '__main__':
    fname = os.path.join(SAMPLE_TRACK_DATA)
    print('Loading TRACK data from file:' , fname)    
    storms = list(load(fname, ex_cols=3, calendar='cftime'))
    print('Number of model storms: ', len(storms))
    
    # Print storm details:
    for storm in storms: 
        #print storm.snbr, storm.genesis_date()
        for ob in storm.obs:
            print(ob.date, ob.lon, ob.lat, ob.vmax, ob.extras['vmax_kts'], ob.mslp, ob.vort)
    print('Number of model storms: ', len(storms))
    
