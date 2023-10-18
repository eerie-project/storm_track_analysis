""" 
Load function to read in Reading Universities TRACK output. 
Should work for individual ensemble member files and also the 
'combine' file. 

Note: All files need to contain the DATES of the storms, not the 
original timestep as output by TRACK.

Note: Any fields in addition to the full field vorticity, mslp
and 925 hPa max winds (e.g. 10m winds or precipitation) are not 
currently stored. If you want these values you need to read in 
the data file and include the variables in the 'extras' dictionary. 


"""
import os.path, sys
import datetime
import numpy as np
try:
    import cftime
except ImportError:
    import netcdftime as cftime
import six, numpy

# Set path for sample model data
SAMPLE_DATA_PATH = os.path.join(os.path.dirname(__file__), 'sample_data')

SAMPLE_TRACK_DATA = os.path.join(SAMPLE_DATA_PATH, 
                                 'combined_ff_trs.vor_10m_fullgrid_N512_xgxqe_L5.new_20002011.date')

WIND_MAX = 400.
WIND_MIN = 0.

import collections

#: Store observations as a named tuple (this enables access to data by its 
#: name instead of position index)
Observation = collections.namedtuple('Observation', ['date', 'lat', 'lon', 'vort', 'vmax', 
                                                     'mslp', 'ace_index', 'extras'])


class Observation(Observation):
    """  Represents a single observation of a model tropical storm. """
    
    def six_hourly_timestep(self):
        """ Returns True if a storm record is taken at 00, 06, 12 or 18Z only """
        return self.date.hour in (0,6,12,18) and self.date.minute == 0 and self.date.second == 0
    
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
            return min(ob.vort for ob in self.obs)
    
    def __len__(self):
        """ The total number of observations for the storm """
        return len(self.obs)
    
    def nrecords(self):
        """ The total number of records/observations for the storm """
        return len(self)
    
    def n_six_hourly_records(self):
        """ The total number of records/observations for the storm ay 00, 06, 12 and 18Z. """
        return len([ob for ob in self.obs if ob.six_hourly_timestep()])    
    
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
                ace_index += numpy.square(ob.extras['vmax_kts'])/10000.
        return round(ace_index, 2)

    def ace_index_no6hrcheck(self):
        """ The accumulated cyclone energy index for the storm. Calculated as the square of the
        storms maximum wind speed every 6 hours (0, 6, 12, 18Z) throughout its lifetime. Observations
        of the storm taken in between these records are not currently used. Returns value rounded to
        2 decimal places. Wind speed units: knots """
        ace_index = 0
        for ob in self.obs:
            ace_index += numpy.square(ob.extras['vmax_kts'])/10000.
            #print('vmax_kts ',ob.extras['vmax_kts'])
        return round(ace_index, 2)
    
    def obs_at_vmax(self):
        """Return the maximum observed vmax Observation instance. If there is more than one obs 
        at vmax then it returns the first instance """
        return max(self.obs, key=lambda ob: ob.vmax)
    
    def obs_at_genesis(self):
        """Returns the Observation instance for the first date that a storm becomes active """       
        for ob in self.obs:
            return ob
        else:
            raise ValueError('model storm was never born :-(')

    def obs_at_lysis(self):
        """Returns the Observation instance for the last date that a storm was active """    
        return [ob for ob in self.obs][-1]

def load(fh, ex_cols=0, calendar=None, mslp_conversion = 1.0, variable_indices={}):
    """
    Reads model tropical storm tracking output from Reading Universities TRACK 
    algorithm. Note: lat, lon, vorticity, maximum wind speed and minimum central
    pressure values are taken from the full resolution field, not T42. The lat/
    lon values correspond to the location of maximum 850 hPa relative vorticity 
    (unless data are unavailable, in which case the original T42 resolution lat/
    lon values are used). 
    
    If you have additional fields/columns after the full field vorticity, max 
    wind and mslp data, then you need to count the number of columms these 
    take up and set this as the ex_cols value (this value does not include 
    &'s or comma's, just data). For example, for a file containing additional 
    10m wind information ex_cols=3 (lat, lon, 10m wind speed value).
    
    Note: if you are using model data which uses a 12 months x 30 daycalendar 
    then you need to set calendar to 'netcdftime'. Default it to use gregorian 
    calendar.

    Inputs:
       variable_indices: dict. Dict entries for w10m, mslp, w925, fullvort 
                         contain the (reversed-order) indix for the column
                         that this variable appears in
    
    Assumption for TRACK TC ascii files (if ex_cols=3):
       last element: 10m windspeed
       second to last: MSLP
       third to last: full windspeed (925 if available, 850 if not)
       fourth to last: full vorticity
       hence, when reversed, the indices for these are: 1,4,7,10

    """
    # indices for variables (when columns are reversed)
    variables = ['w10m', 'mslp', 'w925', 'fullvort']
    var_index = {}
    
    if variable_indices == {}:
        if ex_cols == 3:
            var_index['w10m'] = 1
            var_index['mslp'] = 4
            var_index['w925'] = 7
            var_index['fullvort'] = 10
    else:
        var_index = variable_indices

    for variable in variables:
        if variable not in var_index.keys():
            raise Exception('Variable not in variable_indices '+variable)

    # allow users to pass a filename instead of a file handle.
    if isinstance(fh, six.string_types):
        with open(fh, 'r') as fh:
            for data in load(fh, ex_cols=ex_cols, calendar=calendar, variable_indices=var_index):
                yield data
                
    else:
        # for each line in the file handle            
        for line in fh:
            if line.startswith('TRACK_NUM'):
                split_line = line.split()
                if split_line[2] == 'ADD_FLD':
                    number_fields = int(split_line[3])
                else:
                    raise ValueError('Unexpected line in TRACK output file.')
                
            if line.startswith('TRACK_ID'):
                # This is a new storm. Store the storm number.
                _, snbr, _, _ = line.split()
                snbr =  int(snbr.strip())
                
                # Now get the number of observation records stored in the next line                
                next_line = next(fh)
                if next_line.startswith('POINT_NUM'):
                    _, n_records = next_line.split()
                    n_records = int(n_records)
                else:
                    raise ValueError('Unexpected line in TRACK output file.')
                            
                # Create a new observations list
                storm_obs = []
                
                """ Read in the storm's observations """     
                # For each observation record            
                for _, obs_line in zip(list(range(n_records)), fh):
                    
                    # Get each observation element
                    split_line =  obs_line.strip().split('&')

                    # Get observation date and T42 lat lon location in case higher 
                    # resolution data are not available
                    date, tmp_lon, tmp_lat, _ = split_line[0].split()
                    
                    yr = int(date[0:4])
                    mn = int(date[4:6])
                    dy = int(date[6:8])
                    hr = int(date[8:10])

                    if calendar == 'netcdftime':
                        date = cftime.datetime(yr, mn, dy, hour=hr)
                    elif calendar == '360_day':
                        start_date = cftime.datetime(yr, mn, dy, hour=hr, calendar = calendar)
                        time_unit = 'hours since 1950-01-01 00:00:00'
                        t = cftime.date2num(start_date, time_unit, calendar = calendar)
                        date = cftime.num2date(t, time_unit, calendar = calendar)
                    else:
                        date = datetime.datetime.strptime(date.strip(), '%Y%m%d%H')
                        
                    # Get full resolution mslp - second from last
                    mslp = split_line[::-1][variable_indices['mslp']]
                    mslp = float(mslp)/mslp_conversion
                    if mslp > 1.0e4:
                        mslp /= 100
                    mslp = int(round(mslp))
               
                    # Get full resolution 925hPa maximum wind speed (m/s) - third from last (3 columns for each value)
                    vmax = float(split_line[::-1][variable_indices['w925']])
                    if np.abs(vmax) > 100:
                        print('vmax big ',vmax, yr, mn, dy, hr)
                    if vmax > WIND_MAX or vmax < WIND_MIN:
                        vmax = 0.
                        #continue
                    
                    # Also store vmax in knots (1 m/s = 1.944 kts) to match observations.
                    vmax_kts = vmax * 1.944
                    ace_index = numpy.square(vmax_kts)/10000.
                    
                    # Get full resolution 850 hPa maximum vorticity (s-1)
                    vort = float(split_line[::-1][variable_indices['fullvort']])
                    
                    # Get storm location of maximum vorticity (full resolution field)
                    lat = float(split_line[::-1][variable_indices['fullvort']+1])
                    lon = float(split_line[::-1][variable_indices['fullvort']+2])
                    
                    # wind speed (last column)
                    if ex_cols > 0:
                        w10m = float(split_line[::-1][variable_indices['w10m']])
                    else:
                        w10m = vmax*0.9
                    if w10m > WIND_MAX or w10m < WIND_MIN:
                        w10m = 0.
                    #print('vmax, vmax_kts ',date, w10m, vmax, vmax_kts)

                    # If higher resolution lat/lon data is not available then use lat 
                    # lon from T42 resolution data
                    if lat == 1e12 or lon == 1e12 or lat == 1.0e25 or lon == 1.0e25:
                        lat = float(tmp_lat)
                        lon = float(tmp_lon)
                                        
                    lat = float(tmp_lat)
                    lon = float(tmp_lon)
                                        
                    #print('split_line(reversed) ',split_line[::-1])
                    #print('mslp, w10m, vmax, vort ',date, mslp, '; ', w10m, '; ', vmax, '; ', vort)
                    # Store observations
                    storm_obs.append(Observation(date, lat, lon, vort, vmax, mslp, ace_index,
                                                          extras={'vmax_kts':vmax_kts, 'w10m':w10m}))
                    
                # Yield storm
                yield Storm(snbr, storm_obs, extras={})
                
if __name__ == '__main__':
    fname = os.path.join(SAMPLE_TRACK_DATA)
    print('Loading TRACK data from file:' , fname)    
    storms = list(load(fname, ex_cols=3, calendar='netcdftime'))
    print('Number of model storms: ', len(storms))
    
    # Print storm details:
    for storm in storms: 
        #print storm.snbr, storm.genesis_date()
        for ob in storm.obs:
            print(ob.date, ob.lon, ob.lat, ob.vmax, ob.extras['vmax_kts'], ob.mslp, ob.vort)
    print('Number of model storms: ', len(storms))
    
