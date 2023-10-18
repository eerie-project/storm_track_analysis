import numpy
import os.path
import collections

#: Store observations as a named tuple (this enables access to data by its 
#: name instead of position index)
Observation = collections.namedtuple('Observation', ['date', 'lat', 'lon', 'vort', 'vmax', 
                                                     'mslp', 't63_850', 't63_2', 't63_3', 't63_4', 't63_5', 't63_diff', 'ace_index', 'extras'])


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
    
    @property
    def t63_850_max(self):
        """ The maximum T63 vorticity attained by the storm during its lifetime """
        if self.obs_at_genesis().lat > 0:
            return max(ob.t63_850 for ob in self.obs)
        else:
            return min(ob.t63_850 for ob in self.obs)
    
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
    
