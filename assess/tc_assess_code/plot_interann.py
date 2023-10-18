'''
Plot the TC frequency and ACE for all basins for a given year, including by storm category
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

sys.path.append('/home/h06/hadom/workspace/tenten/storms/assessment')
import tc_assess_code

