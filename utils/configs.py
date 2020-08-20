#!/usr/bin/env python3
"""Configuration specifications"""

import os
from datetime import datetime
import pandas as pd
curr_path = os.path.dirname(__file__) or "."

# --------------------------------------------------------------------------------------
# Variables in this first block **must** be edited to run on your local system
# --------------------------------------------------------------------------------------
# Various PATH configurations.
# DATA : Storage location for python-pickled dictionary files. These contain processed
#        AMDAR data from the MADIS netCDF files.
# NETCDF : Storage location for MADIS netCDF files.
# SOUNDINGS : Output location for sounding .buf files.
# WGET : Path to the wget binary on this system.

DATA = "/Users/leecarlaw/scripts/acars/data"
SOUNDINGS = "/Users/leecarlaw/scripts/acars/soundings"
NETCDF = "/Users/leecarlaw/scripts/acars/netcdf"
WGET = "/usr/local/bin/wget"

# --------------------------------------------------------------------------------------
# Additional editable configurations
# --------------------------------------------------------------------------------------
site_ids = [
    'MDW', 'ORD', 'DAL',
    'MKE',
]

# Number of hours to keep pickled python dictionary files in the top-level /data_store
# directory before sweeping into the archive.
purge_hours = 24

# Various QC Thresholds
T_QC = 7.5
TD_QC = 25.

# --------------------------------------------------------------------------------------
# Variables in this block should not be changed, unless you have a specific reason to
# do so!
# --------------------------------------------------------------------------------------
base_url = "https://madis-data.ncep.noaa.gov/madisNoaa/data/point/acars/netcdf"
archive_url = "https://madis-data.ncep.noaa.gov/madisNoaa/data/archive"
bufkit_url = "https://mtarchive.geol.iastate.edu"
epoch = datetime(1970, 1, 1, 0)
min_pressure_levels = 21
max_pressure_levels = 55

df = pd.read_csv(curr_path + '/L62.txt', delim_whitespace=True)
a_k = df['a'] / 100.
b_k = df['b']
