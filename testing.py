import pandas as pd
pd.options.mode.chained_assignment = None
#pd.set_option("display.max_rows", None, "display.max_columns", None)
import xarray as xr
import numpy as np
import os
from datetime import datetime, timedelta
from collections import defaultdict
import pyproj
proj = pyproj.Proj(init='epsg:3857')

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

from utils.functions import get_iem_bufkit, read_bufkit
from utils.bufkit import to_bufkit_test
from stdatmos import std_atmosphere_pres
from interpolate.tools import calc_distance
#from utils.configs import LOAD_VARS

stdatm = std_atmosphere_pres()
DIR = os.path.dirname(__file__) or "."
def load_metadata(meta_fname=("%s/airports.dat" % DIR)):
    """
    load_meta
    Load in our database of airport codes.
    """
    meta_cols = ['code', 'id', 'synop', 'lat', 'lon', 'elev', 'name']
    meta_types = {'code': int, 'id': str, 'synop': int, 'lat': float, 'lon': float,
                  'elev': int, 'name': str}

    meta_airport = {}
    with open(meta_fname) as fmeta:
        for line in fmeta:
            line_dict = {col: val for col, val in zip(meta_cols, line.strip().split(None, 6)) }

            for col in line_dict.keys():
                line_dict[col] = meta_types[col](line_dict[col])

            id = line_dict.pop('id')
            meta_airport[id] = line_dict
    return meta_airport

def interp_pres(p, pres, field):
    """Generic interpolation routine. Converts pressures into log10 coordinates
    Parameters
    ----------
    p : number, numpy array
        Pressure (hPa) of the level for which the field variable is desired
    pres : numpy array
        The array of pressure
    field : numpy array
        The variable which is being interpolated
    Returns
    -------
    Value of the 'field' variable at the given pressure : number, numpy array
    """
    field_intrp = np.interp(np.log10(p)[::-1], np.log10(pres)[::-1], field[::-1])
    return field_intrp[::-1]


LOAD_VARS = [
    'temperature', 'dewpoint', 'windSpeed', 'windDir', 'altitude', 'invTime', 'latitude',
    'longitude', 'dataType', 'flight', 'orig_airport_id', 'dest_airport_id',
    'sounding_airport_id'
]
MS2KTS = 1.94384
MILETOM = 1609.34
start_time, end_time = datetime(2021, 6, 20, 23), datetime(2021, 6, 21, 6)
DATA_PATH = '/Users/leecarlaw/scripts/acars/netcdf/acars'
DELTA = 30
RADIUS = 100

# Read in hourly files and output merged pandas DataFrame
files = []
while start_time <= end_time:
    timestring = start_time.strftime('%Y%m%d_%H00')
    files.append("%s/%s" % (DATA_PATH, timestring))
    start_time += timedelta(hours=1)
temp = []
for f in files:
    temp_data = {}
    ds = xr.open_dataset(f)
    for var in LOAD_VARS: temp_data[var] = ds[var].values

    df = pd.DataFrame.from_dict(temp_data)
    df['flight'] = df['flight'].astype(int)
    df['ob_time'] = pd.to_datetime(df['invTime'], unit='s')
    df.drop(columns=['invTime'])
    temp.append(df)
df = pd.concat(temp)

# Set the analysis time bounds. Download RAP data and time-interpolate to match at DELTAs
dt_start, dt_end = datetime(2021, 6, 20, 23), datetime(2021, 6, 21, 6, 0)
dt = dt_start
while dt <= dt_end + timedelta(hours=1):
    get_iem_bufkit('MDW', dt, path='./bufkit')
    dt += timedelta(hours=1)

dt = dt_start
obs = []
meta = load_metadata()['MDW']
meta['id'] = 'MDW'
site_loc = proj(meta['lon'], meta['lat'])
site_x, site_y = np.array(site_loc[0]), np.array(site_loc[1])
while dt <= dt_end:
    df['x'], df['y'] = proj(np.array(df['longitude']), np.array(df['latitude']))
    points = list(zip(df['x'], df['y']))
    df['distance'] = calc_distance(points, site_x, site_y)
    dt_1, dt_2 = dt - timedelta(minutes=DELTA/2), dt + timedelta(minutes=DELTA/2)
    subset = df.loc[(df['ob_time'] >= dt_1) & (df['ob_time'] < dt_2) & \
                    (df['distance'] < RADIUS * MILETOM)]
    pressure = stdatm(np.array(subset['altitude'])) / 100.
    out = subset[['temperature', 'dewpoint', 'windSpeed', 'windDir', 'altitude',
                  'longitude', 'latitude', 'distance', 'ob_time']]
    out['pressure'] = pressure
    out['datetime'] = dt
    obs.append(out)

    dt += timedelta(minutes=DELTA)

obs = pd.concat(obs)
obs.rename(columns={'altitude': 'hght'}, inplace=True)
obs['temperature'] = obs['temperature'] - 273.15
obs['dewpoint'] = obs['dewpoint'] - 273.15
obs['windSpeed'] = obs['windSpeed'] * MS2KTS
obs.to_csv('./master.csv')

dt = dt_start
while dt <= dt_end:
    df = obs.loc[obs['datetime'] == dt]
    dt1 = datetime(dt.year, dt.month, dt.day, dt.hour)
    dt2 = dt1 + timedelta(hours=1)
    weight = 1 - (dt - dt1).total_seconds() / 3600. # How much to weight dt1
    rap1 = read_bufkit('./bufkit/rap_MDW_%s.buf' % (dt1.strftime('%Y%m%d%H')), dt1)
    rap2 = read_bufkit('./bufkit/rap_MDW_%s.buf' % (dt2.strftime('%Y%m%d%H')), dt2)
    RAP = {}
    if len(df) > 10:
        df.drop_duplicates(subset='hght', inplace=True)
        df.sort_values('pressure', ascending=False, inplace=True)
        df.reset_index(inplace=True)

        # Interpolate RAP data: 1st to observed pressure levels in logp space and then
        # linearlly in time
        RAP['pressure'] = np.array(df['pressure'])
        RAP['dewpoint1'] = interp_pres(df['pressure'], rap1['pressure'], rap1['dewpoint'])
        RAP['dewpoint2'] = interp_pres(df['pressure'], rap2['pressure'], rap2['dewpoint'])
        RAP['dewpoint'] = RAP['dewpoint1'] * weight + RAP['dewpoint2'] * (1-weight)
        RAP.pop('dewpoint1')
        RAP.pop('dewpoint2')
        RAP = pd.DataFrame.from_dict(RAP)
        df['dewpoint'] = df['dewpoint'].combine_first(RAP['dewpoint'])
        df.drop(columns=['index'], inplace=True)

        i = 0
        while i < len(df):
            diffs = df.iloc[i]['pressure'] - df['pressure']
            df['diffs'] = diffs
            temp = df.loc[(df['diffs'] >= 0) & (df['diffs'] < 10)]

            df.iloc[i] = temp.mean()
            df.drop(temp.index[1:], inplace=True)
            df.loc[i, 'datetime'] = dt
            df.drop(columns=['diffs'], inplace=True)
            i += 1

        df.dropna(how='all', inplace=True)
        df.reset_index(inplace=True)
        df.fillna(-9999.0, inplace=True)
        print(len(df))
        lines = to_bufkit_test(df, meta)
        with open('./output/ACARS_%s.buf' % (dt.strftime("%Y%m%d_%H%M")), 'w') as fsnd:
            fsnd.write("".join(lines))
    dt += timedelta(minutes=DELTA)


##########################################################################################
# Plotting Routines
##########################################################################################
SHAPEFILES = '/Users/leecarlaw/shapefiles'
reader = shpreader.Reader('%s/countyl010g.shp' % (SHAPEFILES))
counties = list(reader.geometries())
COUNTIES = cfeature.ShapelyFeature(counties, ccrs.PlateCarree())
dpi = 100
fig = plt.figure(figsize=(11,8), dpi=dpi)
proj = ccrs.LambertConformal(central_latitude = 42, central_longitude = -89,
                             standard_parallels = (25, 25))
ax = fig.add_axes([0., 0., 1, 0.98], projection=proj)
ax.set_extent([-90.25, -85.75, 40.75, 42.95], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.STATES.with_scale('50m'), linewidth=3, zorder=3)
ax.add_feature(COUNTIES, facecolor='none', edgecolor='darkgrey', linewidth=0.5, zorder=2)
ax.add_feature(cfeature.LAKES.with_scale('50m'), zorder=1)
cax = fig.add_axes([0.03, 0.05, 0.025, 0.875])

dt = dt_start
while dt <= dt_end:
    data = obs.loc[obs['datetime'] == dt]
    dt1 = dt - timedelta(minutes=DELTA/2)
    dt2 = dt + timedelta(minutes=DELTA/2)
    info = dt1.strftime('%b %d %Y %H:%M:%S') + ' - ' + dt2.strftime('%b %d %Y %H:%M:%S')
    t0 = ax.text(0.82, 0.9, "# Obs: %s" % (len(data)), fontsize=20,
                 bbox={'facecolor':'wheat', 'alpha': 0.75, 'edgecolor':'k', 'boxstyle':'round'},
                 transform=ax.transAxes, zorder=1000)
    t1 = ax.annotate(info, xy=(0.995, 1.0), va='bottom', ha='right',
                     xycoords='axes fraction', fontsize=15)
    cf = ax.scatter(data['longitude'], data['latitude'], c=data['hght'],
                    s=50, edgecolor='k', vmin=0, vmax=9000, zorder=999,
                    cmap=plt.cm.magma_r,
                    transform=ccrs.PlateCarree())

    if dt == dt_start:
        cb = plt.colorbar(cf, orientation='vertical', cax=cax,
                          label='Aircraft Altitude (m)')
    plt.savefig('./plots/%s.png' % (dt.strftime('%Y%m%d%H%M')), bbox_inches='tight')
    cf.remove()
    t0.remove()
    t1.remove()
    dt += timedelta(minutes=DELTA)
