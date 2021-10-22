import numpy as np
from datetime import datetime, timedelta
from collections import defaultdict
import os, subprocess, sys
import argparse
import pandas as pd
import xarray as xr
import pyproj

import interpolate.tools as tools
from utils.functions import *
from utils.constants import ZEROCNK, MS2KTS, MILETOM
from utils.configs import *
from utils.bufkit import *

import metpy.calc as mpcalc
from metpy.units import units
from metpy.interpolate.points import barnes_point, cressman_point

from stdatmos import std_atmosphere_pres

DIR = os.path.dirname(__file__) or "."
stdatm = std_atmosphere_pres()
proj = pyproj.Proj(init='epsg:3857')
def load_metadata(meta_fname=("%s/airports.dat" % DIR)):
    """
    load_meta
    Load in our database of airport codes.
    """
    meta_cols = ['code', 'id', 'synop', 'lat', 'lon', 'elev', 'name']
    meta_types = {'code': int, 'id': str, 'synop': int, 'lat': float, 'lon': float, 'elev': int, 'name': str}

    meta_airport = {}
    with open(meta_fname) as fmeta:
        for line in fmeta:
            line_dict = {col: val for col, val in zip(meta_cols, line.strip().split(None, 6)) }

            for col in line_dict.keys():
                line_dict[col] = meta_types[col](line_dict[col])

            id = line_dict.pop('id')
            #code = line_dict.pop('code')
            meta_airport[id] = line_dict
    return meta_airport

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-a', '--analysis', dest='analysis',
                    help="Analysis hour: YYYY-mm-dd/HH")
    ap.add_argument('-n', '--n_hours', dest='n_hours',
                    help="Number of hours behind analysis time to look")
    ap.add_argument('-p', '--data-path', dest='data_path',
                    help="Local path for data read")
    args = ap.parse_args()

    if not args.data_path: args.data_path = NETCDF + '/acars/'
    if not args.n_hours:
        args.n_hours = 1
    else:
        args.n_hours = int(args.n_hours)

    metadata = load_metadata()
    try:
        dt_analysis = datetime.strptime(args.analysis, '%Y-%m-%d/%H')
    except:
        print("Poorly-formatted start and/or end times [YYYY-MM-DD/HH]")
        sys.exit(1)

    files = []
    target_dts = []
    dt = dt_analysis - timedelta(hours=args.n_hours)
    while dt <= dt_analysis:
        filename = "%s/%s" % (args.data_path, dt.strftime("%Y%m%d_%H00"))
        files.append(filename)
        target_dts.append(dt)
        dt += timedelta(hours=1)

    # Loading in the ACARS data
    prof_data = pd.DataFrame()
    for f in files:
        temp_data = {}
        ds = xr.open_dataset(f)
        for var in LOAD_VARS: temp_data[var] = ds[var].values
        temp = pd.DataFrame.from_dict(temp_data)
        prof_data = prof_data.append(temp, ignore_index=True)
    prof_data.dropna(subset=['soundingSecs'], inplace=True)
    prof_data['datetime']  = pd.to_datetime(prof_data['soundingSecs'])

    # Data processing loops
    profiles = defaultdict(list)
    for site in site_ids:
        try:
            meta = metadata[site]
        except:
            pass

        meta['id'] = site
        site_loc = proj(meta['lon'], meta['lat'])
        site_x, site_y = np.array(site_loc[0]), np.array(site_loc[1])
        prof_data['x'], prof_data['y'] = proj(np.array(prof_data['longitude']),
                                              np.array(prof_data['latitude']))
        points = list(zip(prof_data['x'], prof_data['y']))
        prof_data['distance'] = tools.calc_distance(points, site_x, site_y)
        profiles = prof_data[prof_data['distance'] < RADIUS * MILETOM]

        # Compute pressure using the standard atmosphere model. Convert wind speed &
        # direction to u- and v- components. Add in metadata. Units conversions
        wind_speed = (np.array(profiles['windSpeed']) * MS2KTS) * units('kts')
        wind_dir = np.array(profiles['windDir']) * units.deg
        profiles['u'], profiles['v'] = mpcalc.wind_components(wind_speed, wind_dir)
        profiles['pressure'] = stdatm(np.array(profiles['altitude'])) / 100.
        profiles['temperature'] = np.array(profiles['temperature']) - 273.15
        profiles['dewpoint'] = np.array(profiles['dewpoint']) - 273.15
        profiles['hght'] = profiles['altitude']
        master_df = pd.DataFrame.from_dict(profiles)


        observations = defaultdict(list)
        background = defaultdict(list)
        for dt in target_dts:
            print('==============  ', dt, '  ==============')
            cutoff_start = dt - timedelta(minutes=25)
            cutoff_end = dt + timedelta(minutes=25)
            mask  = (master_df['datetime'] > cutoff_start) & \
                    (master_df['datetime'] <= cutoff_end)
            df = master_df[mask]

            rap_data = get_iem_bufkit(site, dt, model='rap', path=DIR + '/bufkit')
            background['datetime'].append(dt)
            observations['datetime'].append(dt)
            for parm in rap_data.keys(): background[parm].append(rap_data[parm])

            obs = defaultdict(list)
            parms = ['dewpoint', 'pressure', 'u', 'v', 'temperature']

            # Scale the vertical threshold limits. We want smaller segments near the
            # surface and larger at the top.
            norm = (rap_data['hght'] - np.min(rap_data['hght']))/np.ptp(rap_data['hght'])
            scale = 200*np.exp(norm*3)

            temp = []
            for k in range(len(rap_data['hght'])):
                height = rap_data['hght'][k]
                obs['hght'].append(height)
                mask = (df['hght'] >= height - scale[k]) & (df['hght'] <= height + scale[k])
                O = df[mask]
                for parm in parms:
                    data = O.dropna(subset=[parm])
                    n_obs = len(data)
                    obs[parm].append(rap_data[parm][k])
                    if n_obs > 1:
                        # Consistency checks
                        med = np.percentile(data[parm], 50)
                        sigma = np.std(data[parm])
                        min_ = med - (1.*sigma)
                        max_ = med + (1.*sigma)
                        qc_mask = (data[parm] > min_) & (data[parm] < max_)
                        n_rejected = len(qc_mask) - qc_mask.sum()
                        passed = data[qc_mask]

                        if len(passed) > 2:
                            coords = list(zip(passed['x'], passed['y']))
                            kappa = tools.calc_kappa(tools.average_spacing(coords))

                            # Create an 'effective distance' from slant ranges
                            height_delta = abs(passed['hght'] - height)
                            eff_distance = np.sqrt(passed['distance']**2 + ((100*height_delta)**2))
                            #weights = np.exp(-passed['distance']**2 / kappa)
                            weights = np.exp(-eff_distance**2 / kappa)
                            total_weights = np.sum(weights)
                            if np.sum(weights) > 0:
                                barnes_val = np.sum(passed[parm] * (weights / total_weights))
                                obs[parm][-1] = barnes_val

            for parm in parms:
                observations[parm].append(obs[parm])
            wdir, wspd = comp2vec(obs['u'], obs['v'])
            observations['windSpeed'].append(wspd)
            observations['windDir'].append(wdir)
            observations['hght'].append(obs['hght'])

        background = pd.DataFrame.from_dict(background)
        observations = pd.DataFrame.from_dict(observations)
        lines = to_bufkit('./', observations, meta, label="ACARS")
        lines.extend(to_bufkit('./', background, meta, label="RAP"))

        out_file = "%s/ACARS_%s.buf" % ('/Users/leecarlaw/bufkit/', site)
        with open(out_file, 'w') as fsnd: fsnd.write("".join(lines))

        """


# Don't run on module import
if __name__ == "__main__":
    main()
