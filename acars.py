#!/usr/bin/env python3
import numpy as np
from datetime import datetime, timedelta
import pickle
import urllib.request as urlreq
from urllib.error import HTTPError
import requests
from collections import defaultdict
import os, subprocess
import glob
import argparse
import math
import time
import pandas as pd

import utils.thermo as thermo
from utils.constants import ZEROCNK
from utils.bufkit_helper import EXTRA_LINES, FILE_HEADER, SFC_HEADER, SFC_VALS
from utils.cron_helper import *
from utils.configs import *

def to_bufkit(output_path, data, rap_data=None, label=None):
    """Writes output to BUFKIT-readble .buf file.

    BUFKIT is extremely unforgiving when it comes to formatting the input files. Many
    of these odds and ends aren't readily apparent, although you'll know you've hit
    them when the entire program crashes with various cryptic VB-based error messages!

    Various parameter constraints I've been able to deduce include (BUFKIT-19):
        :: At least 5 individual time steps (0,1,2,3,4)
        :: Minimum number of pressure levels @ individual timestep = 20 (inclusive)
        :: Maximum number of pressure levels @ individual timestep = 67 (inclusive)
        :: Minimum TOTAL number of pressure levels in a file = 201
        :: ?? Top pressure level must be < 480 mb (divide by zero error otherwise)
        :: ?? Maximum number of 'forecast' times is (about?) 120
        :: Time steps must be consistent throughout the entire file
        :: ?? PMSL and PRES values must exist in the surface variables section

    Parameters
    ----------
    output_path : string
        Path on local system to write BUFKIT file to
    data : pd.DataFrame
        ACARS data stored in a Pandas DataFrame
    rap_data : python dictionary (optional)
        If supplied, augments the ACARS profile with RAP data above the last available
        pressure level. Otherwise, uses the standard data in extra_lines to extend.

    Outputs
    -------
    Writes a .buf file to output_path.
    """
    min_levs = 20
    max_levs = 67
    min_times = 5

    ID = data.iloc[0]['apid']
    n_times = len(data)
    time_store = []
    pres_store = []
    date_last = pd.to_datetime(data.iloc[-1]['dt'])
    knt = 0
    total_lev_knt = 0

    snd_lines = ["%s\r\n" % (label)]
    snd_lines.extend(FILE_HEADER)

    for t in range(n_times):
        pres_store.append(data['pressure'].iloc[t][0])
    sfc_p_std = np.std(pres_store)
    sfc_p_mean = np.average(pres_store)
    p_bound_min = sfc_p_mean - (4 * sfc_p_std)
    p_bound_max = sfc_p_mean + (4 * sfc_p_std)
    # ----------------------------------------------------------------------------------
    # Data for outputting
    # ----------------------------------------------------------------------------------
    for t in range(n_times):
        #date_ob = data.iloc[t]['dt']
        date_ob = date_last - timedelta(minutes=60*(n_times-t-1))
        date = datetime.strftime(date_ob, "%Y%m%d%H%M")
        time_store.append(datetime.strftime(date_ob, "%Y%m%d/%H%M"))
        levs = np.array(data.iloc[t]['pressure'])
        t_out = np.array(data.iloc[t]['temperature'])
        td_out = np.array(data.iloc[t]['dewpoint'])
        wdir_out = np.array(data.iloc[t]['windDir'])
        wspd_out = np.array(data.iloc[t]['windSpeed'])

        # A few last-minute sanity checks
        wspd_out = np.where(wspd_out > 500, 0., wspd_out)
        wdir_out = np.where(wdir_out > 360, 0., wdir_out)
        hght_out = np.array(data.iloc[t]['altitude'])
        out_time = "%s%s%s/%s%s"%(date[2:4],date[4:6],date[6:8],date[8:10],date[10:12])

        snd_lines.extend([
            "\r\n",
            "STID = K%s STNM = %s0 TIME = %s\r\n" % (ID, 9999, out_time),
            "SLAT = %s SLON = %s SELV = %s\r\n" % (data.iloc[0]['lat'],
                                                   data.iloc[0]['lon'],
                                                   data.iloc[0]['elev']),
            "STIM = %s\r\n" % int(knt),
            "\r\n",
            "SHOW = -9999.00 LIFT = -9999.00 SWET = -9999.00 KINX = -9999.00\r\n",
            "LCLP = -9999.00 PWAT = -9999.00 TOTL = -9999.00 CAPE = -9999.00\r\n",
            "LCLT = -9999.00 CINS = -9999.00 EQLV = -9999.00 LFCT = -9999.00\r\n",
            "BRCH = -9999.00\r\n",
            "\r\n",
            "PRES TMPC TMWC DWPC THTE DRCT SKNT OMEG\r\n",
            "CFRL HGHT\r\n"]
        )
        lev_knt = 0
        prof = defaultdict(list)
        for row in range(len(t_out)):
            if levs[0] < p_bound_min:
                timestamp("WARN")
                print("Tossing %s and Flight # %s due to bad sfc pres." % (out_time,
                                                                data.iloc[t]['flight']))
                break
            if lev_knt >= max_levs: break
            thetae = thermo.theta_e(levs[row], t_out[row], td_out[row])

            prof['pres'].append(levs[row])
            prof['thte'].append(thetae)
            prof['tmpc'].append(t_out[row])
            prof['tmwc'].append(t_out[row]+td_out[row] / 2.)
            prof['dwpc'].append(td_out[row])
            prof['wdir'].append(wdir_out[row])
            prof['wspd'].append(wspd_out[row])
            prof['hght'].append(hght_out[row])
        p_sfc = levs[0]
        #pres_store.append(p_sfc)

        # Use the RAP data to extend our profile up to a maximum of 67 pressure levels.
        # If not available, extend using some random data (see utils.bufkit_helper).
        if rap_data is not None:
            num_extra_levs = max_levs - lev_knt
            RAP = rap_data[data.iloc[t]['dt']]
            top_pressure = levs[row]
            rap_idx = np.where(RAP['pressure'] < top_pressure)[0]
            counter = rap_idx[0]
            while (lev_knt < max_levs) and (counter < rap_idx[-1]):
                prof['pres'].append(RAP['pressure'][counter])
                prof['thte'].append(RAP['thte'][counter])
                prof['tmpc'].append(RAP['temperature'][counter])
                prof['tmwc'].append(RAP['tmwc'][counter])
                prof['dwpc'].append(RAP['dewpoint'][counter])
                prof['wdir'].append(RAP['windDir'][counter])
                prof['wspd'].append(RAP['windSpeed'][counter])
                prof['hght'].append(RAP['altitude'][counter])
                lev_knt += 1
                counter += 1
        knt += 1

        # In order to ensure consistent pressure levels which BUFKIT seems to desire,
        # derive these using a hybrid sigma vertical pressure coordinate system based
        # on the ECMWF L62 grid.
        p_levs = []
        for k in range(len(b_k)):
            p = a_k.iloc[k] + b_k.iloc[k] * p_sfc
            p_levs.append(p)
        p_levs = p_levs[::-1]

        PROF = {}
        for var in prof.keys():
            PROF[var] = interp_pres(p_levs, prof['pres'], prof[var])
            PROF[var] = np.where(PROF[var] < -1000., -9999., PROF[var])

        # Appending to the output array
        for row in range(PROF['pres'].shape[0]):
            out_line = "%s %s %s %s %s %s %s %s\r\n" % (round(PROF['pres'][row],2),
                                                        round(PROF['tmpc'][row],2),
                                                        round(PROF['tmwc'][row],2),
                                                        round(PROF['dwpc'][row],2),
                                                        round(PROF['thte'][row],2),
                                                        round(PROF['wdir'][row],2),
                                                        round(PROF['wspd'][row],2),
                                                        -9999.)
            snd_lines.append(out_line)
            out_line = "0.00 %s\r\n" % (round(PROF['hght'][row],2))
            snd_lines.append(out_line)

    # Add the bottom surface quantities
    snd_lines.extend(SFC_HEADER)
    for out_time, out_pres in zip(time_store, pres_store):
        ob_str = "%s0 %s %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\r\n" % (9999,
                                                             out_time,
                                                             out_pres, out_pres,
                                                             0.,
                                                             0., 0., 0.)
        snd_lines.append(ob_str)
        snd_lines.extend(SFC_VALS)
    snd_lines.append("\r\n\r\n")
    return "".join(snd_lines)

def read_bufkit(filename, date):
    """Read BUFKIT data stored on the IEM website

    Parameters
    ----------
    filename : string
        Path to locally stored BUFKIT file
    date: pandas.Timestamp
        Forecast valid time to search for

    Returns
    -------
    bufkit_data : Python dictionary
        Vertical profile of BUFKIT data
    """
    date_str = "%s%s%s/%s00" % (str(date.year).zfill(4)[-2:], str(date.month).zfill(2),
                                str(date.day).zfill(2), str(date.hour).zfill(2))

    with open(filename) as f: lines = f.readlines()
    head = lines[4][0:-14]
    tail = ' ' + date_str + ' \n'
    search_string = head + tail
    idx = lines.index(search_string) + 11

    bufkit_data = {
        'pressure': [],
        'temperature': [],
        'tmwc': [],
        'dewpoint': [],
        'thte': [],
        'windDir': [],
        'windSpeed': [],
        'omeg': [],
        'cfrl': [],
        'altitude': [],
        'apid': [],
    }
    for i in range(idx, len(lines), 2):
        bufkit_data['apid'].append(head[7:11])
        line_1 = lines[i].strip().split(' ')
        line_2 = lines[i+1].strip().split(' ')
        line_3 = lines[i+2]
        bufkit_data['pressure'].append(float(line_1[0]))
        bufkit_data['temperature'].append(float(line_1[1]))

        # Archived data has -9999s at varying levels generally above 150 mb. Not a
        # big deal since we won't have too much data above this level...
        bufkit_data['tmwc'].append(float(line_1[2]))
        bufkit_data['dewpoint'].append(float(line_1[3]))
        bufkit_data['thte'].append(float(line_1[4]))
        bufkit_data['windDir'].append(float(line_1[5]))
        bufkit_data['windSpeed'].append(float(line_1[6]))
        bufkit_data['omeg'].append(float(line_1[7]))
        bufkit_data['cfrl'].append(float(line_2[0]))
        bufkit_data['altitude'].append(float(line_2[1]))
        if line_3 == '\n': break
    return bufkit_data

def test_url(url):
    """Test for online file existence. Check the cheeky variable name!"""
    ru = requests.head(url)
    return ru.ok

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-a', '--archive', dest='archive', action='store_true',
                    help="Set for archive functionality.")
    args = ap.parse_args()
    args.data_path = DATA
    args.bufkit_path = SOUNDINGS

    if not args.archive:
        data_path = "%s/data_store/" % (args.data_path)
    else:
        data_path = "%s/data_store_archive/" % (args.data_path)

    # Read in the stored data files
    f_list = sorted(glob.glob(data_path + '/*'))
    prof_data = defaultdict(list)
    for f in f_list:
        try:
            with open(f, 'rb') as fp: temp = pickle.load(fp)
            prof_data[temp['apid']].append(temp)
        except:
            print("Error reading %s" % (f))

    # ----------------------------------------------------------------------------------
    # Pre-processing data for writing to BUFKIT files. Limit to those with surface
    # pressures above 900 mb and at least 20 levels
    # ----------------------------------------------------------------------------------
    for site in prof_data.keys():
        if (site in site_ids) and (prof_data[site][-1]['lat'] > 24) \
                              and (prof_data[site][-1]['lat'] < 51):
            timestamp("INFO")
            print("========================== %s ==========================" % (site))
            arrs = []
            for t in range(len(prof_data[site])):
                data = prof_data[site][t]
                if min_pressure_levels < len(data['pressure']) < max_pressure_levels \
                    and max(data['pressure'] > 900. ):
                    arrs.append(data)
            df = pd.DataFrame(arrs)
            df.drop_duplicates('dt', inplace=True)

            # BUFKIT expects evenly-spaced data in time.
            df['ob_time'] = df['dt']
            df['dt'] = df['dt'].dt.round('60min')

            # --------------------------------------------------------------------------
            # Downloading data here. Using unique times to only download data when
            # necessary. Step back in time while incrementing the forecast hour if files
            # not available.
            # --------------------------------------------------------------------------
            rap_data = {}
            unique_times = df.pivot_table(index=['dt']).index
            for t in range(len(unique_times)):
                fcst_offset = 1
                file_exists = False

                # Determine if this most likely file exists locally
                date = unique_times[t] - pd.to_timedelta(fcst_offset, unit='h')
                rap_datestring = "%s%s%s%s" % (str(date.year).zfill(4),
                                               str(date.month).zfill(2),
                                               str(date.day).zfill(2),
                                               str(date.hour).zfill(2))
                rap_file = "%s/RAP/%s_%s" % (args.bufkit_path, site, rap_datestring)

                # Otherwise look for it online
                if not os.path.exists(rap_file):
                    while not file_exists and fcst_offset <= 20:
                        date = unique_times[t] - pd.to_timedelta(fcst_offset, unit='h')
                        URL = "%s/%s/%s/%s/%s/%s/rap/rap_k%s.buf" % (bufkit_url,
                                                            str(date.year).zfill(4),
                                                            str(date.month).zfill(2),
                                                            str(date.day).zfill(2),
                                                            "bufkit",
                                                            str(date.hour).zfill(2),
                                                            site.lower())
                        file_exists = test_url(URL)
                        fcst_offset += 1
                    rap_datestring = "%s%s%s%s" % (str(date.year).zfill(4),
                                                   str(date.month).zfill(2),
                                                   str(date.day).zfill(2),
                                                   str(date.hour).zfill(2))
                    rap_file = "%s/RAP/%s_%s" % (args.bufkit_path, site, rap_datestring)
                    timestamp("INFO")
                    print("Downloading RAP BUFKIT data from IEM...")
                    print(URL)
                    urlreq.urlretrieve(URL, rap_file)
                rap_data[unique_times[t]] = read_bufkit(rap_file, unique_times[t])

            # --------------------------------------------------------------------------
            # Quality-control checks
            # Interpolate RAP sounding to the ACARS profile (in logp space).
            # --------------------------------------------------------------------------
            for index, row in df.iterrows():
                tmpc = interp_pres(row['pressure'], rap_data[row['dt']]['pressure'],
                                   rap_data[row['dt']]['temperature'])
                dwpc = interp_pres(row['pressure'], rap_data[row['dt']]['pressure'],
                                   rap_data[row['dt']]['dewpoint'])
                tmpc_deltas = np.fabs(tmpc - row['temperature'])

                # Replace ACARS data with BUFKIT data where QC checks fail
                df.at[index, 'temperature'] = np.where(tmpc_deltas > T_QC, tmpc,
                                                       row['temperature'])
                dwpc_deltas = np.fabs(dwpc - row['dewpoint'])
                df.at[index, 'dewpoint'] = np.where(dwpc_deltas > TD_QC, dwpc,
                                                    row['dewpoint'])
                #idx_1 = np.where(tmpc_deltas > T_QC)
                #idx_2 = np.where(dwpc_deltas > TD_QC)
                #qc_t = len(idx_1[0]) / len(row['pressure'])*100
                #qc_td = len(idx_2[0]) / len(row['pressure'])*100

                #if qc_t > 0. or qc_td > 0.:
                #    print(row['dt'])
                #    print("---->Replacing %s pct of temperature values" % (qc_t))
                #    print("---->Replacing %s pct of dewpoint values" % (qc_td))
                # Currently no QC checks for winds. May add this in later.
            # --------------------------------------------------------------------------
            # Output to file.
            #
            # For multiple profiles at a given hour, determine the mean (or median?).
            # The number of profiles can then be output at the top left of BUFKIT
            # --------------------------------------------------------------------------
            n_ensembles = df.pivot_table(index=['dt'], aggfunc='size').max()
            snd_str = ""
            #snd_str = to_bufkit(args.bufkit_path, df, rap_data=rap_data)
            #RAP = pd.DataFrame.from_dict(rap_data)
            #snd_str += to_bufkit(args.bufkit_path, RAP, rap_data=rap_data)
            #date_str = datetime.strftime(temp_df.iloc[i]['dt'], '%Y%m%d%H%M')
            for i in range(n_ensembles):
                temp_df = pd.DataFrame()
                for time in range(len(unique_times)):
                    prof = df.loc[(df.dt == unique_times[time])]
                    n_profs = len(prof)
                    temp_df = temp_df.append(prof.iloc[np.clip(i, 0, n_profs)-1])
                snd_str += to_bufkit(args.bufkit_path, temp_df, rap_data=rap_data,
                                     label=str(i))
            out_file = "%s/ACARSvapor_%s.buf" % (args.bufkit_path, site)
            with open(out_file, 'w') as fsnd: fsnd.write(snd_str)

# Don't run on module import
if __name__ == "__main__":
    main()
