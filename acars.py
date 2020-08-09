#from __future__ import print_function
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

# --------------------------------------------------------------------------------------
#
# Some configuration information. Move this to a separate module file
#
# --------------------------------------------------------------------------------------
bufkit_url = 'https://mtarchive.geol.iastate.edu'
#site_ids = [
#    'MDW', 'MKE'
#]
site_ids = ['MDW']

# QC thresholds (degrees C)
T_QC = 5.
TD_QC = 15.

def to_bufkit(output_path, data, rap_data=None, label='NONE'):
    """Writes output to BUFKIT-readble .buf file.

    BUFKIT is extremely unforgiving when it comes to formatting the input files. Many of
    these odds and ends aren't readily apparent, although you'll know you've hit them
    when the entire program crashes with various cryptic VB-based error messages!

    Various parameter constraints I've been able to deduce include (BUFKIT-19):
        :: At least 5 individual time steps (0,1,2,3,4)
        :: Minimum number of pressure levels @ individual timestep = 20 (inclusive)
        :: Maximum number of pressure levels @ individual timestep = 67 (inclusive)
        :: Minimum TOTAL number of pressure levels in a file = 201
        :: Top pressure level must be < 480 mb (divide by zero error otherwise)
        :: ?? Maximum number of 'forecast' times is (about?) 120
        :: ?? 'Forecast' time steps only as small as every 15 minutes and only ever
           quarter hour
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
    Writes a .buf file to \output_path\.
    """
    min_levs = 20
    max_levs = 67
    min_times = 5

    ID = data.iloc[0]['apid']
    n_times = len(data)
    time_store = []
    pres_store = []
    date_last = pd.to_datetime(data.iloc[-1]['dt']).round('30min')
    knt = 0
    total_lev_knt = 0

    snd_lines = ["%s\r\n" % (label)]
    snd_lines.extend(FILE_HEADER)
    # ----------------------------------------------------------------------------------
    # Data for outputting
    # ----------------------------------------------------------------------------------
    for t in range(n_times):
        #date_ob = date_last - timedelta(hours=n_times-t-1)
        date_ob = date_last - timedelta(minutes=30*(n_times-t-1))
        date = datetime.strftime(date_ob, "%Y%m%d%H%M")
        time_store.append(datetime.strftime(date_ob, "%Y%m%d/%H%M"))
        levs = np.array(data.iloc[t]['pressure'])
        t_out = np.array(data.iloc[t]['temperature'])
        td_out = np.array(data.iloc[t]['dewpoint'])
        wdir_out = np.clip(np.array(data.iloc[t]['windDir']), 1, 359)
        wspd_out = np.clip(np.array(data.iloc[t]['windSpeed']), 0.5, 1000)
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
        for row in range(len(t_out)):
            if lev_knt >= max_levs: break
            thetae = thermo.equivalent_potential_temperature(levs[row],
                                                             t_out[row],
                                                             td_out[row])
            out_line = "%s %s %s %s %s %s %s %s\r\n" % (round(levs[row],2),
                                                        round(t_out[row],2),
                                                        round((t_out[row]+td_out[row])/2.,2),
                                                        round(td_out[row],2),
                                                        round(thetae,2),
                                                        round(wdir_out[row],2),
                                                        round(wspd_out[row],2),
                                                        -9999.)
            snd_lines.append(out_line)
            out_line = "0.00 %s\r\n" % (round(hght_out[row],2))
            snd_lines.append(out_line)
            lev_knt += 1
            total_lev_knt += 1
        pres_store.append(levs[0])

        # Use the RAP data to extend our profile up to 67 pressure levels. If not
        # available, extend using some random data (see utils.bufkit_helper)
        line_knt = 0
        counter = 0
        if rap_data is not None:
            RAP = rap_data[data.iloc[t]['dt']]
            top_pressure = levs[row]
            rap_idx = np.where(RAP['pressure'] < top_pressure)[0]
            while (lev_knt < max_levs) and (counter < rap_idx[-1]):
                counter = rap_idx[0] + line_knt
                out_line = "%s %s %s %s %s %s %s %s\r\n" % (RAP['pressure'][counter],
                                                            RAP['temperature'][counter],
                                                            RAP['tmwc'][counter],
                                                            RAP['dewpoint'][counter],
                                                            RAP['thte'][counter],
                                                            RAP['windDir'][counter],
                                                            RAP['windSpeed'][counter],
                                                            -9999.)
                snd_lines.append(out_line)
                out_line = "0.00 %s\r\n" % (RAP['altitude'][counter])
                snd_lines.append(out_line)
                lev_knt += 1
                line_knt += 1
        else:
            while (lev_knt < max_levs):
                snd_lines.extend(EXTRA_LINES[line_knt:line_knt+2])
                line_knt += 2
                lev_knt += 1
        knt += 1

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

    if not is_bufkitdata_bad(snd_lines):
        snd_str = "".join(snd_lines)
        return snd_str
    else:
        return ""

def is_bufkitdata_bad(data):
    """Runs various checks to **TRY** to ensure this file won't crash BUFKIT. I'm not
    sure we'll ever be 100% sure BUFKIT won't die, however!
    """
    ret = False
    indices = [index for index, value in enumerate(data) if value == \
               'PRES TMPC TMWC DWPC THTE DRCT SKNT OMEG\r\n']
    sum = 0
    if len(indices) >= 5:
        for i in range(len(indices)-1):
            start = indices[i]+2
            end = indices[i+1]-12
            num_plevs = int(((end - start) / 2) + 1)
            sum = sum + num_plevs
        idx = data.index('STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS\r\n')
        num_plevs = int((idx-2-indices[i+1])/2)
        sum = sum + num_plevs
    else:
        ret = True
    if sum < 201:
        ret = True
    return ret

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
        'altitude': []
    }
    for i in range(idx, len(lines), 2):
        line_1 = lines[i].strip().split(' ')
        line_2 = lines[i+1].strip().split(' ')
        line_3 = lines[i+2]
        bufkit_data['pressure'].append(float(line_1[0]))
        bufkit_data['temperature'].append(np.clip(float(line_1[1]), -150, 70))
        bufkit_data['tmwc'].append(float(line_1[2]))

        # Archived data has -9999s above about 150 mb. Not a big deal since we won't
        # have too much data above this level.
        bufkit_data['dewpoint'].append(np.clip(float(line_1[3]), -100, 40))
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
    ap.add_argument('-d', '--data-path', dest="data_path", default="/Users/leecarlaw/scripts/acars/data")
    ap.add_argument('-b', '--bufkit-path', dest="bufkit_path", default="/Users/leecarlaw/scripts/acars/soundings")
    args = ap.parse_args()
    #data_path = "%s/data_store_archive/" % (args.data_path)
    data_path = "%s/data_store/" % (args.data_path)

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
            print("============ %s =============" % (site))
            arrs = []
            for t in range(len(prof_data[site])):
                data = prof_data[site][t]
                if max(data['pressure']) > 900 and (min(data['pressure']) < 450) \
                                               and (len(data['pressure']) > 20):
                    arrs.append(data)
            df = pd.DataFrame(arrs)
            df.drop_duplicates('dt', inplace=True)
            df['ob_time'] = df['dt']
            df['dt'] = df['dt'].dt.round('30min')

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
                if not os.path.exists(rap_file):
                    print(URL)
                    urlreq.urlretrieve(URL, rap_file)
                rap_data[unique_times[t]] = read_bufkit(rap_file, unique_times[t])

            # --------------------------------------------------------------------------
            # Quality-control checks
            #
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

            # --------------------------------------------------------------------------
            # Output to file.
            #
            # For multiple profiles at a given hour, determine the mean (or median?).
            # The number of profiles can then be output at the top left of BUFKIT
            # --------------------------------------------------------------------------
            n_ensembles = df.pivot_table(index=['dt'], aggfunc='size').max()
            snd_str = ""
            snd_str = to_bufkit(args.bufkit_path, df, rap_data=rap_data)
            '''
            for i in range(n_ensembles):
                temp_df = pd.DataFrame()
                for time in range(len(unique_times)):
                    prof = df.loc[(df.dt == unique_times[time])]
                    n_profs = len(prof)
                    temp_df = temp_df.append(prof.iloc[np.clip(i, 0, n_profs)-1])
                snd_str += to_bufkit(args.bufkit_path, temp_df, rap_data=bufkit_data,
                                     label=str(i+1))
            '''
            out_file = "%s/ACARSvapor_%s.buf" % (args.bufkit_path, site)
            with open(out_file, 'w') as fsnd: fsnd.write(snd_str)

if __name__ == "__main__":
    main()
