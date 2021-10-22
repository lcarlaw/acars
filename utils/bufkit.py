import pandas as pd
import numpy as np
from datetime import datetime
from collections import defaultdict
import os

import metpy.calc as mpcalc
from metpy.units import units

import utils.thermo as thermo
from interpolate.tools import interp_pres

curr_path = os.path.dirname(__file__) or "."
# --------------------------------------------------------------------------------------
#
# Standard header and other BUFKIT-related file data that we don't need floating around
# in the main codes.
#
# --------------------------------------------------------------------------------------
#df = pd.read_csv(curr_path + '/L62.txt', delim_whitespace=True)
#a_k = df['a'] / 100.
#b_k = df['b']
FILE_HEADER = [
    "SNPARM = PRES;TMPC;TMWC;DWPC;THTE;DRCT;SKNT;OMEG;CFRL;HGHT\r\n",
    "STNPRM = SHOW;LIFT;SWET;KINX;LCLP;PWAT;TOTL;CAPE;LCLT;CINS;EQLV;LFCT;BRCH\r\n"
]

SFC_HEADER = [
    "STN YYMMDD/HHMM PMSL PRES SKTC STC1 SNFL WTNS\r\n",
    "P01M C01M STC2 LCLD MCLD HCLD\r\n",
    "SNRA UWND VWND R01M BFGR T2MS\r\n",
    "Q2MS WXTS WXTP WXTZ WXTR USTM\r\n",
    "VSTM HLCY SLLH WSYM CDBP VSBK\r\n",
    "TD2M\r\n"
]

SFC_VALS = [
    "%s 0.00 0.00 0.00 0.00 0.00\r\n" % (0.00),
    "-9999.00 0.00 0.00 -9999.00 -9999.00 0.00\r\n",
    "0.00 0.00 0.00 0.00 0.00 0.00\r\n",
    "0.00 0.00 -0.00 999.00 -9999.00 0.00\r\n",
    "0.00\r\n"
]

# In the cases where we need additional "faked" data to avoid crashing BUFKIT. I can't
# imagine a scenario where a jet provides data above 16,000 m. This data is only used
# if RAP data is not supplied to augment the vertical ACARS profile.
EXTRA_LINES = [
    "107.80 -70.06 -70.14 -102.17 383.53 117.26 7.22 -9999.0\r\n",
    "0.00 16244.36\r\n",
    "94.90 -71.86 -71.93 -102.76 394.21 90.70 15.94 -9999.0\r\n",
    "0.00 16998.72\r\n",
    "83.20 -73.36 -73.43 -103.36 406.25 101.31 23.79 -9999.0\r\n",
    "0.00 17771.12\r\n",
    "72.90 -71.56 -71.66 -103.96 425.67 111.55 24.87 -9999.0\r\n",
    "0.00 18547.52\r\n",
    "63.70 -66.66 -66.87 -104.57 453.13 114.48 23.92 -9999.0\r\n",
    "0.00 19353.28\r\n",
    "55.50 -57.96 -58.63 -105.19 491.16 90.00 13.61 -9999.0\r\n",
    "0.00 20294.94\r\n",
    "48.20 -57.16 -57.99 -105.81 513.23 96.17 14.47 -9999.0\r\n",
    "0.00 21184.93\r\n",
    "41.80 -55.06 -56.25 -106.44 539.73 93.18 14.02 -9999.0\r\n",
    "0.00 22090.05\r\n",
    "36.10 -54.16 -55.64 -107.07 565.11 71.82 13.71 -9999.0\r\n",
    "0.00 23027.91\r\n",
    "31.20 -52.16 -54.19 -107.70 594.51 74.51 20.37 -9999.0\r\n",
    "0.00 23967.30\r\n",
    "26.80 -51.06 -53.60 -108.35 623.96 88.51 22.36 -9999.0\r\n",
    "0.00 24953.14\r\n",
    "22.90 -50.16 -53.25 -109.01 655.25 91.85 24.12 -9999.0\r\n",
    "0.00 25977.62\r\n",
    "19.50 -48.06 -52.17 -109.69 692.47 85.66 28.27 -9999.0\r\n",
    "0.00 27031.67\r\n",
    "16.50 -45.26 -50.89 -110.38 735.32 87.32 33.28 -9999.0\r\n",
    "0.00 28139.24\r\n"
]

def to_bufkit(out_path, data, meta, rap_data=None, label=None):
    """
    Create BUFKIT-readable text.

    BUFKIT is extremely unforgiving when it comes to formatting the input files. Many
    of these odds and ends aren't readily apparent, although you'll know you've hit
    them when the entire program crashes with various cryptic VB-based error messages!

    Various parameter constraints I've been able to deduce include (BUFKIT-19):
        [x] :: At least 5 individual time steps (0,1,2,3,4)
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

    Return
    -------
    Returns BUFKIT data as a string.

    """

    min_levs = 20
    max_levs = 67
    min_times = 5
    ID = meta['id']
    n_times = len(data)
    time_store = []
    pres_store = []
    date_last = pd.to_datetime(data.iloc[-1]['datetime'])
    knt = 0
    total_lev_knt = 0

    snd_lines = ["%s\r\n" % (label)]
    snd_lines.extend(FILE_HEADER)

    # ----------------------------------------------------------------------------------
    # Data for outputting
    # ----------------------------------------------------------------------------------
    for t in range(n_times):
        date_ob = data.iloc[t]['datetime']
        date = datetime.strftime(date_ob, "%Y%m%d%H%M")
        time_store.append(datetime.strftime(date_ob, "%Y%m%d/%H%M"))
        levs = np.array(data.iloc[t]['pressure'])
        t_out = np.array(data.iloc[t]['temperature'])
        td_out = np.array(data.iloc[t]['dewpoint'])
        wdir_out = np.array(data.iloc[t]['windDir'])
        wspd_out = np.array(data.iloc[t]['windSpeed'])
        hght_out = np.array(data.iloc[t]['hght'])

        # A few last-minute sanity checks
        td_out = np.clip(td_out, -100, t_out-0.15)

        out_time = "%s%s%s/%s%s"%(date[2:4],date[4:6],date[6:8],date[8:10],date[10:12])
        snd_lines.extend([
            "\r\n",
            "STID = K%s STNM = %s0 TIME = %s\r\n" % (ID, 9999, out_time),
            "SLAT = %s SLON = %s SELV = %s\r\n" % (meta['lat'],
                                                   meta['lon'],
                                                   meta['elev']),
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
            if lev_knt >= max_levs: break
            if levs[row] > 1:
                prof['pres'].append(levs[row])
                prof['tmpc'].append(t_out[row])
                prof['dwpc'].append(td_out[row])
                prof['wdir'].append(wdir_out[row])
                prof['wspd'].append(wspd_out[row])
                prof['hght'].append(hght_out[row])
                prof['tmwc'].append(-9999.0)
        p_sfc = levs[0]
        pres_store.append(p_sfc)

        # Derived parameter calculations
        prof['thte'] = mpcalc.equivalent_potential_temperature(prof['pres']*units.hPa,
                                                   prof['tmpc']*units.degC,
                                                   prof['dwpc']*units.degC).magnitude

        # Since this iterative calculation is so slow, I've left the wetbulb temperature
        # out of the ACARS sounding for the time being.
        #prof['tmwc'] = mpcalc.wet_bulb_temperature(prof['pres']*units.hPa,
        #                                           prof['tmpc']*units.degC,
        #                                           prof['dwpc']*units.degC).magnitude

        # Use the RAP data to extend our profile up to a maximum of 67 pressure levels.
        # If not available, extend using some random data (see utils.bufkit_helper).
        if rap_data is not None:
            num_extra_levs = max_levs - lev_knt
            RAP = rap_data[data.iloc[t]['datetime']]
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
                prof['hght'].append(RAP['hght'][counter])
                lev_knt += 1
                counter += 1

        knt += 1

        '''
        # In order to ensure consistent pressure levels which BUFKIT seems to desire,
        # derive these using a hybrid sigma vertical pressure coordinate system based
        # on the ECMWF L62 grid.
        p_levs = []
        for k in range(len(b_k)):
            p = a_k.iloc[k] + b_k.iloc[k] * p_sfc
            p_levs.append(p)
        p_levs = p_levs[::-1]
        print(p_levs)

        prof = {}
        for var in prof.keys():
            prof[var] = interp_pres(p_levs, prof['pres'], prof[var])
            prof[var] = np.where(prof[var] < -1000., -9999., prof[var])
        '''

        # Appending to the output array
        for row in range(len(prof['pres'])):
            out_line = "%s %s %s %s %s %s %s %s\r\n" % (round(prof['pres'][row],2),
                                                        round(prof['tmpc'][row],2),
                                                        round(prof['tmwc'][row],2),
                                                        round(prof['dwpc'][row],2),
                                                        round(prof['thte'][row],2),
                                                        round(prof['wdir'][row],2),
                                                        round(prof['wspd'][row],2),
                                                        0.)
            snd_lines.append(out_line)
            out_line = "0.00 %s\r\n" % (round(prof['hght'][row],2))
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
    return snd_lines


def to_bufkit_test(data, meta):
    min_levs = 20
    max_levs = 67
    min_times = 5
    ID = meta['id']
    n_times = len(data)
    time_store = []
    pres_store = []
    date_last = pd.to_datetime(data.iloc[-1]['datetime'])
    knt = 0
    total_lev_knt = 0

    snd_lines = ["%s\r\n" % ("TEST")]
    snd_lines.extend(FILE_HEADER)

    # ----------------------------------------------------------------------------------
    # Data for outputting
    # ----------------------------------------------------------------------------------

    date_ob = data['datetime'].iloc[0]
    date = datetime.strftime(date_ob, "%Y%m%d%H%M")
    time_store.append(datetime.strftime(date_ob, "%Y%m%d/%H%M"))
    levs = np.array(data['pressure'])
    t_out = np.array(data['temperature'])
    td_out = np.array(data['dewpoint'])
    wdir_out = np.array(data['windDir'])
    wspd_out = np.array(data['windSpeed'])
    hght_out = np.array(data['hght'])

    td_out = np.clip(td_out, -99., t_out)

    out_time = "%s%s%s/%s%s"%(date[2:4],date[4:6],date[6:8],date[8:10],date[10:12])
    snd_lines.extend([
        "\r\n",
        "STID = K%s STNM = %s0 TIME = %s\r\n" % (ID, 9999, out_time),
        "SLAT = %s SLON = %s SELV = %s\r\n" % (meta['lat'],
                                               meta['lon'],
                                               meta['elev']),
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
        if lev_knt >= max_levs: break
        if levs[row] > 1:
            prof['pres'].append(levs[row])
            prof['tmpc'].append(t_out[row])
            prof['dwpc'].append(td_out[row])
            prof['wdir'].append(wdir_out[row])
            prof['wspd'].append(wspd_out[row])
            prof['hght'].append(hght_out[row])
            prof['tmwc'].append(-9999.0)
    p_sfc = levs[0]
    pres_store.append(p_sfc)

    # Derived parameter calculations
    prof['thte'] = mpcalc.equivalent_potential_temperature(prof['pres']*units.hPa,
                                               prof['tmpc']*units.degC,
                                               prof['dwpc']*units.degC).magnitude
    prof['thte'] = np.nan_to_num(prof['thte'], nan=-9999.0)
    knt += 1

    # Appending to the output array
    for row in range(len(prof['pres'])):
        out_line = "%s %s %s %s %s %s %s %s\r\n" % (round(prof['pres'][row],2),
                                                    round(prof['tmpc'][row],2),
                                                    round(prof['tmwc'][row],2),
                                                    round(prof['dwpc'][row],2),
                                                    round(prof['thte'][row],2),
                                                    round(prof['wdir'][row],2),
                                                    round(prof['wspd'][row],2),
                                                    0.)
        snd_lines.append(out_line)
        out_line = "0.00 %s\r\n" % (round(prof['hght'][row],2))
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
    return snd_lines
