import numpy as np
import os
from datetime import datetime, timedelta
import metpy.calc as mpcalc
from metpy.units import units

import urllib.request as urlreq
from urllib.error import HTTPError

from utils.configs import bufkit_url

def get_iem_bufkit(site, target_dt, model='rap', path='./'):
    """
    Download and parse BUFKIT data from the IEM archive

    Parameters
    ----------
    site : string
        Site ID name. Example: 'ORD'
    target_dt : datetime
        Valid forecast time
    model : string [Optional. Default 'rap']
        rap or hrrr
    path : string [Optional. Default './']
        Location on the filesystem to store downoaded text files

    Returns
    -------
    data : dictionary
        Parsed BUFKIT data for a single forecast time (target_dt)

    """

    cycle_datestring = (target_dt - timedelta(hours=0)).strftime("%Y%m%d%H")
    filename = "%s/%s_%s_%s.buf" % (path, model, site, cycle_datestring)
    if not os.path.exists(filename):
        URL = "%s/%s/%s/%s/bufkit/%s/%s/%s_k%s.buf" % (bufkit_url,
                                                       cycle_datestring[0:4],
                                                       cycle_datestring[4:6],
                                                       cycle_datestring[6:8],
                                                       cycle_datestring[8:10],
                                                       model, model,
                                                       site.lower())
        print("DOWNLOADING", (URL))
        urlreq.urlretrieve(URL, filename)
    data = read_bufkit(filename, target_dt)
    return data

def mag(u, v):
    """Compute the magnitude of a vector from its components

    Parameters
    ----------
    u : number, array_like
        U-component of the wind
    v : number, array_like
        V-component of the wind

    Returns
    -------
    mag : number, array_like
        The magnitude of the vector (units are the same as input)
    """
    return np.sqrt(u**2 + v**2)

def comp2vec(u, v):
    '''
    Convert U, V components into direction and magnitude

    Parameters
    ----------
    u : number, array_like
        U-component of the wind
    v : number, array_like
        V-component of the wind

    Returns
    -------
    wdir : number, array_like (same as input)
        Angle in meteorological degrees
    wspd : number, array_like (same as input)
        Magnitudes of wind vector (input units == output units)

    '''

    u = np.array(u).astype(np.float64)
    v = np.array(v).astype(np.float64)
    wdir = np.degrees(np.arctan2(-u, -v))

    if wdir.shape:
        wdir[wdir < 0] += 360
        wdir[np.fabs(wdir) < 1e-10] = 0.
    else:
        if wdir < 0:
            wdir += 360
        if np.fabs(wdir) < 1e-10:
            wdir = 0.
    return wdir, mag(u, v)

def read_bufkit(filename, date):
    """
    Read BUFKIT data stored on the IEM website. Currently searches for data valid at a
    certain time, specified by the date argument.

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
        'hght': [],
        'apid': [],
        'u': [],
        'v': []
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
        bufkit_data['hght'].append(float(line_2[1]))
        if line_3 == '\n': break

    U, V = mpcalc.wind_components(bufkit_data['windSpeed']*units('m/s'),
                                  bufkit_data['windDir']*units.deg)
    bufkit_data['u'] = U.magnitude
    bufkit_data['v'] = V.magnitude
    return bufkit_data

def test_url(url):
    """
    Test for online file existence.

    Parameters
    ----------
    url : string
        URL we're testing for

    """
    ru = requests.head(url)
    return ru.ok
