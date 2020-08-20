"""download_acars.py

This script controls the downloading of ACARS/AMDAR data from MADIS and subsequent
output into Python-pickled dictionaries for storage on the local file system.

Much of this is based heavily on Tim Supinie's acars GitHub repository which can be
viewed here (https://github.com/tsupinie/acars). Several alterations and additions have
been made, however, to allow for downloading of the full VAPOR-AMDAR dataset (user-
restricted), as well as the addition of an archived "look-back" functionality (back to
2001).

The user must first acquire an account with the AMDAR/ESRL folks to access the
full suite of data (https://amdar.noaa.gov/FAQ.html). A username and password will then
need to be stored in a .env file within the main script directory for this script to
run properly.
"""

import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import pickle
import urllib.request as urlreq
from urllib.error import HTTPError
import base64
import os, subprocess, errno
import gzip
import glob
import re
import argparse
import pandas as pd
try:
    from BytesIO import BytesIO
except ImportError:
    from io import BytesIO
from stdatmos import std_atmosphere_pres
from utils.constants import ZEROCNK
from utils.cron_helper import *
from utils.configs import *

# No username and passwords here!! Create a .env file containing USERNAME and PASSWORD
# variables
from dotenv import load_dotenv
load_dotenv()
username = os.environ.get("USERNAME")
password = os.environ.get("PASSWORD")

stdatm = std_atmosphere_pres()
meta_path = os.path.dirname(__file__) or "."

class ACARSProfile(object):
    def __init__(self, apid, apcode, valid_dt, vapor, **prof_vars):
        self.apid = apid
        self.apcode = apcode
        self.dt = valid_dt
        self.prof_vars = prof_vars
        self.prof_vars['pressure'] = stdatm(self.prof_vars['altitude'])
        self.vapor = vapor
        self._sort()

    def apply_qc(self, meta_airport):
        # Check for a missing surface height
        if type(self.prof_vars['altitude'][0]) == type(np.ma.masked):
            print("Skipping profile: surface height at '%s' is qc'ed" % self.apid)
            return False

        # Check for distance from the claimed source airport
        ap_lat = meta_airport[self.apcode]['lat']
        ap_lon = meta_airport[self.apcode]['lon']
        ap_elev = meta_airport[self.apcode]['elev']
        ap_synop = meta_airport[self.apcode]['synop']

        self.lat = ap_lat
        self.lon = ap_lon
        self.elev = ap_elev
        self.synop = ap_synop

        dist = np.hypot(ap_lat - self.prof_vars['latitude'],
                        ap_lon - self.prof_vars['longitude'])
        if dist.min() > 2:
            print("Skipping profile: claims to be from '%s', but data are too far away"\
                  % self.apid)
            return False

        # This seems to do reasonably well filtering out just VAPOR AMDAR soundings.
        if self.vapor:
            if type(np.sum(self.prof_vars['dewpoint'])) == type(np.ma.masked):
                return False

        # Remove duplicate heights or out-of-order pressures
        bad_hghts = np.append(False, np.isclose(np.diff(self.prof_vars['altitude']), 0))
        bad_press = np.append(False, np.diff(self.prof_vars['pressure']) >= 0)
        keep = np.where(~(bad_hghts | bad_press))

        for var, val in self.prof_vars.items():
            self.prof_vars[var] = val[keep]

        # Check for number of data points
        if len(self.prof_vars['altitude']) < 3:
            return False
        return True

    def append(self, other):
        if self.apid != other.apid or self.dt != other.dt:
            raise ValueError("Profile id and time in append() must be the same")

        for var, vals in other.prof_vars.items():
            self.prof_vars[var] = np.ma.append(self.prof_vars[var], vals)

        self._sort()

    def _sort(self):
        sort_idxs = np.argsort(self.prof_vars['altitude'])
        for var, vals in self.prof_vars.items():
            self.prof_vars[var] = vals[sort_idxs]

def load_profiles(fname, supplemental, meta_airport, vapor=False):
    def _load_profiles(fname):
        load_vars = ['temperature', 'dewpoint', 'soundingSecs', 'sounding_airport_id',
                     'latitude', 'longitude', 'windSpeed', 'windDir', 'altitude',
                     'dataType', 'flight']
        nc = Dataset(fname)
        profile_data = {var: nc.variables[var][:] for var in load_vars}
        nc.close()

        # Split the profile arrays wherever the valid time changes. I guess this will
        # mess up if two adjacent profiles happen to be valid at the same time, but I'm
        # not sure that throws out too many good profiles.
        splits = np.where(np.diff(profile_data['soundingSecs']))[0] + 1
        profile_data_split = {var: np.split(prof_var, splits) for var, prof_var in \
                              profile_data.items()}
        profiles = []
        for vals in zip(*(profile_data_split[var] for var in load_vars)):
            val_dict = dict(zip(load_vars, vals))
            if type(val_dict['soundingSecs'][0]) == type(np.ma.masked):
                continue

            time = val_dict.pop('soundingSecs')
            apid = val_dict.pop('sounding_airport_id')
            prof_dt = epoch + timedelta(seconds=time[0])
            try:
                prof_id = meta_airport[apid[0]]['id']
            except KeyError:
                continue

            prof = ACARSProfile(prof_id, apid[0], prof_dt, vapor, **val_dict)
            profiles.append(prof)
        return profiles

    profiles = _load_profiles(fname)
    for fname in supplemental:
        supp_profiles = _load_profiles(fname)
        for supp_prof in supp_profiles:
            for prof in profiles:
                if supp_prof.apid == prof.apid and supp_prof.dt == prof.dt:
                    prof.append(supp_prof)
                    break
    return profiles

def load_meta(meta_fname=("%s/airport_info.dat" % meta_path)):
    """
    load_meta

    Load in our database of airport codes.
    """

    meta_cols = ['code', 'id', 'synop', 'lat', 'lon', 'elev', 'name']
    meta_types = {
        'code': int,
        'id': str,
        'synop': int,
        'lat': float,
        'lon': float,
        'elev': int,
        'name': str
    }
    meta_airport = {}
    with open(meta_fname) as fmeta:
        for line in fmeta:
            line_dict = {col: val for col, val in zip(meta_cols,
                                                      line.strip().split(None, 6)) }
            for col in line_dict.keys():
                line_dict[col] = meta_types[col](line_dict[col])

            code = line_dict.pop('code')
            meta_airport[code] = line_dict
    return meta_airport

def get_times(data_path):
    """This function figures out what times need to be downloaded from the server.
    Profiles are in files by hour, and the files are updated as new data comes in,
    so if we download a file once, we won't necessarily get all the data. So this
    also keeps track of which files were last accessed when. It does this by touching
    a file on disk (in marker_path) when it last accessed a file. If the file on the
    server is newer, the function knows that the file needs to be downloaded again.
    The marker files are deleted when they become older than the oldest file on the
    server.

    Parameters
    ----------
    data_path : string
        Full path to the netCDF directory on the local system

    Returns
    -------
    times_to_dl : list
        List of Python datetime objects indicating which times to download from MADIS
    """

    def touch(fname, times=None):
        with open(fname, 'a'):
            os.utime(fname, times)

    marker_path = "%s/markers" % data_path

    # Figure out the files and their update times from the MADIS server
    req = urlreq.Request(base_url)
    base64string = base64.b64encode(bytes('%s:%s' % (username, password), 'ascii'))
    req.add_header("Authorization", "Basic %s" % base64string.decode('utf-8'))
    result = urlreq.urlopen(req)
    txt = result.read().decode('utf-8')

    files = re.findall(">([\d]{8}_[\d]{4}).gz<", txt)
    update_times = re.findall("([\d]{2}-[\w]{3}-[\d]{4} [\d]{2}:[\d]{2})", txt)

    # If a MADIS server file is newer than the corresponding marker file, add
    # it to the list of times to download
    times_to_dl = []
    for file_time_str, update_time_str in zip(files, update_times):
        marker_fname = "%s/%s.txt" % (marker_path, file_time_str)

        update_time = datetime.strptime(update_time_str, '%d-%b-%Y %H:%M')
        file_time = datetime.strptime(file_time_str, '%Y%m%d_%H%M')

        if not os.path.exists(marker_fname) or (os.path.exists(marker_fname)
            and epoch + timedelta(seconds=os.path.getmtime(marker_fname))<update_time):
            #last_touch = epoch + timedelta(seconds=os.path.getmtime(marker_fname))
            #update_delta = (datetime.utcnow() - last_touch).total_seconds()
            touch(marker_fname)
            #if update_delta < update_window:
            times_to_dl.append(file_time)

    # Check for old marker files and delete them if they're older than the oldest file
    # on the MADIS server
    earliest_time = datetime.strptime(files[0], '%Y%m%d_%H%M')
    for fname in glob.glob("%s/*.txt" % marker_path):
        ftime = datetime.strptime(os.path.basename(fname), "%Y%m%d_%H%M.txt")
        if ftime < earliest_time:
            os.unlink(fname)

    return times_to_dl

def dl_profiles(data_path, dt, realtime=True):
    """Download netCDF files from the MADIS server. They're actually gzipped on the
    remote server, so unzip them in memory and write plain netCDF.

    Parameters
    ----------
    data_path : string
        Full path to the netCDF directory on the local system
    dt : datetime object
        Python datetime object indicating which time is being downloaded
    realtime : boolean [Optional: Default = True]
        Boolean flag indicating if this is a realtime request or not. An archived
        request will direct to a separate URL.

    Returns
    -------
    Downloads the requested AMDAR NetCDF files to data_path
    """

    if realtime:
        dl_url = base_url
    else:
        dl_url = "%s/%s/%s/%s/point/acars/netcdf" % (archive_url,
                                                    str(dt.year).zfill(4),
                                                    str(dt.month).zfill(2),
                                                    str(dt.day).zfill(2))
    url = "%s/%s.gz" % (dl_url, dt.strftime('%Y%m%d_%H%M'))

    # To send a username and password to restricted-access page, need to pass
    # to wget or curl.
    fname = "%s/%s.gz" % (data_path, dt.strftime("%Y%m%d_%H%M"))
    arg = "%s -q --user=%s --password=%s %s -O %s" % (WGET, username, password, url,
                                                      fname)

    # ----------------------------------------------------------------------------------
    # Actual downloading takes place here.
    # ----------------------------------------------------------------------------------
    os.system(arg)
    os.system("gunzip -f %s" % (fname))
    fname = "%s/%s" % (data_path, dt.strftime("%Y%m%d_%H%M"))
    return fname

def apply_granularity(profiles, granularity):
    unique_profiles = {}
    for profile in profiles:
        ap_code = profile.apid
        snd_time = (profile.dt - epoch).total_seconds()

        snd_time = granularity * np.round(snd_time / granularity)
        key = (snd_time, ap_code)
        if key not in unique_profiles or \
                                    len(unique_profiles[key].prof_vars['altitude']) < \
                                    len(profile.prof_vars['altitude']):
            unique_profiles[key] = profile
    return list(unique_profiles.values())

def main():
    ap = argparse.ArgumentParser()
    #ap.add_argument('-d', '--data-path', dest="data_path", default=NETCDF)
    #ap.add_argument('-o', '--output-path', dest="output_path", default=DATA)
    ap.add_argument('-a', '--archive-mode', dest="archive_mode", help="YYYY-MM-DD/HH")
    ap.add_argument('-v', '--vapor', dest='vapor', action='store_true',
                    help="Set for VAPOR soundings only. Otherwise, grab everything")
    args = ap.parse_args()
    meta = load_meta()
    time_gran = 600
    dts = get_times(args.data_path)
    args.data_path = DATA
    args.output_path = NETCDF

    # ----------------------------------------------------------------------------------
    # For archive "lookback" functionality...
    # ----------------------------------------------------------------------------------
    if args.archive_mode is None:
        current_time = datetime.utcnow()
        realtime = True
    else:
        current_time = datetime.strptime(args.archive_mode, "%Y-%m-%d/%H")
        dts = []
        for t in range(12):
            dts.append(current_time - timedelta(hours=t))
        realtime = False

    timestamp("INFO")
    print("Initializing new download request -----------------------------")

    if args.archive_mode is None:
        tree_path = "%s/data_store/" % (args.output_path)
    else:
        tree_path = "%s/data_store_archive/" % (args.output_path)
    try:
        os.mkdir(tree_path)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    dt_fnames = {}
    for dt in dts:
        dt_fnames[dt] = dl_profiles(args.data_path, dt, realtime=realtime)
        dt_prev = dt - timedelta(hours=1)
        dt_next = dt + timedelta(hours=1)
        if dt_prev not in dts:
            try:
                temp = dl_profiles(args.data_path, dt_next, realtime=realtime)
                if temp[1] == 0:
                    dt_fnames[dt_next] = temp[0]
            except:
                timestamp("WARN")
                print("No file for %s" % dt_prev.strftime("%d %b %Y %H%M UTC"))
        if dt_next not in dts:
            try:
                temp = dl_profiles(args.data_path, dt_next, realtime=realtime)
                if temp[1] == 0:
                    dt_fnames[dt_next] = temp[0]
            except:
                timestamp("WARN")
                print("No file for %s" % dt_next.strftime("%d %b %Y %H%M UTC"))
    for dt in dts:
        timestamp("INFO")
        print("New profiles for %s" % dt.strftime("%H%M UTC %d %b"))
        dt_prev = dt - timedelta(hours=1)
        dt_next = dt + timedelta(hours=1)
        supplemental = []
        if dt_prev in dt_fnames:
            supplemental.append(dt_fnames[dt_prev])
        if dt_next in dt_fnames:
            supplemental.append(dt_fnames[dt_next])
        profiles = load_profiles(dt_fnames[dt], supplemental, meta, args.vapor)
        profiles_qc = [ profile for profile in profiles if profile.apply_qc(meta) ]
        profiles_gran = apply_granularity(profiles_qc, time_gran)

        for profile in profiles_gran:
            data = {}
            data['pressure'] = np.array(profile.prof_vars['pressure']) / 100.
            data['altitude'] = np.array(profile.prof_vars['altitude'])
            data['temperature'] = np.array(profile.prof_vars['temperature']) - ZEROCNK
            data['dewpoint'] = np.array(profile.prof_vars['dewpoint']) - ZEROCNK
            data['windDir'] = np.array(profile.prof_vars['windDir'])
            data['windSpeed'] = np.array(profile.prof_vars['windSpeed'])
            data['apid'] = profile.apid
            data['synop'] = profile.synop
            data['dt'] = profile.dt
            data['elev'] = profile.elev
            data['lat'] = profile.lat
            data['lon'] = profile.lon
            if len(profile.prof_vars['flight'][0]) > 0:
                try:
                    flight_num = b''.join(profile.prof_vars['flight'][0])
                    flight_num = int(flight_num.decode('utf-8'))
                except:
                    flight_num = int(profile.prof_vars['flight'][0])
            else:
                flight_num = -9999
            data['flight'] = flight_num

            # Construct the file name (using the time granularity)
            dt_sec = round((profile.dt - epoch).total_seconds() / \
                            time_gran) * time_gran
            dt_round = epoch + timedelta(seconds=dt_sec)
            fname = "%s/%s_%s" % (tree_path, profile.apid,
                                     dt_round.strftime("%Y%m%d%H%M"))

            #exist_size = os.path.getsize(fname) if os.path.exists(fname) else 0
            # Check for a smaller file that already exists. The idea is to avoid
            # writing over a "good" file with a "bad" file, where file size is used
            # as a proxy for "goodness". This may not be the best proxy, though.
            #if len(profile.prof_vars['pressure']) < exist_size:
            #    print("Skipping profile %s" % fname)
            #else:
            with open(fname, 'wb') as fsnd: pickle.dump(data, fsnd)

        # Loop through the temporary data_store directory and move files older than
        # purge_hours to the main storage tree
        purge_dt = current_time - timedelta(hours=purge_hours)
        f_list = glob.glob(tree_path + '/*')
        for file_ in f_list:
            file_dt = datetime.strptime(file_[-12:], '%Y%m%d%H%M')
            store_path = "%s/%s" % (args.output_path, file_dt.strftime("%Y/%m/%d/%H"))
            short_fname = file_[file_.rfind('/')+1:]
            try:
                os.makedirs(store_path)
            except OSError:
                pass
            if file_dt < purge_dt:
                timestamp("INFO")
                print("    Moving %s to %s" % (short_fname, store_path))
                print(file_, store_path + '/' + short_fname)
                os.rename(file_, store_path + '/' + short_fname)

if __name__ == "__main__":
    main()
