#from __future__ import print_function
import numpy as np
from netCDF4 import Dataset
from datetime import datetime, timedelta
import pickle
import urllib.request as urlreq
from urllib.error import HTTPError
import base64
import re
import os, subprocess, errno
import gzip
import glob
import argparse
import math
import time
import pandas as pd
from random import choice
from string import ascii_uppercase
try:
    from BytesIO import BytesIO
except ImportError:
    from io import BytesIO
from stdatmos import std_atmosphere_pres
from utils.constants import ZEROCNK

# No username and passwords here!!
from dotenv import load_dotenv
load_dotenv()
password = os.environ.get("PASSWORD")
username = os.environ.get("USERNAME")

_purge_hours = 24
update_window = 6 * 3600
_epoch = datetime(1970, 1, 1, 0)
missing = np.nan
_base_url = "https://madis-data.ncep.noaa.gov/madisNoaa/data/point/acars/netcdf"
_archive_url = "https://madis-data.ncep.noaa.gov/madisNoaa/data/archive"
_stdatm = std_atmosphere_pres()
_meta_path = os.path.dirname(__file__)
if _meta_path == "":
    _meta_path = "."

class ACARSProfile(object):
    def __init__(self, apid, apcode, valid_dt, **prof_vars):
        self.apid = apid
        self.apcode = apcode
        self.dt = valid_dt
        self.prof_vars = prof_vars

        self.prof_vars['pressure'] = _stdatm(self.prof_vars['altitude'])

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

        # For the time being, we're really only interested in VAPOR soundings.
        # These are of dataType 3 and 4 (vapor measuring and edr, rh, ice measuring
        # aircraft
        #if self.prof_vars['dataType'][0] not in [3.,6.]:
        #    return False
        #if self.prof_vars['dataType'][0] != 3:
        #    return False

        # This seems to do reasonably well filtering out just VAPOR AMDAR soundings.
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

def load_profiles(fname, supplemental, meta_airport):
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
            prof_dt = _epoch + timedelta(seconds=time[0])
            try:
                prof_id = meta_airport[apid[0]]['id']
            except KeyError:
                continue

            prof = ACARSProfile(prof_id, apid[0], prof_dt, **val_dict)
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

def load_meta(meta_fname=("%s/airport_info.dat" % _meta_path)):
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
    """
    get_times

    This function figures out what times need to be downloaded from the server.
    Profiles are in files by hour, and the files are updated as new data comes in,
    so if we download a file once, we won't necessarily get all the data. So this
    also keeps track of which files were last accessed when. It does this by touching
    a file on disk (in marker_path) when it last accessed a file. If the file on the
    server is newer, the function knows that the file needs to be downloaded again.
    The marker files are deleted when they become older than the oldest file on the
    server.
    """
    def touch(fname, times=None):
        with open(fname, 'a'):
            os.utime(fname, times)

    marker_path = "%s/markers" % data_path

    # Figure out the files and their update times from the MADIS server
    req = urlreq.Request(_base_url)
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
            and _epoch + timedelta(seconds=os.path.getmtime(marker_fname))<update_time):
            #last_touch = _epoch + timedelta(seconds=os.path.getmtime(marker_fname))
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


def dl_profiles(path, dt, realtime=True):
    """
    dl_profiles

    Download netCDF files from the MADIS server. They're actually gzipped on the remote
    server, so unzip them in memory and write plain netCDF.
    """
    if realtime:
        dl_url = _base_url
    else:
        dl_url = "%s/%s/%s/%s/point/acars/netcdf" % (_archive_url,
                                                    str(dt.year).zfill(4),
                                                    str(dt.month).zfill(2),
                                                    str(dt.day).zfill(2))
    url = "%s/%s.gz" % (dl_url, dt.strftime('%Y%m%d_%H%M'))

    # To send a username and password to restricted-access page, need to pass
    # to wget or curl
    fname = "%s/%s.gz" % (path, dt.strftime("%Y%m%d_%H%M"))
    arg = "wget -q --user=%s --password=%s %s -O %s" % (username, password, url, fname)

    # ----------------------------------------------------------------------------------
    # Actual downloading takes place here.
    # ----------------------------------------------------------------------------------
    os.system(arg)
    os.system("gunzip -f %s" % (fname))
    time.sleep(0.1)
    fname = "%s/%s" % (path, dt.strftime("%Y%m%d_%H%M"))
    return fname

def apply_granularity(profiles, granularity):
    unique_profiles = {}
    for profile in profiles:
        ap_code = profile.apid
        snd_time = (profile.dt - _epoch).total_seconds()

        snd_time = granularity * np.round(snd_time / granularity)
        key = (snd_time, ap_code)
        if key not in unique_profiles or \
                                    len(unique_profiles[key].prof_vars['altitude']) < \
                                    len(profile.prof_vars['altitude']):
            unique_profiles[key] = profile
    return list(unique_profiles.values())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-d', '--data-path', dest="data_path", default="/Users/leecarlaw/scripts/acars/netcdf")
    ap.add_argument('-o', '--output-path', dest="output_path", default="/Users/leecarlaw/scripts/acars/data")
    #ap.add_argument('-v', '--vapor', dest='vapor', action='store_true')
    args = ap.parse_args()
    meta = load_meta()
    time_gran = 600

    print("Retrieving profile times ...")
    dts = get_times(args.data_path)

    # ----------------------------------------------------------------------------------
    # For archive "lookback" functionality...
    # ----------------------------------------------------------------------------------
    #end_time = datetime(year=2020, month=8, day=8, hour=11)
    current_time = datetime.utcnow()
    #current_time = datetime(year=2019, month=10, day=21, hour=2)
    #dts = []
    #for t in range(12):
    #    dts.append(current_time - timedelta(hours=t))
    #realtime = False
    realtime = True

    if realtime:
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
    print("Downloading profiles ...")
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
                print("No file for %s" % dt_prev.strftime("%d %b %Y %H%M UTC"))
        if dt_next not in dts:
            try:
                temp = dl_profiles(args.data_path, dt_next, realtime=realtime)
                if temp[1] == 0:
                    dt_fnames[dt_next] = temp[0]
            except:
                print("No file for %s" % dt_next.strftime("%d %b %Y %H%M UTC"))
    for dt in dts:
        print("New profiles for %s" % dt.strftime("%H%M UTC %d %b"))
        dt_prev = dt - timedelta(hours=1)
        dt_next = dt + timedelta(hours=1)
        print("Parsing profiles ...")
        supplemental = []
        if dt_prev in dt_fnames:
            supplemental.append(dt_fnames[dt_prev])
        if dt_next in dt_fnames:
            supplemental.append(dt_fnames[dt_next])
        profiles = load_profiles(dt_fnames[dt], supplemental, meta)
        profiles_qc = [ profile for profile in profiles if profile.apply_qc(meta) ]
        profiles_gran = apply_granularity(profiles_qc, time_gran)
        print("Dumping files ...")
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
            dt_sec = round((profile.dt - _epoch).total_seconds() / \
                            time_gran) * time_gran
            dt_round = _epoch + timedelta(seconds=dt_sec)
            #random_str = ''.join(choice(ascii_uppercase) for i in range(6))
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
            #with open(fname + '.txt', 'w') as fsnd: print(data, file=fsnd)

        # Loop through the temporary data_store directory and move files older than
        # _purge_hours to the main storage tree
        purge_dt = current_time - timedelta(hours=_purge_hours)
        f_list = glob.glob(tree_path + '/*')
        for file_ in f_list:
            #idx = file_.rfind('_')
            #file_dt = datetime.strptime(file_[idx-12:idx], '%Y%m%d%H%M')
            file_dt = datetime.strptime(file_[-12:], '%Y%m%d%H%M')
            store_path = "%s/%s" % (args.output_path, file_dt.strftime("%Y/%m/%d/%H"))
            short_fname = file_[file_.rfind('/')+1:]
            try:
                os.makedirs(store_path)
            except OSError:
                pass
            if file_dt < purge_dt:
                print("    Moving %s to %s" % (short_fname, store_path))
                print(file_, store_path + '/' + short_fname)
                os.rename(file_, store_path + '/' + short_fname)

if __name__ == "__main__":
    main()
