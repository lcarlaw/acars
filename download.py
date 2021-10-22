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
try:
    from BytesIO import BytesIO
except ImportError:
    from io import BytesIO
from stdatmos import std_atmosphere_pres
from utils.cron_helper import *
from utils.configs import *

# No username and passwords here!! Create a .env file containing USERNAME and PASSWORD
# variables
from dotenv import load_dotenv
load_dotenv()
username = os.environ.get("USERNAME")
password = os.environ.get("PASSWORD")

def get_times(data_path, feed='acars'):
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
    feed : string [Optional]
        Which MADIS feed to download. Options are: acars, metar. Defaults to acars

    Returns
    -------
    times_to_dl : list
        List of Python datetime objects indicating which times to download from MADIS
    """

    def touch(fname, times=None):
        with open(fname, 'a'):
            os.utime(fname, times)

    marker_path = "%s/markers" % data_path
    try:
        os.mkdir(marker_path)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass

    # Figure out the files and their update times from the MADIS server
    URL = "%s/%s/netcdf" % (base_url, feed)
    req = urlreq.Request(URL)
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
            touch(marker_fname)
            times_to_dl.append(file_time)

    # Check for old marker files and delete them if they're older than the oldest file
    # on the MADIS server
    earliest_time = datetime.strptime(files[0], '%Y%m%d_%H%M')
    for fname in glob.glob("%s/*.txt" % marker_path):
        ftime = datetime.strptime(os.path.basename(fname), "%Y%m%d_%H%M.txt")
        if ftime < earliest_time:
            os.unlink(fname)

    return times_to_dl

def get_madis(data_path, dt, feed='acars', realtime=True):
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
        dl_url = "%s/%s/netcdf" % (base_url, feed)
    else:
        dl_url = "%s/%s/%s/%s/point/%s/netcdf" % (archive_url,
                                                  str(dt.year).zfill(4),
                                                  str(dt.month).zfill(2),
                                                  str(dt.day).zfill(2),
                                                  feed)
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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-a', '--archive-mode', dest="archive_mode", help="YYYY-MM-DD/HH")
    ap.add_argument('-f', '--feed', dest="feed", help="MADIS feed to retrieve. [acars, metar]")
    ap.add_argument('-n', '--num-hours', dest="num_hours", help="If archive-mode set, number of hours to search backwards")
    ap.add_argument('-p', '--output-path', dest='output_path', help="Local path to store data")
    args = ap.parse_args()

    if args.output_path is None: args.output_path = NETCDF
    if args.feed is None: args.feed = 'acars'
    args.output_path = "%s/%s" % (args.output_path, args.feed)
    dts = get_times(args.output_path, feed=args.feed)

    # ----------------------------------------------------------------------------------
    # For archive "lookback" functionality...
    # ----------------------------------------------------------------------------------
    if args.archive_mode is None:
        current_time = datetime.utcnow()
        realtime = True
    else:
        current_time = datetime.strptime(args.archive_mode, "%Y-%m-%d/%H")
        if args.num_hours is None: args.num_hours = 6
        dts = []
        for t in range(int(args.num_hours)):
            dts.append(current_time - timedelta(hours=t))
        realtime = False

    timestamp("INFO")
    print("Initializing new download request --------------------------")

    dt_fnames = {}
    for dt in dts:
        print(">>>> ", dt)
        get_madis(args.output_path, dt, feed=args.feed, realtime=realtime)

if __name__ == "__main__":
    main()
