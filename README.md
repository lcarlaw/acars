
# acars
Python code for downloading ACARS data and converting to BUFKIT- and SHARPpy-readable files. In addition, this code attempts to augment the at-times-questionable ACARS data with short-term forecast data from the hourly-cycled Rapid Refresh model.  

## Basic Setup Notes
Out of the box, this code will not work since we have to access a user-restricted MADIS feed to download the ACARS/AMDAR data. A username/password pair are stored in a secret `.env` file and read in at execution time. See the [ESRL FAQ here](https://amdar.noaa.gov/FAQ.html) for more information and application.

### Creating the base environment
The setup here proceeds using Anaconda, as well as assuming a completely vanilla Python3 install.  This code will not work with Python2--you should probably upgrade to version 3 at this point.  I've also edited my `~/.condarc` file to add conda-forge to the default channels.

```
conda create --name acars python=3.7
conda activate acars

conda install netcdf4 numpy pandas requests
pip install python-dotenv

[OPTIONAL for crons and uploading to Google Drive]:
conda install pydrive
```
#### Google Drive API
This is optional, but will allow for automated uploads to Google Drive. Follow the [initialization steps outlined here](https://mihevc.org/2016/02/04/crontabed-pydrive-uploader.html). This is currently going to the `lcarlaw.data.store@gmail.com` account.

## Usage
Two scripts are employed for general use: one for downloading and the MADIS netCDF files, and another for constructing the BUFKIT-readable soundings.

You'll need to edit a few variables in the `utils/configs.py` script. These are indicated clearly within the comments, but include things such as setting the `${DATA}`, `${SOUNDINGS}`, and `${NETCDF}` paths.

### Data Download
All arguments are optional:

```
python download_acars.py [-a ARCHIVE_MODE] [-v]
```
`ARCHIVE_MODE`:  Specifying a date time of the form `YYYY-MM-DD/HH` which cause the script to download archived ACARS data (the archive goes back to about 2001). Archived output will be dumped into `${DATA}/data_store_archive`
`-v`: Optional switch to download only VAPOR profiles.

### BUFKIT/SHARPpy Output

```
python acars.py [-a]
```
`-a`: Set this flag if this is an archived run. Script will then search `${DATA}/data_store_archive` instead.
