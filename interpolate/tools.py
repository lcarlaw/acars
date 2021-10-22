"""Tools to assist with various interpolation routines"""

from scipy import spatial
from scipy.spatial.distance import cdist
import numpy as np

def calc_distance(points, lon, lat):
    """Search for the nearest grid point using a KDTree

    Parameters
    ----------
    points : list
        List of lists indicating lon/lat pais: [[LON1, LAT1], [LON2, LAT2]]
    lat : np.array
        2-d array of gridded latitudes
    lon : np.arrat
        2-d array of gridded longitudes
    """

    lonlat = np.column_stack((lon.ravel(), lat.ravel()))
    tree = spatial.cKDTree(lonlat)
    dist, idx = tree.query(points, k=1)
    #print(tree.query_ball_point(points, 1, p=np.inf))
    return dist

def interp_pres(p, pres, field):
    """
    Generic interpolation routine. Converts pressures into log10 coordinates

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

def calc_kappa(spacing, kappa_star=5.052):
    r"""Calculate the kappa parameter for barnes interpolation.

    Parameters
    ----------
    spacing: float
        Average spacing between observations
    kappa_star: float
        Non-dimensional response parameter. Default 5.052.

    Returns
    -------
        kappa: float

    """
    return kappa_star * (2.0 * spacing / np.pi)**2

def average_spacing(points):
    """Calculate the average spacing to the nearest other point.

    Parameters
    ----------
    points : (M, N) array_like
         M points in N dimensional space

    Returns
    -------
        The average distance to the nearest neighbor across all points

    """
    dist_matrix = cdist(points, points)
    diag = np.arange(len(dist_matrix))
    dist_matrix[diag, diag] = np.nan
    return np.nanmin(dist_matrix, axis=0).mean()
