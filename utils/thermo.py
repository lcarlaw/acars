import numpy as np
from numpy.core.numeric import normalize_axis_index
from .constants import *

sat_pressure_0c = 6.112
c1 = 0.0498646455
c2 = 2.4082965
c3 = 7.07475
c4 = 38.9114
c5 = 0.0915
c6 = 1.2035

def mixing_ratio_from_specific_humidity(specific_humidity):
    """Calculate the mixing ratio from specific humidity.

    Parameters
    ----------
    specific_humidity: `pint.Quantity`
        Specific humidity of air

    Returns
    -------
    `pint.Quantity`
        Mixing ratio

    Notes
    -----
    Formula from [Salby1996]_ pg. 118.

    .. math:: w = \frac{q}{1-q}

    * :math:`w` is mixing ratio
    * :math:`q` is the specific humidity

    See Also
    --------
    mixing_ratio, specific_humidity_from_mixing_ratio
    """
    return specific_humidity / (1 - specific_humidity)

def saturation_mixing_ratio(tot_press, temperature):
    """Calculate the saturation mixing ratio of water vapor.

    This calculation is given total pressure and the temperature. The
    implementation uses the formula outlined in [Hobbs1977]_ pg.73.

    Parameters
    ----------
    tot_press: `pint.Quantity`
        Total atmospheric pressure
    temperature: `pint.Quantity`
        air temperature

    Returns
    -------
    `pint.Quantity`
        The saturation mixing ratio, dimensionless

    """
    return mixing_ratio(saturation_vapor_pressure(temperature), tot_press)


def mixing_ratio(part_press, tot_press):
    """Calculate the mixing ratio of a gas.

    This calculates mixing ratio given its partial pressure and the total pressure of
    the air. There are no required units for the input arrays, other than that they have
    the same units.

    Parameters
    ----------
    part_press : `pint.Quantity`
        Partial pressure of the constituent gas
    tot_press : `pint.Quantity`
        Total air pressure
    molecular_weight_ratio : `pint.Quantity` or float, optional
        The ratio of the molecular weight of the constituent gas to that assumed for
        air. Defaults to the ratio for water vapor to dry air
        (:math:`\epsilon\approx0.622`).

    Returns
    -------
    `pint.Quantity`
        The (mass) mixing ratio, dimensionless (e.g. Kg/Kg or g/g)

    Notes
    -----
    This function is a straightforward implementation of the equation given in many
    places, such as [Hobbs1977]_ pg.73:

    .. math:: r = \epsilon \frac{e}{p - e}

    See Also
    --------
    saturation_mixing_ratio, vapor_pressure
    """
    return eps * part_press / (tot_press - part_press)

def saturation_vapor_pressure(temperature):
    """Calculate the saturation water vapor (partial) pressure."""
    # Converted from original in terms of C to use kelvin. Using raw absolute
    # values of C in a formula plays havoc with units support.
    return sat_pressure_0c * np.exp(
        17.67 * (temperature - 273.15) / (temperature - 29.65)
    )

def vapor_pressure(pressure, mixing):
    """Calculate water vapor (partial) pressure.
    Given total `pressure` and water vapor `mixing` ratio, calculates the partial
    pressure of water vapor."""
    return pressure * mixing / (mpconsts.epsilon + mixing)

def theta_e(pressure, temperature, dewpoint):
    """Calculate equivalent potential temperature.

    This calculation must be given an air parcel's pressure, temperature, and dewpoint.
    The implementation uses the formula outlined in [Bolton1980]_:

    First, the LCL temperature is calculated:

    .. math:: T_{L}=\frac{1}{\frac{1}{T_{D}-56}+\frac{ln(T_{K}/T_{D})}{800}}+56

    Which is then used to calculate the potential temperature at the LCL:

    .. math:: \theta_{DL}=T_{K}\left(\frac{1000}{p-e}\right)^k
              \left(\frac{T_{K}}{T_{L}}\right)^{.28r}

    Both of these are used to calculate the final equivalent potential temperature:

    .. math:: \theta_{E}=\theta_{DL}\exp\left[\left(\frac{3036.}{T_{L}}
                                              -1.78\right)*r(1+.448r)\right]

    Parameters
    ----------
    pressure: numpy array
        Total atmospheric pressure [hPa]
    temperature: numpy array
        Temperature of parcel [C]
    dewpoint: numpy array
        Dewpoint of parcel [C]

    Returns
    -------
    th_e : numpy array
        The equivalent potential temperature of the parcel [Kelvin]

    Notes
    -----
    [Bolton1980]_ formula for Theta-e is used, since according to
    [DaviesJones2009]_ it is the most accurate non-iterative formulation
    available.

    """
    t = temperature + ZEROCNK
    td = dewpoint + ZEROCNK
    p = pressure
    e = saturation_vapor_pressure(td)
    r = saturation_mixing_ratio(p, td)
    t_l = 56 + 1. / (1. / (td - 56) + np.log(t / td) / 800.)
    th_l = t * (1000 / (p - e)) ** ROCP * (t / t_l) ** (0.28 * r)
    th_e = th_l * np.exp((3036. / t_l - 1.78) * r * (1 + 0.448 * r))
    return th_e

def wet_bulb_temperature(pressure, temperature, dewpoint):
    """Calculate the wet-bulb temperature using Normand's rule.

    This function calculates the wet-bulb temperature using the Normand method. The LCL
    is computed, and that parcel brought down to the starting pressure along a moist
    adiabat. The Normand method (and others) are described and compared by [Knox2017]_.

    Parameters
    ----------
    pressure : `pint.Quantity`
        Initial atmospheric pressure
    temperature : `pint.Quantity`
        Initial atmospheric temperature
    dewpoint : `pint.Quantity`
        Initial atmospheric dewpoint

    Returns
    -------
    `pint.Quantity`
        Wet-bulb temperature

    See Also
    --------
    lcl, moist_lapse
    """
    it = np.nditer([pressure, temperature, dewpoint, None],
                   op_dtypes=['float', 'float', 'float', 'float'],
                   flags=['buffered'])

    for press, temp, dewp, ret in it:
        press = press
        temp = temp
        dewp = dewp
        lcl_pressure, lcl_temperature = lcl(press, temp, dewp)
        moist_adiabat_temperatures = moist_lapse(concatenate([lcl_pressure, press]),
                                                 lcl_temperature)
        ret[...] = moist_adiabat_temperatures[-1].magnitude

    # If we started with a scalar, return a scalar
    if it.operands[3].size == 1:
        return it.operands[3][0] * moist_adiabat_temperatures.units
    return it.operands[3] * moist_adiabat_temperatures.units
