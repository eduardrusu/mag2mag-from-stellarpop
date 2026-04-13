"""
A module of tools to work with SED, spectral, and filter operations.
"""

import glob
import os
from math import log10, pi

from numpy import loadtxt
from scipy.interpolate import splev, splint, splrep

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
filterpath = os.path.join(BASE_DIR, "filters/")
SEDpath = os.path.join(BASE_DIR, "templates/")

filter_list = glob.glob(filterpath + "*res")
filter_list.sort()
filter_list = [filter.split("/")[-1].split(".res")[0] for filter in filter_list]
sed_list = glob.glob(SEDpath + "*sed")
sed_list.sort()
sed_list = [sed.split("/")[-1].split(".sed")[0] for sed in sed_list]

BC03factor = 3.826e33 / (4 * pi * 3.08568e19**2)


def filterfromfile(file):
    """
    Create a filter model from a file.
    """
    f = open(filterpath + file + ".res")
    filter = loadtxt(f)
    f.close()
    return splrep(filter[:, 0], filter[:, 1], k=1, s=0)


def get_sed(name):
    """
    Returns a model of the SED, a tuple of (wave,data)
    """
    f = open(SEDpath + name + ".sed")
    sed = loadtxt(f)
    f.close()
    return sed[:, 0], sed[:, 1]


def ab_filter_magnitude(filter, spectrum, redshift):
    """
    Determines the AB magnitude (up to a constant) given an input filter, SED,
        and redshift. Does not include cosmological dimming, so should only
        be used for relative magnitudes.
    """
    sol = 299792452.0  # speed of light in vacuum in m/s

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Convert to f_nu
    data = data * wave**2 / (sol * 1e10)

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= 1.0 + redshift
    wmin, wmax = filter[0][0], filter[0][-1]
    cond = (wave >= wmin) & (wave <= wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond], filter)

    freq = sol * 1e10 / wave[cond]
    data = data[cond] * (1.0 + redshift)

    # Flip arrays
    freq = freq[::-1]
    data = data[::-1]
    response = response[::-1]

    # Integrate
    observed = splrep(freq, response * data / freq, s=0, k=1)
    flux = splint(freq[0], freq[-1], observed)

    bp = splrep(freq, response / freq, s=0, k=1)
    bandpass = splint(freq[0], freq[-1], bp)

    return -2.5 * log10(flux / bandpass) - 48.6


def vega_filter_magnitude(filter, spectrum, redshift):
    """
    Determines the Vega magnitude (up to a constant) given an input filter,
        SED, and redshift.
    """

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= 1.0 + redshift
    data /= 1.0 + redshift
    wmin, wmax = filter[0][0], filter[0][-1]
    cond = (wave >= wmin) & (wave <= wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond], filter)

    # Determine the total observed flux (without the bandpass correction)
    observed = splrep(wave[cond], (response * data[cond]), s=0, k=1)
    flux = splint(wmin, wmax, observed)

    # Determine the magnitude of Vega through the filter
    vwave, vdata = get_sed("Vega")
    cond = (vwave >= wmin) & (vwave <= wmax)
    response = splev(vwave[cond], filter)
    vega = splrep(vwave[cond], response * vdata[cond], s=0, k=1)
    vegacorr = splint(wmin, wmax, vega)

    return -2.5 * log10(flux / vegacorr)  # +2.5*log10(1.+redshift)


def filter_magnitude(filter, spectrum, redshift, zp):

    wave = spectrum[0].copy()
    data = spectrum[1].copy()

    # Redshift the spectrum and determine the valid range of wavelengths
    wave *= 1.0 + redshift
    data /= 1.0 + redshift
    wmin, wmax = filter[0][0], filter[0][-1]
    cond = (wave >= wmin) & (wave <= wmax)

    # Evaluate the filter at the wavelengths of the spectrum
    response = splev(wave[cond], filter)

    # Determine the total observed flux (without the bandpass correction)
    observed = splrep(wave[cond], (response * data[cond]), s=0, k=1)
    flux = splint(wmin, wmax, observed)
    flux = (response * data[cond]).sum()
    return -2.5 * log10(flux) + zp
