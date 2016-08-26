import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.interpolate import splrep
from scipy.interpolate import splev
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import spline
from scipy.signal import fftconvolve

ZERO_TOLERANCE = 1e-6


def flatten_2d(arr):
    """
    Flattens 2-dim array

    :param arr: 2d array
    :return:
    """
    newarr = []
    if any([isinstance(subarr, (list, tuple)) for subarr in arr]):
        for subarr in arr:
            if isinstance(subarr, (tuple, list)):
                newarr.extend(subarr)
            else:
                newarr.append(subarr)
        return newarr
    else:
        return arr


def instrumental_broadening(wave, flux, width=0.25, width_type='sigma', interpolate_back=True):
    """
    A convolution of a spectrum with a normal distribution.

    :param: wave:
    :param: flux:
    :param width:
    :param width_type:
    :return:
    """
    # print "Computing instr. broadening."
    # If there is no broadening to apply, don't bother
    if width < ZERO_TOLERANCE:
        return flux

    # Convert user input width type to sigma (standard devation)
    width_type = width_type.lower()
    if width_type == 'fwhm':
        sigma = width / 2.3548
    elif width_type == 'sigma':
        sigma = width
    else:
        raise ValueError(("Unrecognised width_type='{}' (must be one of 'fwhm'"
                          "or 'sigma')").format(width_type))

    # Make sure the wavelength range is equidistant before applying the
    # convolution
    delta_wave = np.diff(wave).min()
    range_wave = wave.ptp()
    n_wave = int(range_wave / delta_wave) + 1
    wave_ = np.linspace(wave[0], wave[-1], n_wave)
    # flux_ = np.interp(wave_, wave, flux)
    flux_ = interpolate_spec(wave, flux, wave_)
    dwave = wave_[1] - wave_[0]
    n_kernel = int(2 * 4 * sigma / dwave)

    # The kernel might be of too low resolution, or the the wavelength range
    # might be too narrow. In both cases, raise an appropriate error
    if n_kernel == 0:
        raise ValueError(("Spectrum resolution too low for "
                          "instrumental broadening (delta_wave={}, "
                          "width={}").format(delta_wave, width))
    elif n_kernel > n_wave:
        raise ValueError(("Spectrum range too narrow for "
                          "instrumental broadening"))

    # Construct the broadening kernel
    wave_k = np.arange(n_kernel) * dwave
    wave_k -= wave_k[-1] / 2.
    kernel = np.exp(- (wave_k) ** 2 / (2 * sigma ** 2))
    kernel /= sum(kernel)

    # Convolve the flux with the kernel
    flux_conv = fftconvolve(1 - flux_, kernel, mode='same')

    # And interpolate the results back on to the original wavelength array,
    # taking care of even vs. odd-length kernels
    if n_kernel % 2 == 1:
        offset = 0.0
    else:
        offset = dwave / 2.0

    if interpolate_back:
        flux = np.interp(wave + offset, wave_, 1 - flux_conv, left=1, right=1)
        # flux = interpolate_spec(wave_, 1-flux_conv, wave+offset)
        # Return the results.
        return flux


def interpolate_block(x, block, xnew):
    """
    Interpolates in each line of a 2d array.

    :param x: independent variable
    :type x: numpy.float64
    :param block: 2d array for each column f(x)= block[i]
    :type block: numpy.float64
    :param xnew: point at which it is interpolated
    :type xnew: float
    :return:
    """
    intens = np.zeros(len(block[0]))
    n = len(block[:, 0])

    # set up the order of interpolation
    if n > 4:
        k = 3
    else:
        k = n - 1
    # k=3

    # TODO Can thius be done faster with bisplrep and bisplev
    # do the interpolation
    for i in range(0, len(block[0])):
        y = block[:, i]

        tck = splrep(x, y, k=k)
        intens[i] = splev(xnew, tck, der=0)

    return intens


def interpolate_block_faster(x, block, xnew):
    """
    Interpolation of teh spectra... hopefully faster?

    :param x:
    :param block:
    :param xnew:
    :return:
    """

    # length of the datablock
    nx = len(block[0])
    ny = len(x)
    # print x

    if (ny > 3) & (ny < 6):
        ky = 3
    elif ny > 5:
        ky = 5
    else:
        ky = ny - 1

    # print ky

    f = RectBivariateSpline(x, np.arange(nx), block, kx=ky, ky=1)
    intens = f(xnew, np.arange(nx))[0]

    return intens


def interpolate_spec(wave0, intens0, wave1):
    """
    Defines a function intens0 = f(wave0) and
    than interpolates in it at wave1.

    :param wave0: initial wavelength array
    :type wave0: numpy.float64
    :param intens0: initial intensity array
    :type intens0: numpy.float64
    :param wave1: wavelength array at which we interpolate
    :type wave1: numpy.float64
    :return intens1: final intensity array
    :rtype intens1: numpy.float64
    """
    tck = splrep(wave0, intens0, k=3)
    intens1 = splev(wave1, tck)

    return intens1


def is_within_interval(v, arr):
    """
    Tests whether value v lies within interval [min(arr); max(arr)]

    :param v: tested values
    :type v: numpy.float64
    :param arr: tested array
    :type v: numpy.float64
    :return:
    :param:
    :type: bool
    """
    # print v, max(arr), min(arr)
    if (v - max(arr) > ZERO_TOLERANCE) | (min(arr) - v > ZERO_TOLERANCE):
        return False
    else:
        return True


def generate_least_number(l):
    """
    Goes over integer in list and finds the
    smallest integer not in the list.

    :param l: the list
    :return: int the smallest integer
    """
    num = 0
    while num in l:
        num += 1
    return num


def keys_to_lowercase(d):
    """
    Converts dictionary keys to lowercase

    :param d the converted dictionary
    :return: dnew
    """

    dnew = {}
    for key in d.keys():
        keynew = key.lower()
        dnew[keynew] = d[key]

    return dnew


def parlist_to_list(l, property='value'):
    """
    Converts a list of Parameter class to a
    regular list - only the property is returned

    :param l:
    :param prop:
    :return:
    """
    ol = []
    for par in l:
        ol.append(par[property])

    return ol


def sum_dict_keys(d):
    """
    Sums dictionary key records.

    :param d: the dictionary
    :return: s the sum
    """
    s = 0.0
    for key in d.keys():
        s += d[key]
    return s


def read_text_file(f):
    """
    Reads ascii file f.

    :param f: the file
    :type f: str
    :return lines: list of all lines within file f
    :rtype: list
    """

    ifile = open(f, 'r')
    lines = ifile.readlines()
    ifile.close()

    return lines


def renew_file(f):
    """
    Deletes an existing file.

    :param f:
    :return:
    """
    ofile = open(f, 'w')
    ofile.close()


def rotate_spectrum(wave, intens, vrot, epsilon=0.6, interpolate_back=True):
    """
    Rotates a spectrum represented by arrays wave and intes to the prjected
    rotational velocity vrot.

    :param wave: wavelength array
    :type wave: numpy.float64
    :param intens: intensity array
    :type intens: numpy.float64
    :param vrot: projected rotational velocity in km/s
    :type vrot: float
    :param epsilon: Coefficient of linear limb-darkening.
    :type epsilon: float
    :param interpolate_back: interpolate the spectrum back to the original wavelength sampling
    :type interpolate_back: bool
    :return intens: the rotated spectrum in the original wavelength sanmpling
    :rtype intens: numpy.float64
    :return intens_conv: the rotated spectrum equidistant in rv
    :rtype intens_conv: numpy.float64
    :return wave_conv: the wavelength array equidistant in rv
    :rtype wave_conv: numpy.float64
    """
    if vrot > ZERO_TOLERANCE:
        # we need it equidistant in RV
        wave_log = np.log(wave)
        rv = np.linspace(wave_log[0], wave_log[-1], len(wave))
        step = rv[1] - rv[0]

        # interpolate
        intens_rv = interpolate_spec(wave_log, intens, rv)

        # scale rotational velocity with light speed
        vrot = 1000 * vrot / c.value

        # get the kernel
        # velocity vector
        n = int(np.ceil(2 * vrot / step))
        rv_ker = np.arange(n) * step
        rv_ker = rv_ker - rv_ker[-1] / 2.
        y = 1 - (rv_ker / vrot) ** 2

        # the kernel
        kernel = (2 * (1 - epsilon) * np.sqrt(y) + np.pi * epsilon / 2. * y) / (np.pi * vrot * (1 - epsilon / 3.0))
        kernel = kernel / kernel.sum()

        # convolve the flux
        intens_conv = fftconvolve(1 - intens_rv, kernel, mode='same')
        if n % 2 == 1:
            rv = np.arange(len(intens_conv)) * step + rv[0]
        else:
            rv = np.arange(len(intens_conv)) * step + rv[0] - step / 2.

        wave_conv = np.exp(rv)

        # interpolate back
        if interpolate_back:
            intens = interpolate_spec(wave_conv, 1 - intens_conv, wave)
            return intens
        else:
            return 1 - intens_conv, wave_conv


def shift_spectrum(wave, RV):
    """
    Doppler-shifts spectrum.
    :param wave: original wavelength array
    :type wave: numpy.float64
    :param RV: radial velocity in km/s
    :type RV: float
    :return new_wave: shifted wavelength array
    :rtype new_wave: numpy.float64

    """
    # shifts the wavelengths
    new_wave = wave * (1 + RV * 1000 / c.value)

    return new_wave


def select_index_for_multiple_keywords(d, **kwargs):
    """
    From a dictionary of lists selects
    one index meeting all requirements.

    :param kwargs:
    :return:
    """
    keys = d.keys()
    length = len(d[keys[0]])

    for i in range(0, length):
        for k in keys:
            if d[k] == kwargs[k] and k == keys[-1]:
                return i
    return -1


def string2bool(s):
    """
    Converts string to boolean.

    :param s:
    :return:
    """
    if s.lower() in ['true', '1']:
        return True
    else:
        return False


def write_numpy(f, cols, fmt):
    """
    An example of lack of brain of the main developer of this "code".

    :param f: outputfile or handler
    :param cols: block of data to be writte
    :param fmt: format of the blocs
    :return: None
    """

    np.savetxt(f, cols, fmt=fmt)
