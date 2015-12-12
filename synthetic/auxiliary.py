import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.interpolate import splrep
from scipy.interpolate import splev
from scipy.interpolate import bisplrep
from scipy.interpolate import bisplev
from scipy.interpolate import RectBivariateSpline
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


def interpolate_block(x, block, xnew):
    """
    Interpolates in each column of a 2d matrix.
    input:
	x.. abscissa
	block.. 2d matrix of 
    output:
	intens.. interpolated
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
        # f = interp1d(np.array(x), np.array(y), kind='cubic')
        # intens[i] = f(xnew)
        # intens[i] = spline(x, y, xnew, order=k)

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

    if ny > 3:
        ky = 3
    else:
        ky = ny - 1

    f = RectBivariateSpline(x, np.arange(nx), block, kx=ky, ky=2)
    intens = f(xnew, np.arange(nx))[0]

    return intens

def interpolate_spec(wave0, intens0, wave1):
    """
    input:
      wave0, intens0.. spectrum to be interpolated
      wave1.. new set of wavelengths
    output:
      intens1.. new set of intensities
    """
    tck = splrep(wave0, intens0)
    intens1 = splev(wave1, tck)

    return intens1


def is_within_interval(v, arr):
    """
    Constrols that value v
    lies within a specified array.
    input:
      v.. tested value
      arr.. tested array
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
    Reads an ASCII file and returns its string rep.
    Input:
      f.. the file to be read
    Output:
      lines.. string representation of the file
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

def rotate_spectrum(wave, intens, vrot, fwhm=0.0, epsilon=0.6):
    """
    input:
      wave.. input wavelength
      intens.. input intensities
      vrot.. vsini in km/s
      fwhm.. instrumental broadening profile
    output
      intens1.. rotated intensities
    """

    if vrot > ZERO_TOLERANCE:
        # we need it equidistant in RV
        # lambda_n = lambda_0*(1+RV/c)**n
        # n = len(wave)
        # RV = ((wave[-1]/wave[0])**(1/(n-1))-1)/c.value
        # indices = np.arange(n)
        # wave_rv = wave[0]*(1-RV/c.value)**indices
        wave_log = np.log(wave)
        rv = np.linspace(wave_log[0], wave_log[-1], len(wave))
        step = rv[1] - rv[0]

        # interpolate
        intens_rv = interpolate_spec(wave_log, intens, rv)

        # scale rotational velocity with light speed
        vrot = 1000 * vrot / c.value

        # get the kernel
        # velocity vector
        n = int(2 * vrot / step)
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
        intens = interpolate_spec(wave_conv, 1 - intens_conv, wave)

    return intens


def shift_spectrum(wave, RV):
    """
    input:
      wave.. old wavelength array
      RV.. radial velocity in km/s
    output:
      new_wave.. shifted array of wavelengths
    """
    # shifts the wavelengths
    new_wave = wave * (1 + RV * 1000 / c.value)

    return new_wave

def string2bool(s):
    """
    Converts string to boolean
    :param s:
    :return:
    """
    if s.lower() in ['true', '1']:
        return True
    else:
        return False

def write_numpy(f, cols, fmt):
    """
    :param f: outputfile or handler
    :param cols: block of data to be writte
    :param fmt: format of the blocs
    :return:
    """

    np.savetxt(f, cols, fmt=fmt)
