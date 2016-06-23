#!/usr/bin/env python

"""
ready_phoenix.py
Convert Phoenix synthetic spectra from FITS to DAT.

"""

__author__ = "Miroslav Broz (miroslav.broz@email.cz)"
__version__ = "Jun 23rd 2016"

import sys
import numpy as np

from scipy.interpolate import splrep, splev
from astropy.io import fits
from pyterpol.synthetic.auxiliary import instrumental_broadening

def read_wave(filename):
    """Read wavelength data"""
    hdu = fits.open(filename)
    wave = hdu[0].data

    hdu.info()
    print("")
    print(repr(hdu[0].header))
    print("")
    print(wave)
    print("")
    hdu.close()
    return wave

def fits2dat(filename, wave):
    """Convert Phoenix synthetic spectrum from FITS to DAT."""
    hdu = fits.open(filename)
    intens = hdu[0].data

    hdu.info()
    print("")
    print(repr(hdu[0].header))
    print("")
    print(intens)
    print("")
    hdu.close()
    
#    np.savetxt(filename[:-4]+'dat.0', np.column_stack([wave, intens]), fmt="%.6e %.10e")  # dbg

    # convolution (to 1 Angstrom)
    step = 1.0
    intens = instrumental_broadening(wave, intens, width=step)
    s = np.array(zip(wave, intens))

    print(intens)
    print("")
    print(s)
    print("")

    # spline interpolation
    wmin = s[0,0]
    wmax = s[-1,0]
    wnew = np.arange(wmin, wmax, step)
    
    tck = splrep(s[:,0], s[:,1])
    s_new = splev(wnew, tck)
    intens = s_new

    # elliminate negatives!
    for i in xrange(0,len(intens)):
        if intens[i] < 0.0:
            intens[i] = 0.0
        intens[i] = intens[i]*1.e-8  # erg s^-1 cm^-2 cm^-1 -> erg s^-1 cm^-2 A^-1 (as in POLLUX)

    # save spectra
    out = filename[:-4]+'vis.dat'
    np.savetxt(out, np.column_stack([wnew, intens]), fmt="%.6e %.10e")

    sys.exit(1)  # dbg

def main():
    """Convert all files"""
    if len(sys.argv) > 1:
        inlist = sys.argv[1:]
    else:
        inlist = np.loadtxt("inlist", dtype=str)

    wave = read_wave("WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")

    for filename in inlist:
        fits2dat(filename, wave)

if __name__ == "__main__":
    main()


