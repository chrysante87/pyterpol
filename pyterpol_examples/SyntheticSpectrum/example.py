"""
This is a tutorial script how to handle the class Synthetic Spectrum.
"""
import pyterpol
import numpy as np
import matplotlib.pyplot as plt

# load the spectrum using numpuy
wave, intens = np.loadtxt('grid.dat', unpack=True, usecols=[0,1])

# The synthetic spectrum can be created either from arrays
ss = pyterpol.SyntheticSpectrum(wave=wave, intens=intens)

# or just loaded from an array
ss = pyterpol.SyntheticSpectrum(f='grid.dat')

# usually your spectra have additional properties
ss = pyterpol.SyntheticSpectrum(f='grid.dat', teff=18000, logg=4.5)

# which you can then easily review
print ss

# or change
ss['teff'] = 12000
print ss

# Or we can directly view the spectrum
ss.plot(savefig=True, figname='original.png')

# In order to rotate the spectrum or shift in rv
# one uses the method get_spectrum
newwave, newintens = ss.get_spectrum(rv=100, vrot=50)

# lets wrap this into SyntheticSpectrum and plot it
nss = pyterpol.SyntheticSpectrum(wave=newwave, intens=newintens)
nss.plot(savefig=True, figname='adjusted.png')
plt.show()

