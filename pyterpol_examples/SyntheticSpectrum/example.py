"""
This is a tutorial script how to handle the class Synthetic Spectrum.
"""
import pyterpol
import numpy as np
import matplotlib.pyplot as plt

# Load the spectrum using the library numpy
wave, intens = np.loadtxt('grid.dat', unpack=True, usecols=[0,1])

# The synthetic spectrum can be created either from arrays
ss = pyterpol.SyntheticSpectrum(wave=wave, intens=intens)

# or just loaded from an array
ss = pyterpol.SyntheticSpectrum(f='grid.dat')

# usually your spectra have additional properties.
# Note that parameters passed during the contruction
# of teh class do not have impact on the underlying spectrum
ss = pyterpol.SyntheticSpectrum(f='grid.dat', teff=18000, logg=4.5, idiotic_parameter='pes')

# which you can then easily review
print ss

# or change
ss['teff'] = 12000

# and review the changes
print ss

# Or we can directly view the spectrum
ss.plot(savefig=True, figname='original.png')

# the spectrum can be rotated
newwave, newintens = ss.get_spectrum(vrot=50)

# shifted in rv
newwave, newintens = ss.get_spectrum(rv=50)

# shrinked
newwave, newintens = ss.get_spectrum(lr=0.7)

# or transformed to KOREL
newwave, newintens = ss.get_spectrum(korel=True)

# or all together
newwave, newintens = ss.get_spectrum(vrot=50, lr=0.7, rv=50, korel=True)

# lets wrap this into SyntheticSpectrum and plot it
nss = pyterpol.SyntheticSpectrum(wave=newwave, intens=newintens)
nss.plot(savefig=True, figname='adjusted.png')
plt.show()

