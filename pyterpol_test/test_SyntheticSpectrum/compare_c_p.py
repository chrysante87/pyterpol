import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
from pyterpol.synthetic.makespectrum import SyntheticSpectrum

gridSpec = 'grid.dat'

# basic gridspectrum used by interpol 
c0_wave, c0_int = np.loadtxt(gridSpec, unpack=True, usecols=[0,1])

# shifted and rotated gridspectrum
crs_wave, crs_int = np.loadtxt('output012', unpack=True, usecols=[0,1])

syspe = SyntheticSpectrum(f=gridSpec)
py0_wave, py0_int = syspe.get_spectrum()
pyrs_wave, pyrs_int = syspe.get_spectrum(rv=30, vrot=50)

# get padded spectrum
wave, intens = syspe.get_spectrum()
p_wave, p_intens = syspe.pad_continuum(wave, intens, 10)


plt.subplot(211)
plt.plot(c0_wave, c0_int, 'k-', lw=1)
plt.plot(crs_wave, crs_int, 'k-', lw=2)
plt.plot(py0_wave, py0_int+0.5, 'r-', lw=1)
plt.plot(pyrs_wave, pyrs_int+0.5, 'r-', lw=2)
plt.plot(p_wave, p_intens+0.5, 'b-')


plt.subplot(212)
plt.plot(c0_wave, c0_int-py0_int, 'y-', lw=1)
plt.plot(c0_wave, crs_int-pyrs_int+0.1, 'y-', lw=2)
plt.show()