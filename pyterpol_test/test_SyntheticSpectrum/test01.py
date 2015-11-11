import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
from pyterpol.synthetic.makespectrum import SyntheticSpectrum

syspe = SyntheticSpectrum(f='temp.dat', teff=10000, logg=4.2)

# check normal output
print syspe.get_spectrum()

# choose a wave region
# this one should raise warning
wave = np.linspace(6500, 6600, 800) 
base_wave, base_intens = syspe.get_spectrum(wave)

# get a shifted spectrum
wave = np.linspace(6500, 6600, 800)
rv_wave, rv_intens =  syspe.get_spectrum(wave, rv=50)
np.savetxt('shifted.dat', np.column_stack([rv_wave, rv_intens]), fmt="%20.10e")

# get rotated spectrum
wave = np.linspace(6500, 6600, 800)
rot_wave, rot_intens =  syspe.get_spectrum(wave, vrot=10)
print rot_wave, rot_intens
np.savetxt('rotated.dat', np.column_stack([rv_wave, rv_intens]), fmt="%20.10e")


plt.plot(base_wave, base_intens, 'k-')
plt.plot(rv_wave, rv_intens, 'r-')
plt.show()