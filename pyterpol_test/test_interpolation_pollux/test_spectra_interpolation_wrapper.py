import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol
import time

# create the grid - custom does nothing
sygri = pyterpol.SyntheticGrid('POLLUX', debug=True)

# interpolate at
wmin = 6500
wmax = 6650
keys = ['Z', 'LOGG', 'TEFF']
values = [1.0, 4.0, 12000]
params = {k:v for k,v in zip(keys, values)}
name = '_'.join([k+'_'+str(v) for k,v in zip(keys, values)])

plt.figure(figsize=(12,5), dpi=100)
spectrum = sygri.get_synthetic_spectrum(params, wmin, wmax)
w,i = spectrum.get_spectrum()
plt.plot(w, i, 'k-')
plt.savefig(name+'.png')
plt.close()



