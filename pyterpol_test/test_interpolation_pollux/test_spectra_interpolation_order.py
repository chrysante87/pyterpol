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
values = [1.0, 3.8, 12200]
params = {k:v for k,v in zip(keys, values)}
name = '_'.join([k+'_'+str(v) for k,v in zip(keys, values)])

ords = np.arange(2,6)

for order in ords:
    spectrum = sygri.get_synthetic_spectrum(params, wmin, wmax, order=order)
    w,i = spectrum.get_spectrum()
    
    #plot the spectrum 
    plt.figure(figsize=(12,5), dpi=100)
    plt.plot(w, i, 'k-')
    plt.savefig(name+"_order_%s" % str(order)+'.png')
    plt.close()
    
    # save the spectrum in text file
    np.savetxt(name+"_order_%s" % str(order)+'.dat', np.column_stack([w,i]), fmt="%.4f %.6f")




