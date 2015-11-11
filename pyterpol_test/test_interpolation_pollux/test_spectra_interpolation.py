import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol
import time

# create the grid - custom does nothing
sygri = pyterpol.SyntheticGrid('POLLUX', debug=True)

# interpolate at
keys_0 = ['Z', 'LOGG', 'TEFF']
values = [[1.0, 4.0, 12000],[1.0, 4.0, 12200],[1.0, 3.8, 12200]]
for vals in values:    
    props = {}
    for k, v in zip(keys_0, vals):
        props[k] = v
    wmin = 6500
    wmax = 6650
    
    # select parameters
    parlist, vals, keys = sygri.select_and_verify_parameters(order=3, **props)
    
    # select corresponding spectra
    spectra = sygri.get_spectra_for_interpolation(parlist, keys, wmin=wmin, wmax=wmax)
    
    # now does the interpolation
    intens = sygri.interpolate_spectra(parlist, spectra, vals) 
    
    first={}
    for val, key in zip(parlist[0], keys):
        first[key] = val
    # choose one representative 
    wmin, wmax, step  = sygri.get_all(**first)[0].get_size()
    
    # get the wavelength vector
    wave = np.arange(wmin, wmax+step/2., step)
    
    # plot the result
    name = "_".join([key+'_'+str(props[key]) for key in props.keys()])
    plt.plot(wave, intens)
    plt.savefig(name+'.png')
    plt.close()
    
    # save data
    np.savetxt(name+'.dat', np.column_stack([wave, intens]), fmt="%.4f %.6f")

