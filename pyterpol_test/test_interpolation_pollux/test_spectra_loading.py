import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol
import time

# create the grid - custom does nothing
sygri = pyterpol.SyntheticGrid('POLLUX')

# try to load a spectrum
keys = ['Z', 'LOGG', 'TEFF']
vals = [1.0, 4.0, 13000]
wmin = 6400
wmax = 6600

# create a dictionary for future look up
specs= {}
for k, v in zip(keys, vals):
    specs[k] = v
    
# load the spectra
t0 = time.time()
print sygri.get_spectra_for_interpolation([vals], keys, wmin=wmin, wmax=wmax)
t1 = time.time()
print t1-t0

# load them again
t0 = time.time()
print sygri.get_spectra_for_interpolation([vals], keys, wmin=wmin, wmax=wmax)
t1 = time.time()
print t1-t0

# check that the spectrum remained loaded - and it did :D
spectrum = sygri.get_all(**specs)[0]
print spectrum.loaded

