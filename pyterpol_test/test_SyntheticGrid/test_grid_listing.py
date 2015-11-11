import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid()

# read properties
sygri.read_list_from_file('gridlist', columns=['FILENAME', 'TEFF', 'LOGG'], family='BSTAR')

# test narrowing down the synthetic spectra list
for x in sygri.get_all(teff=20000):
    print x
for x in sygri.get_all(teff=23000, logg=4.5):
    print x    
for x in sygri.get_all(teff=25000, logg=6.0):
    print x

# checks narrowing down of the grid with > <   
for x in sygri.narrow_down_grid(teff=20000):
    print x

for x in sygri.narrow_down_grid(teff=25000, logg=(3.0,4.0)):
    print x 
    
for x in sygri.narrow_down_grid(teff=25000, logg=(6.0,7.0)):
    print x   
    