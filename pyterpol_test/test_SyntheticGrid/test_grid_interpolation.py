import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid()

# read properties
sygri.read_list_from_file('b1.dat', columns=['FILENAME', 'TEFF', 'LOGG', 'Z'], family='BSTAR')
sygri.read_list_from_file('p1.dat', columns=['FILENAME', 'TEFF', 'LOGG', 'Z'], family='POLLUX')

# test of dealing with degeneracies - this should return two pectra
sl = sygri.get_all(z=1.0, teff=15000, logg=4.5)
#print sl

# resolve degeneracy - should raise an exception - checked
#sl = sygri.resolve_degeneracy(sl)

## so we set the grid order and now it should return one spectrum - checked
sygri.set_grid_order(['BSTAR', 'POLLUX'])
sl = sygri.resolve_degeneracy(sl)
#print sl['family']

# this should create a list with intensities of individual
# spectra that will be used for interpolation
parlist, vals, keys = sygri.select_and_verify_parameters(teff=15000, logg=2.75, z=1.5, order=2)
for row in parlist:
    print row
#spectra = sygri.get_spectra_for_interpolation(parlist, ['logg', 'z', 'teff'])
#for spec in spectra[:10]:
    #print spec
    
#print len(parlist), len(spectra)

#try to interpolate the spectra
#print sygri.interpolate_spectra(parlist, spectra, [3.5, 1.0, 15100])

