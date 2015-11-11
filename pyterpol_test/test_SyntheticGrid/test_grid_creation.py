import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid()

# this one works - checked
sygri.read_list_from_file('gridlist', columns=['FILENAME', 'TEFF', 'LOGG'], family='BSTAR')
print sygri
print sygri.get_available_values('TEFF')
print sygri.get_available_values('FAMILY')
print sygri.get_available_values('LOGG')
sygri.clear_all()

# this one should also -checked
sygri.read_list_from_file('gridlist2', columns=['FILENAME', 'TEFF', 'LOGG', 'FAMILY'])
#print sygri.SyntheticSpectraList
print sygri.get_available_values('TEFF')
print sygri.get_available_values('FAMILY')
print sygri.get_available_values('LOGG')
sygri.clear_all()

# this one should raise warning - family is assigned automatically - checked
sygri.read_list_from_file('gridlist', columns=['FILENAME', 'TEFF', 'LOGG'])
#print sygri.SyntheticSpectraList
print sygri.get_available_values('TEFF')
print sygri.get_available_values('FAMILY')
print sygri.get_available_values('LOGG')
sygri.clear_all()

# this one should raise error - checked
sygri.read_list_from_file('gridlist', columns=['FILEAME', 'TEFF', 'LOGG'])
#print sygri.SyntheticSpectraList
print sygri.get_available_values('TEFF')
print sygri.get_available_values('FAMILY')
print sygri.get_available_values('LOGG')
sygri.clear_all()