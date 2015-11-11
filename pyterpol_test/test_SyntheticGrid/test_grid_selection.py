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

# check choosing algorithm
#print sygri.get_available_values('teff', logg=2.0, z=2.00)

# test narrowing down the synthetic spectra list
values = sygri.select_parameters(logg=4.0, teff=14000, z=1.5)
for val in values:
    print val
    
# an independent attempt without recursion
for row in sygri.parameterList:
    print row
    
# prints column description
print sygri.columns

print sygri.get_available_values_fast('teff', logg=2.0, z=2.00)
print sygri.get_available_values_fast('logg', z=2.00)