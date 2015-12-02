import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid(debug=True, mode='default')
parlis = sygri.select_parameters(order=4, **{'logg':3.5, 'z':1.5, 'teff':20000.})
print parlis
parlis = sygri.deselect_exact(parlis,  **{'logg':3.5, 'z':1.5, 'teff':20000.})
print parlis
