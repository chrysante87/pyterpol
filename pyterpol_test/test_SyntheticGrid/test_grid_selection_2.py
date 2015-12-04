import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid(debug=True, mode='default')
parlis = sygri.select_parameters(order=4, **{'logg':3.9999999999999982, 'z':1.0, 'teff':16000.})
print parlis
parlis = sygri.deselect_exact(parlis,  **{'logg':3.9999999999999982, 'z':1.0, 'teff':16000.})
print parlis
