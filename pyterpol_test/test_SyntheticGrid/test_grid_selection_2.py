import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol

# create the grid
sygri = pyterpol.SyntheticGrid(debug=True, mode='default')
parlis = sygri.select_parameters(**{'logg':4.25, 'z':1.0, 'teff':7714.})
print parlis

