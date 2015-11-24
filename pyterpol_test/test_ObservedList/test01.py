"""
Test of the submodule ObservedList in fitting.
"""
import numpy as np
from pyterpol.fitting.interface import ObservedList

# define the type
ol = ObservedList()

# add some observations
ol.add_observation(filename='o.asc', group=dict(rv=0, teff=0))
ol.add_observation(filename='o.asc', group=dict(rv=0, teff=1))
ol.add_observation(filename='o.asc')
ol.add_observation(filename='o.asc')

# readout_all groups
ol.set_groups()

# list all defined spectra
print ol