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

# check that everything has loaded up
print ol.observedSpectraList['spectrum'][0]

# look at the outcome
print ol.get_defined_groups()