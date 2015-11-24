"""
Test of the submodule ObservedList in fitting.
"""
import numpy as np
from pyterpol.fitting.interface import ObservedList

# define the type
ol = ObservedList(debug=True)

# add some observations
ol.add_one_observation(filename='o.asc', group=dict(rv=0, teff=0))
ol.add_one_observation(filename='o.asc', group=dict(rv=0, teff=1))
ol.add_one_observation(filename='o.asc')
ol.add_one_observation(filename='o.asc')

# list the class
print ol

# clear the class
ol.clear_all()

# build a list of observations
obs = [
    dict(filename='o.asc', group=dict(rv=0, teff=0)),
    dict(filename='o.asc', group=dict(rv=0, teff=1)),
    dict(filename='o.asc'),
    dict(filename='o.asc')
]

# attach the spectra again
ol.add_observations(obs)

# list the class
print ol

# check the list of the observed spectra
print ol.observedSpectraList