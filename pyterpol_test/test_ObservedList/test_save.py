"""
Test of the submodule ObservedList in fitting.
Testing querying of spectra.
"""
import numpy as np
from pyterpol.fitting.interface import ObservedList

# define the type
ol = ObservedList()

# build a list of observations
obs = [
    dict(filename='o.asc', group=dict(rv=0, teff=0)),
    dict(filename='o.asc', group=dict(rv=0, teff=1)),
    dict(filename='o.asc'),
    dict(filename='o.asc', component='primary', korel=True, group=dict(logg=2))
]

# attach the spectra again
ol.add_observations(obs)

# save the class
ol.save('save_ol.txt')