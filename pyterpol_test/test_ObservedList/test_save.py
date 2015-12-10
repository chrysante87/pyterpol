"""
Test of the observedlist save/load
"""
import numpy as np
from pyterpol.fitting.interface import ObservedList

# define the type
ol = ObservedList()

# build a list of observations
obs = [
    dict(filename='o.asc', group=dict(rv=[0, 1], teff=0), error=0.01),
    dict(filename='o.asc', group=dict(rv=0, teff=1)),
    dict(filename='o.asc'),
    dict(filename='o.asc', component='primary', korel=True, group=dict(logg=2))
]

# attach the spectra again
ol.add_observations(obs)

# save the class
ol.save('save_ol.txt')
print ol
# clear it
ol.clear_all()
print ol

# an
ol.load('save_ol.txt')
print ol