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
    dict(filename='o.asc', component=1, korel=True)
]

# attach the spectra again
ol.add_observations(obs)
print ol

# try some queries - first groups
osl = ol.get_spectra(verbose=True, teff=1)                  # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

osl = ol.get_spectra(verbose=True, rv=0)                    # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

# now wavelengths
osl = ol.get_spectra(verbose=True, wmin=4300)               # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

osl = ol.get_spectra(verbose=True, wmin=4300, wmax=4500)    # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

# osl = ol.get_spectra(verbose=True, wmin=4300, wmax=5000)    # -- correctly incorrect
# print ObservedList(observedSpectraList=osl, debug=True)

# components for a change
osl = ol.get_spectra(verbose=True, component='ALL')         # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

# korel
osl = ol.get_spectra(verbose=True, korel=True)         # -- correct
print ObservedList(observedSpectraList=osl, debug=True)

# and some mixture
osl = ol.get_spectra(verbose=True, component='ALL', rv=1)   # -- correct
print ObservedList(observedSpectraList=osl, debug=True)