"""
Basic testing of the RegionList class
"""

import pyterpol

# type
rl = pyterpol.RegionList()

# try to add some regions - empty
rl.add_region()
rl.add_region(wmin=4300., wmax=4500., group=0)
print rl

# clear it and star over
rl.clear_all()

rl.add_region(wmin=4200, wmax=4600)
rl.add_region(wmin=4200, wmax=4500, group=10)
rl.clear_all()

# try to set up the regions from an observed file
ol = pyterpol.ObservedList()

# build a list of observations
obs = [
    dict(filename='o.asc', group=dict(rv=0, teff=0)),
    dict(filename='o.asc', group=dict(rv=0, teff=1)),
    dict(filename='o.asc'),
    dict(filename='o.asc', component='primary', korel=True, group=dict(logg=2))
]
ol.add_observations(obs)

# get the spectra wrapped in ObservedSpectrum class
observed = ol.get_spectra()

# setup the regions from spectra
limits = rl.get_regions_from_obs(observed)

