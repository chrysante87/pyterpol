"""
Tryout communication ObservedList -> RegionList -> StarList
"""

import pyterpol

debug = True

# 1) Observations
ol = pyterpol.ObservedList(debug=debug)
obslist = [
    dict(filename='o.asc', component='secondary', error=0.01),
    dict(filename='o.asc', component='primary', error=0.01),
           ]
ol.add_observations(obslist)
print ol

#2) Region
rl = pyterpol.RegionList(debug=debug)

# try to define region from observed
rl.get_regions_from_obs(ol.get_spectra())
print rl


# add regions
rl.add_region(identification='myfirst', wmin=0., wmax=1.)
rl.add_region(identification='mysecond', component='primary', wmin=10., wmax=20.)
rl.add_region(identification='mysecond', component='secondary', wmin=10., wmax=20., groups=dict(teff=1))
print rl


# 3) StarList
sl = pyterpol.StarList(debug=debug)






