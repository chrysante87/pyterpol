"""
Make a comparison with some data created with the old version.
KOREL format of the ouput data is used.
"""

import pyterpol

# 1) setup the model
sl = pyterpol.StarList()
sl.add_component('primary', teff=11257., logg=4.43, vrot=28.8, lr=0.744, rv=-17.94, z=1.000)
sl.add_component('secondary', teff=7714., logg=4.25, vrot=26.42, lr=0.256, rv=-16.73, z=1.000)

# setup the data
obs = [
    dict(filename='output000', korel=True, component='primary',),
    dict(filename='output001', korel=True, component='secondary')
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

# create interface
itf = pyterpol.Interface(ol=ol, sl=sl, debug=True)
itf.setup()

# populate the comparisons
itf.populate_comparisons()

# plot the comparisons
itf.plot_all_comparisons()