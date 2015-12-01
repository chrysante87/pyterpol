"""
Make a comparison with some data created with the old version.
KOREL format of the ouput data is used.
"""

import pyterpol

# 1) setup the model
sl = pyterpol.StarList()
sl.add_component('primary', teff=32554., logg=3.77, vrot=64.5, lr=0.449, rv=6.685, z=1.003)
sl.add_component('secondary', teff=31205., logg=3.36, vrot=213.0, lr=0.551, rv=6.685, z=1.003)

# 2) setup regions - skip on this, we want to compare only
# rl = pyterpol.RegionList(wmin=4810, wmax=5090)

# 3) setup observed data
ol = pyterpol.ObservedList()
obs = [
    dict(filename='d', error=0.01)
]
ol.add_observations(obs)

# get the interface
itf = pyterpol.Interface(ol=ol, sl=sl)

# set grid properties
itf.set_grid_properties(order=2)

# setut the grid
itf.setup()

# do the comparisons
itf.populate_comparisons()
itf.plot_all_comparisons()
