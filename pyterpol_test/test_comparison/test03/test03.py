"""
Try to do some comparisons with various systems.
"""

import pyterpol
import matplotlib.pyplot as plt

# ======================== Binary wih KOREL file ===========================

debug=False
# 1) define some observations
ol = pyterpol.ObservedList(debug=debug)
obs = [
       dict(filename='o.asc',  error=0.01, korel=True, component='primary', group={'rv':0}),
       dict(filename='o2.asc', error=0.01, korel=True, component='secondary', group={'rv':1})
       ]
ol.add_observations(obs)

# 2) define fitted regions
rl = pyterpol.RegionList(debug=debug)
rl.add_region(wmin=4300, wmax=4380)
rl.add_region(wmin=4460, wmax=4500)

# 3) define a star
sl = pyterpol.StarList(debug=debug)
sl.add_component(component='primary', teff=20000., logg=4.5, rv=0., vrot=150., lr=0.5, z=1.0)
sl.add_component(component='secondary', teff=10000., logg=4.5, rv=0., vrot=10., lr=0.5, z=1.0)

# 4) define interface
itf = pyterpol.Interface(ol=ol, rl=rl, sl=sl, debug=debug)
itf.setup()
print itf
print itf.rel_rvgroup_region
print itf.list_comparisons()

itf.populate_comparisons()
itf.plot_all_comparisons()
