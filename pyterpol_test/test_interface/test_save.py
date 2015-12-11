"""
Testing the comparison of synthetic and observed spectra.
"""

import pyterpol
import matplotlib.pyplot as plt


debug=False
# 1) define some observations
ol = pyterpol.ObservedList(debug=debug)
obs = [dict(filename='o.asc', error=0.01, component='primary', korel=True),
       dict(filename='o2.asc', error=0.01, component='secondary', korel=True)]
ol.add_observations(obs)

# 2) define fitted regions
rl = pyterpol.RegionList(debug=debug)
rl.add_region(wmin=4300, wmax=4380)
rl.add_region(wmin=4460, wmax=4500, groups={'teff':1})

# 3) define a star
sl = pyterpol.StarList(debug=debug)
sl.add_component(component='primary', teff=20000., logg=4.5, rv=-30., vrot=150., lr=0.5, z=1.0)
sl.add_component(component='secondary', teff=20000., logg=4.5, rv=30., vrot=10., lr=0.5, z=1.0)

# 4) define interface
itf = pyterpol.Interface(ol=ol, rl=rl, sl=sl, debug=False)
itf.setup()
itf.set_parameter(parname='rv', fitted=True)
itf.choose_fitter('nlopt_nelder_mead', ftol=1e-5)
print itf

# 5) save the interface
itf.save('itf_save.txt')
itf.clear_all()
# print itf

# 6) load the interface
itf.load('itf_save.txt')
itf.populate_comparisons()
itf.plot_all_comparisons()
print itf.fitter.fitted_pars_identification

