import numpy as np
import pyterpol

# 1) Generate region
rl = pyterpol.RegionList()
rl.add_region(wmin=6324, wmax=6424)
rl.add_region(wmin=4380, wmax=4497)

# 2) Generate observed data
ol = pyterpol.ObservedList()
obs = [
    dict(filename='DE_blue02_n.dat', component='secondary', korel=True, error=0.01),
    dict(filename='DE_red02_n.dat', component='secondary', korel=True, error=0.01)
]
# append the observations
ol.add_observations(obs)

# 3) Generate components
sl = pyterpol.StarList()
# add components
sl.add_component(component='secondary', teff=4200., logg=1.70, vrot=20., lr=0.2)

# 4) construct the interface
itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)
itf.set_grid_properties(order=2)
itf.setup()

# 5) adjust parameters
itf.set_parameter(component='secondary', parname='teff', fitted=True, vmin=4005., vmax=5000.)
itf.set_parameter(component='secondary', parname='logg', fitted=True, vmin=1.0, vmax=2.5)
itf.set_parameter(component='secondary', parname='vrot', fitted=True, vmin=10., vmax=30.)
itf.set_parameter(component='secondary', parname='lr', fitted=True, vmin=0.05, vmax=0.4)
itf.set_parameter(component='secondary', parname='lr', group=1, fitted=True, value=0.10, vmin=0.05, vmax=0.4)
itf.set_parameter(component='secondary', parname='rv', fitted=True, vmin=-20.0, vmax=20.0)

# 6) choose a fitting environment
fitparams = itf.get_fitted_parameters()
itf.choose_fitter('nlopt_nelder_mead', fitparams=fitparams, xtol=1e-4)
itf.populate_comparisons()
print itf.list_comparisons()
l = itf.get_comparisons()

# 7) run fitting
itf.run_fit(l=l)
#
# 8) plot result
print itf
print itf.list_comparisons()
itf.plot_all_comparisons()
itf.write_synthetic_spectra()





