"""
Test computation of chi2 - fitting of one RV
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=6320, wmax=6380)
rl.add_region(wmin=6500, wmax=6600)

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=18000., logg=4.5, rv=10., z=1.0, vrot=50.0, lr=0.3)
sl.add_component(component='secondary', teff=25000., logg=4.5, rv=10., z=1.0, vrot=150.0, lr=0.7)

obs = [
    dict(filename='a', error=0.001, group=dict(rv=1)),
    dict(filename='b', error=0.001, group=dict(rv=2)),
    dict(filename='c', error=0.001, group=dict(rv=3))
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

# setup the class
itf = pyterpol.Interface(sl=sl, ol=ol, rl=rl)
itf.setup()


# setup fitted parameterss
itf.set_parameter(parname='rv', group=3, fitted=True)
fitpars =  itf.get_fitted_parameters()

# choose a fitter
itf.choose_fitter('np_nelder_mead', fitparams=fitpars)
print itf

# first of all reduce the comparison list
l = itf.get_comparisons(rv=3)

# have a look at the chi-2
init_pars = pyterpol.parlist_to_list(fitpars)
init_chi2 = itf.compute_chi2(init_pars, l=l)
"Initial settings:",  init_pars, init_chi2

# plot initial comparison
itf.plot_all_comparisons(l=l)






