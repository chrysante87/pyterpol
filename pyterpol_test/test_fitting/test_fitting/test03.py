"""
Test computation of chi2 - fitting of 6 RVs on three spectra.
This also shows that simplex wont get out of the local minimum easily,
so more powerful fitting environments are needed if we want
to get over big obstacles.
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=6320, wmax=6380)
rl.add_region(wmin=6500, wmax=6600)

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=18000., logg=4.2, rv=-100.0, z=1.0, vrot=50.0, lr=0.3)
sl.add_component(component='secondary', teff=20000., logg=4.1, rv=100.0, z=1.0, vrot=120.0, lr=0.7)

obs = [
    dict(filename='a', error=0.001, group=dict(rv=1)),
    dict(filename='b', error=0.001, group=dict(rv=2)),
    dict(filename='c', error=0.001, group=dict(rv=3))
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

# setup the class
itf = pyterpol.Interface(sl=sl, ol=ol, rl=rl, debug=False)
itf.setup()


# setup fitted parameters
# itf.set_parameter(parname='rv', group=3, fitted=True)
itf.set_parameter(parname='logg', component='secondary', fitted=True)
itf.set_parameter(parname='teff', component='secondary', fitted=True)
itf.set_parameter(parname='vrot', component='secondary', fitted=True)

fitpars = itf.get_fitted_parameters()
print itf
#
# # choose a fitter
itf.choose_fitter('sp_nelder_mead', fitparams=fitpars)
# print itf

# first of all reduce the comparison list
l = itf.get_comparisons(rv=3)

# # have a look at the chi-2
init_pars = pyterpol.parlist_to_list(fitpars)
init_chi2 = itf.compute_chi2(init_pars, l=l)
print "Initial settings:",  init_pars, init_chi2

# plot initial comparison
itf.plot_all_comparisons(l=l, figname='fit_more_spectra')
print itf.list_comparisons(l=l)

# do the fitting
itf.run_fit(l=l)

#evaluate final parameters
final_pars = pyterpol.parlist_to_list(itf.get_fitted_parameters())
final_chi2 = itf.compute_chi2(final_pars, l=l)
print "Final settings:", final_pars, final_chi2

# 3 plot initial comparison
itf.plot_all_comparisons(l=l, figname='final_spectra')
itf.accept_fit()






