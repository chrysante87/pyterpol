"""
Test computation of chi2 - fitting of 6 RVs on three spectra.
This also shows that simplex wont get out of the local minimum easily,
so more powerful fitting environments are needed if we want
to get over big obstacles.
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=6325, wmax=6375)
rl.add_region(wmin=6540, wmax=6600)

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=17000., logg=4.0, rv=-100.0, z=1.0, vrot=40.0, lr=0.35)
sl.add_component(component='secondary', teff=26000., logg=4.0, rv=100.0, z=1.0, vrot=140.0, lr=0.65)

obs = [
    dict(filename='a', error=0.001, group=dict(rv=1)),
    dict(filename='b', error=0.001, group=dict(rv=2)),
    dict(filename='c', error=0.001, group=dict(rv=3))
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

# setup the class
itf = pyterpol.Interface(sl=sl, ol=ol, rl=rl, debug=False, spectrum_by_spectrum=['rv'])
itf.set_grid_properties(order=2)
itf.setup()
print itf
print itf.list_comparisons()


# setup fitted parameters

itf.set_parameter(parname='logg', component='secondary', fitted=True, vmin=3.5, vmax=4.5)
itf.set_parameter(parname='teff', component='secondary', fitted=True, vmin=20000., vmax=28000.)
itf.set_parameter(parname='vrot', component='secondary', fitted=True, vmin=100., vmax=170.)
itf.set_parameter(parname='rv', component='secondary', fitted=True, vmin=-100., vmax=100.)
itf.set_parameter(parname='lr', component='secondary', fitted=True, vmin=0.5, vmax=0.8)
itf.set_parameter(parname='logg', component='primary', fitted=True, vmin=3.5, vmax=4.5)
itf.set_parameter(parname='teff', component='primary', fitted=True, vmin=15000., vmax=20000.)
itf.set_parameter(parname='vrot', component='primary', fitted=True, vmin=40., vmax=80.)
itf.set_parameter(parname='rv', component='primary', fitted=True, vmin=-120., vmax=120.)
itf.set_parameter(parname='lr', component='primary', fitted=True, vmin=0.2, vmax=0.5)

fitpars = itf.get_fitted_parameters()
#
# # # choose a fitter
itf.choose_fitter('nlopt_nelder_mead', fitparams=fitpars, ftol=1e-4)
print itf

# # first of all reduce the comparison list
l = itf.get_comparisons()
#
# # have a look at the chi^2
init_pars = pyterpol.parlist_to_list(fitpars)
init_chi2 = itf.compute_chi2(init_pars, l=l)
print "Initial settings:",  init_pars, init_chi2
#
# # plot initial comparison
itf.plot_all_comparisons(l=l, figname='initial')
#
# # # do the fitting
itf.run_fit(l=l)
# #
# # # evaluate final parameters
final_pars = pyterpol.parlist_to_list(itf.get_fitted_parameters())
final_chi2 = itf.compute_chi2(final_pars, l=l)
print "Final settings:", final_pars, final_chi2
# # #
# # # 3 plot initial comparison
itf.plot_all_comparisons(l=l, figname='final_spectra')
itf.accept_fit()







