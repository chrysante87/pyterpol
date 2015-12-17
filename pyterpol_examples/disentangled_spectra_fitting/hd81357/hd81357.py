"""
Real life demonstration. HD81357 is an interacting binary.
Its secondary is a Roche-lobe filling star, which is probably
losing its mass. We obtained disentangled spectra of the secondary
in two spectral regions. Here is an estimate of its radiative properties.
"""

# import numpy as np
import pyterpol

## 1) Create RegionList. This step is mandatory.
rl = pyterpol.RegionList()

# Add some regions - two in our case
rl.add_region(wmin=6324, wmax=6424, groups=dict(lr=0))
rl.add_region(wmin=4380, wmax=4497, groups=dict(lr=1))

# 2) Create ObservationList
ol = pyterpol.ObservedList()

# attach the disentangled spectra - in case of KOREL data it is mandatory to specify
# to which component the spectrum belongs and flag it as a KOREL spectrum
obs = [
    dict(filename='DE_blue02_n.dat', component='secondary', korel=True, error=0.01),
    dict(filename='DE_red02_n.dat', component='secondary', korel=True, error=0.01)
]
ol.add_observations(obs)

# 3) Create StarList
sl = pyterpol.StarList()

# add components
# value of each component is passed at the additionb
sl.add_component(component='secondary', teff=4200., logg=1.70, vrot=20., lr=0.2)

# 4) Create the Interface
# Interface serves a wrapper to the three different Lists, fitting environment
# and the synthetic grids. Fitting environment can be added later
itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)

# define grid properties - we want the program to use the
# cubic spline.
itf.set_grid_properties(order=4)

# Once you setup the Interface, you should not change any List
# using anything other than methods defined in the Interface. 
# Doing so, may lead to unpredictable consequences.
itf.setup()

# 4) set the parameters - setting up boundaries is very important, it not only
# speeds up the computation, but also prevents the code from running out of
# the grid
itf.set_parameter(component='secondary', parname='teff', fitted=True, vmin=4005., vmax=6000.)
itf.set_parameter(component='secondary', parname='vrot', fitted=True, vmin=10., vmax=30.)
itf.set_parameter(component='secondary', parname='lr', fitted=True, vmin=0.02, vmax=0.4)
itf.set_parameter(component='secondary', parname='lr', group=1, fitted=True, value=0.10, vmin=0.05, vmax=0.4)
itf.set_parameter(component='secondary', parname='rv', fitted=True, vmin=-20.0, vmax=20.0)

# 6) choose a fitting environment - in this case it is nelder mead and
# the tolerated relative change of chi^2
itf.choose_fitter('nlopt_nelder_mead', ftol=1e-6)

# check that everything is set as intended - we have
# two relative luminosities and two radial velocities
# and that's what we wanted
print itf

"""
==============================================StarList==============================================
Component: secondary
name: rv value: 0.0 vmin: -20.0 vmax: 20.0 fitted: True group: 1 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -20.0 vmax: 20.0 fitted: True group: 2 _typedef: <type 'float'>
name: teff value: 4200.0 vmin: 4005.0 vmax: 6000.0 fitted: True group: 0 _typedef: <type 'float'>
name: vrot value: 20.0 vmin: 10.0 vmax: 30.0 fitted: True group: 0 _typedef: <type 'float'>
name: logg value: 1.7 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.2 vmin: 0.02 vmax: 0.4 fitted: True group: 0 _typedef: <type 'float'>
name: lr value: 0.1 vmin: 0.05 vmax: 0.4 fitted: True group: 1 _typedef: <type 'float'>
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: <type 'float'>
=============================================RegionList=============================================
Region name: region00: (wmin, wmax) = (6324, 6424):
component: all groups: {'lr': 0}
Region name: region01: (wmin, wmax) = (4380, 4497):
component: all groups: {'lr': 1}
============================================ObservedList============================================
List of all attached spectra:
filename: DE_blue02_n.dat component: secondary korel: True loaded: True hasErrors: True global_error: 0.01 group: {'rv': [2]} (min, max): (4377.0, 4500.0)
filename: DE_red02_n.dat component: secondary korel: True loaded: True hasErrors: True global_error: 0.01 group: {'rv': [1]} (min, max): (6321.0, 6426.96)
===============================================Fitter===============================================
Fitter: nlopt_nelder_mead optional_arguments: {'ftol': 1e-06}
Initial parameters:(lr, g.): (0.2, 0); (lr, g.): (0.1, 1); (rv, g.): (0.0, 1); (rv, g.): (0.0, 2); (teff, g.): (4200.0, 0);
(vrot, g.): (20.0, 0);
====================================================================================================

"""

# get the initial chi-square
init_pars = itf.get_fitted_parameters(attribute='value')
init_chi2 = itf.compute_chi2(init_pars)
print "The initial chi-square: %f" % (init_chi2)

"""
The initial chi-square: 20739.073943
"""

# 7) run fitting
itf.run_fit()

# get teh final chi-square
final_pars = itf.get_fitted_parameters(attribute='value')
final_chi2 = itf.compute_chi2(final_pars)
print "The final chi-square: %f" % (final_chi2)

"""
The final chi-square: 5110.224473
"""

# write the fit result
itf.write_fitted_parameters(outputname='result.dat')

# 8) when the fit is done, save the file
itf.save('hd81357.sav')

# 9) Have a look at the comparisons
itf.plot_all_comparisons(figname='finalfit')

# 10) Export the disentangle spectra
itf.write_synthetic_spectra(outputname='final_spectra')

# 11) Have a look how everything converged
itf.plot_convergence(figname='covergence_hd81357.png')

# 12) Have look at uncertainty of the fit
itf.plot_covariances(nbin=20, parameters=['lr', 'teff', 'vrot', 'logg'])





