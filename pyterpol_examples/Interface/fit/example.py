"""
This tutorial serves as demonstration of how to fit observed
spectra with Pyterpol.

Our observed spectra were created with the old C++ version of the
code. We have three spectra of a binary consisting of
primary: teff = 25000, g = 4.2, , vrot = 150, lr = 0.7, z = 1.0
secondary: teff = 18000, g = 4.2, , vrot = 50, lr = 0.3, z = 1.0
and various radial velocities. They look as if they were
observed spectra.

No we will pick up the interface where we ended and
fit the data.
"""

import pyterpol
import numpy as np

# Create the fitting environment and load the last
# session.
itf = pyterpol.Interface.load('setup.itf')

# review the loaded interface - if we compare it with the
# previous example, we see that everything loaded as it
# should
print itf

"""
==============================================StarList==============================================
Component: primary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: None
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: None
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: None
name: teff value: 25000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: None
name: vrot value: 150.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: None
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: None
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: None
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: None
Component: secondary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: None
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: None
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: None
name: teff value: 18000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: None
name: vrot value: 50.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: None
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: None
name: lr value: 0.3 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: None
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: None
=============================================RegionList=============================================
Region name: region00: (wmin, wmax) = (6330.0, 6375.0):
component: all groups: {'lr': 0}
Region name: region01: (wmin, wmax) = (6500.0, 6600.0):
component: all groups: {'lr': 0}
============================================ObservedList============================================
List of all attached spectra:
filename: a component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 1} (min, max): (6250.0, 6799.9499999999998)
filename: b component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 2} (min, max): (6250.0, 6799.9499999999998)
filename: c component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 3} (min, max): (6250.0, 6799.9499999999998)
===============================================Fitter===============================================
Fitter: None optional_arguments: {}
Initial parameters:
====================================================================================================
"""

# the second step is to set what will be fitted. Do not forget to
# set the boundaries vmin, vmax too.

# we can do it parameter by parameter
itf.set_parameter(group=1, component='primary', parname='rv', fitted=True)

# or we can set all at once
itf.set_parameter(group=1, component='primary', parname='rv', fitted=True, vmin=-120., vmax=120.)

# or set everything for primary
itf.set_parameter(component='primary', parname='rv', fitted=True, vmin=-120., vmax=120.)

# or set everything for every rv in the StarList
itf.set_parameter(parname='rv', fitted=True, vmin=-200., vmax=200.)

# now lets set the fitter - one even does not have to construct the class
# it is sufficient to choose the fitter - lets take nelder and mead
# preferably nelder and mead, because they use boundaries
# for simplex it is good to set the initial uncertainty
init_step = 50*np.ones(6)
itf.choose_fitter('nlopt_nelder_mead', init_step=init_step, ftol=1e-6)

# lets review the whole session
print itf

"""
==============================================StarList==============================================
Component: primary
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 1 _typedef: None
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 2 _typedef: None
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 3 _typedef: None
name: teff value: 25000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: None
name: vrot value: 150.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: None
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: None
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: None
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: None
Component: secondary
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 1 _typedef: None
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 2 _typedef: None
name: rv value: 0.0 vmin: -120.0 vmax: 120.0 fitted: True group: 3 _typedef: None
name: teff value: 18000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: None
name: vrot value: 50.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: None
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: None
name: lr value: 0.3 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: None
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: None
=============================================RegionList=============================================
Region name: region00: (wmin, wmax) = (6330.0, 6375.0):
component: all groups: {'lr': 0}
Region name: region01: (wmin, wmax) = (6500.0, 6600.0):
component: all groups: {'lr': 0}
============================================ObservedList============================================
List of all attached spectra:
filename: a component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 1} (min, max): (6250.0, 6799.9499999999998)
filename: b component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 2} (min, max): (6250.0, 6799.9499999999998)
filename: c component: all korel: False loaded: True hasErrors: True global_error: 0.001 group: {'rv': 3} (min, max): (6250.0, 6799.9499999999998)
===============================================Fitter===============================================
Fitter: nlopt_nelder_mead optional_arguments: {}
Initial parameters:(rv, g.): (0.0, 1); (rv, g.): (0.0, 2); (rv, g.): (0.0, 3); (rv, g.): (0.0, 1); (rv, g.): (0.0, 2);
(rv, g.): (0.0, 3);
====================================================================================================
"""

# check the initial chi-square - first we have to get the fitted parameters
# and convert list of Parameters -> list of floats
init_pars =  [par['value'] for par in itf.get_fitted_parameters()]

# or we can let the function do it for us
init_pars = itf.get_fitted_parameters(attribute='value')

# or we can let the function  to do it for us
init_chi2 = itf.compute_chi2(init_pars)
print "Initial chi-square: %f" % init_chi2

# finally run the fitting
itf.run_fit()

# check the final chi-square
final_pars = itf.get_fitted_parameters(attribute='value')
final_chi2 = itf.compute_chi2(final_pars)
print "Final chi-square (nlopt_nelder_mead): %f" % final_chi2

# and plot everything
itf.plot_all_comparisons(figname='final_nm')

# It is not surprising that the fitting failed - Why?!
# for radial velocities one is in general far from the
# global minimum - the variation can ce high, so
# it is better to get the firt estimate with a global method
# like differential evolution
itf.choose_fitter('sp_diff_evol')
itf.run_fit()

# check the final chi-square
final_pars = itf.get_fitted_parameters(attribute='value')
final_chi2 = itf.compute_chi2(final_pars)
print "Final chi-square: %f (sp_diff_evol)" % final_chi2

# lets see the difference
itf.plot_all_comparisons(figname='final_de')

# The message here is that before one really tries
# to fit radiative rpoerties, it is better to do the
# fitting of RVs first. Since we are not using
# any previous information on the RVs (the orbital
# solution is not attached) it is better to
# use global method - especially for large parameter space





