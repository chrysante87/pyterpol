"""
This tutorial serves as demonstration of how to fit observed
spectra with Pyterpol.

Our observed spectra were created with the old C++ version of the
code. We have three spectra of a binary consisting of
primary: teff = 25000, g = 4.2, , vrot = 150, lr = 0.7, z = 1.0
secondary: teff = 18000, g = 4.2, , vrot = 50, lr = 0.3, z = 1.0
and various radial velocities. They look as if they were
observed spectra.


We have successfully fitted the data with the differential
evolution algorithm form the SciPy library. Our next step is
to get the output from fitting.
"""

import pyterpol

# First load the session
itf = pyterpol.Interface.load('fitted.itf')

# check that everything has loaded correctly
print itf

"""
==============================================StarList==============================================
Component: primary
name: rv value: 49.9857247022 vmin: -120.0 vmax: 120.0 fitted: True group: 1 _typedef: None
name: rv value: 19.9864936135 vmin: -120.0 vmax: 120.0 fitted: True group: 2 _typedef: None
name: rv value: 100.009478284 vmin: -120.0 vmax: 120.0 fitted: True group: 3 _typedef: None
name: teff value: 25000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: None
name: vrot value: 150.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: None
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: None
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: None
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: None
Component: secondary
name: rv value: -49.9460982465 vmin: -120.0 vmax: 120.0 fitted: True group: 1 _typedef: None
name: rv value: -19.9589330606 vmin: -120.0 vmax: 120.0 fitted: True group: 2 _typedef: None
name: rv value: -99.9753261321 vmin: -120.0 vmax: 120.0 fitted: True group: 3 _typedef: None
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
Fitter: sp_diff_evol optional_arguments: {}
Initial parameters:
====================================================================================================
"""

# write a dictionary of parameters and their errors
itf.write_fitted_parameters(outputname='result.dat')

# first we would like to see how our comparisons look like
# naming the figures using 'figname' is not mandatory, nut
# it is advised.
itf.plot_all_comparisons()

# we may want to export the synthetic spectra
# we can write one component in one region -
# this will export a synthetic spectrum for each
# rv_group
itf.write_synthetic_spectra(component='primary', region='region00', outputname='primary')

# or we can write everything together
itf.write_synthetic_spectra()

# convergence can be plotted - by default chi^2.
itf.plot_convergence(figname='convergence_chi.png')

# interface plots the data from the fit log, so it is
# better to save it - also even if our model/fitlog changed
# we can still plot the convergence, stored witihn a fitlog
itf.plot_convergence(parameter='all', f='fit.log', figname='convergence_parameters.png')

# and we can also plot covariance, which will tell us
# what is the uncertainty of the fit - we are interested in rv
# This will plot covariances between rvs for group 1s
itf.plot_covariances(parameters=['rv'], groups=[1], figname='rv_g_1')

# Again it is not necessary to us the registered fitlog
itf.plot_covariances(f='fit.log', parameters=['rv'], groups=[2], figname='rv_g_2')





