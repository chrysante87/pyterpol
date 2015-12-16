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
itf.set_parameter(parname='rv', fitted=True, vmin=-120., vmax=120.)




