"""
This tutorial serves as demonstration of how to set up an Interface.
Our observed spectra were created with the old C++ version of the
code. We have three spectra of a binary consisting of
primary: teff = 25000, g = 4.2, , vrot = 150, lr = 0.7, z = 1.0
secondary: teff = 18000, g = 4.2, , vrot = 50, lr = 0.3, z = 1.0
and various radial velocities. They look as if they were
observed spectra.

We will make advantage of the default behavior.
"""

import pyterpol

# 1) First we create a starlist
sl = pyterpol.StarList()
sl.add_component(component='primary', teff=25000., logg=4.2, vrot=150., lr=0.7, z=1.0, rv=0.0)
sl.add_component(component='secondary', teff=18000., logg=4.2, vrot=50., lr=0.3, z=1.0, rv=0.0)

# 2) Now think of regions where we might want to do
# the comparison
rl = pyterpol.RegionList()

# the silicon lines
rl.add_region(wmin=6330, wmax=6375)

# Halpha
rl.add_region(wmin=6500, wmax=6600)

# 3) Now attach the data
ol = pyterpol.ObservedList()
obs = [
    dict(filename='a'),
    dict(filename='b'),
    dict(filename='c'),
]
ol.add_observations(obs)

# 4) create the interface
itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)
itf.setup()

# review the class - this is a nice example of how the
# default groups are assigned. Both components have now
# six rvs and two lrs. The explanation is simple - we have
# three observed spectra and two regions. There is one
# relative luminosity for each spectrum. and one for
# each spectrum and each region
print itf

"""
==============================================StarList==============================================
Component: primary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 4 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 5 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 6 _typedef: <type 'float'>
name: teff value: 25000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: <type 'float'>
name: vrot value: 150.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: <type 'float'>
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 1 _typedef: <type 'float'>
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: <type 'float'>
Component: secondary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 4 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 5 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 6 _typedef: <type 'float'>
name: teff value: 18000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: <type 'float'>
name: vrot value: 50.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: <type 'float'>
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.3 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.3 vmin: 0.0 vmax: 1.0 fitted: False group: 1 _typedef: <type 'float'>
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: <type 'float'>
=============================================RegionList=============================================
Region name: region00: (wmin, wmax) = (6330, 6375):
component: all groups: {'lr': 0}
Region name: region01: (wmin, wmax) = (6500, 6600):
component: all groups: {'lr': 1}
============================================ObservedList============================================
List of all attached spectra:
filename: a component: all korel: False loaded: True hasErrors: False global_error: None group: {'rv': [1, 4]} (min, max): (6250.0, 6799.9499999999998)
filename: b component: all korel: False loaded: True hasErrors: False global_error: None group: {'rv': [2, 5]} (min, max): (6250.0, 6799.9499999999998)
filename: c component: all korel: False loaded: True hasErrors: False global_error: None group: {'rv': [3, 6]} (min, max): (6250.0, 6799.9499999999998)
===============================================Fitter===============================================
Fitter: None optional_arguments: {}
Initial parameters:
====================================================================================================
"""

# since our 'observed spectra' are just a model spectra
# the radial velocity and relative luminosity is the
# same for each spectrum, so we might set, that
# relative luminosity is the same for each region
# and radil velocity is the same for each spectrum.

# we have groups for the task - clear the ObservedList
# and RegionList
rl.clear_all()
ol.clear_all()

# add the regions again and set a group in relative luminosity
# for both
rl.add_region(wmin=6330, wmax=6375, groups=dict(lr=0))
rl.add_region(wmin=6500, wmax=6600, groups=dict(lr=0))

# set a radial velocity group for each spectrum
# and add some errors, so we do not have to listen to
# the errros all the time
obs = [
    dict(filename='a', group=dict(rv=1), error=0.001),
    dict(filename='b', group=dict(rv=2), error=0.001),
    dict(filename='c', group=dict(rv=3), error=0.001),
]
ol.add_observations(obs)

# create the Interface again
itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)
itf.setup()

# review - it - we can now see that there is only
# one relative luminosity and three radial velocities
# for each component.
print itf

"""
==============================================StarList==============================================
Component: primary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: <type 'float'>
name: teff value: 25000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: <type 'float'>
name: vrot value: 150.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: <type 'float'>
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.7 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: <type 'float'>
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: <type 'float'>
Component: secondary
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 1 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 2 _typedef: <type 'float'>
name: rv value: 0.0 vmin: -1000.0 vmax: 1000.0 fitted: False group: 3 _typedef: <type 'float'>
name: teff value: 18000.0 vmin: 6000.0 vmax: 50000.0 fitted: False group: 0 _typedef: <type 'float'>
name: vrot value: 50.0 vmin: 0.0 vmax: 500.0 fitted: False group: 0 _typedef: <type 'float'>
name: logg value: 4.2 vmin: 0.0 vmax: 5.0 fitted: False group: 0 _typedef: <type 'float'>
name: lr value: 0.3 vmin: 0.0 vmax: 1.0 fitted: False group: 0 _typedef: <type 'float'>
name: z value: 1.0 vmin: 0.0 vmax: 2.0 fitted: False group: 0 _typedef: <type 'float'>
=============================================RegionList=============================================
Region name: region00: (wmin, wmax) = (6330, 6375):
component: all groups: {'lr': 0}
Region name: region01: (wmin, wmax) = (6500, 6600):
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

# lets save the class - it will create a text file, from
# which the interface can be easily loaded
itf.save('setup.itf')

# and have a look what does the comparisons look like -
# figname serves only as a prefix in this case.
itf.plot_all_comparisons(figname='initial')
