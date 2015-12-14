"""
Demonstration of how to fit KOREL spectra with pyterpol. 
This example shows how to fit disentangled spectra of 
HD81357. There are two spectra available, one in red 
part, the second one the blue part of the visible 
spectrum.
"""

import numpy as np
import pyterpol

## 1) Create RegionList. This step is mandatory.
rl = pyterpol.RegionList()

# Add some regions. 
rl.add_region(wmin=6324, wmax=6424)
rl.add_region(wmin=4380, wmax=4497)

# 2) Create ObservationList
ol = pyterpol.ObservedList()

# attach some data to the list - in case of KOREL data it is mandatory to specify
# to which component the spectrum belongs and flag it as a KOREL spectrum
obs = [
    dict(filename='DE_blue02_n.dat', component='secondary', korel=True, error=0.01),
    dict(filename='DE_red02_n.dat', component='secondary', korel=True, error=0.01)
]
ol.add_observations(obs)

# 3) Create StarList
sl = pyterpol.StarList()

# add components
# one can pass values of each parameter, but only those. Each parameter has 
# many more attributes, which have to be taken care of before the fitting 
# starts.
sl.add_component(component='secondary', teff=4200., logg=1.70, vrot=20., lr=0.2)

# 5) Create the Interface
# Interface serves a wrapper to the three different Lists, fitting environment
# and the synthetic grids.
itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)

# Once you setup the Interface, you should not change any List
# using anything other than methods defined in the Interface. 
# Doing so, may led to unpredictable consequences.
itf.setup()

# 6) set the parameters - setting up boundaries is very important, it not only
# speeeds up the computtation, but also prevents the code from running out of 
# the grid
itf.set_parameter(component='secondary', parname='teff', fitted=True, vmin=4005., vmax=5000.)
# itf.set_parameter(component='secondary', parname='logg', fitted=True, vmin=1.0, vmax=2.5)
itf.set_parameter(component='secondary', parname='vrot', fitted=True, vmin=10., vmax=30.)
itf.set_parameter(component='secondary', parname='lr', fitted=True, vmin=0.05, vmax=0.4)
itf.set_parameter(component='secondary', parname='lr', group=1, fitted=True, value=0.10, vmin=0.05, vmax=0.4)
itf.set_parameter(component='secondary', parname='rv', fitted=True, vmin=-20.0, vmax=20.0)

# 6) choose a fitting environment - in this case it is nelder mead and
# the tolerated relative change of chi^2
fitparams = itf.get_fitted_parameters()
itf.choose_fitter('nlopt_nelder_mead', ftol=1e-6)

# 7) run fitting
itf.run_fit()

# 8) when the fit is done, save the file
# the fit is evaluated in file hd81357_plot.py
itf.save('hd81357.sav')




