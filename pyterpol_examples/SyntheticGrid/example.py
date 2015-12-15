"""
This script serves a demonstration of the class SyntheticGrid.
"""
# import the library
import pyterpol
import matplotlib.pyplot as plt

# The handling of the synthetic grid is shadowed from the user,
# therefore the interaction of the user with the grid should
# restrict to only few methods.

# How to create a grid? Ups - we have forgotten, which modes are available.as
# So we create a default grid and have a look at the grids that are available
# In general user should use default grid, because it spans over all
# implemented grids
sg = pyterpol.SyntheticGrid()

# The method list_modes returns string, so it has to be printed
print sg.list_modes()

# Now we know the modes, so we can either create the grid again
sg = pyterpol.SyntheticGrid(mode='bstar', debug=True)

# or just set mode for the existing one - BSTAR will be our
# exemplary grid.
sg.set_mode(mode='bstar')

# we set up a grid, so we can interpolate.
# synthetic spectra should be queried with
# the method get_synthetic_spectrum - lets do some calls
# Grid parameters have to be wrapped using
# following dictionary
pars = dict(teff=18200, logg=4.3, z=1.2)

# We should also pass some boundaries, unless we want
# to get the whole wavelength range of the grid
spectrum1 = sg.get_synthetic_spectrum(pars, [4250, 4500])
print len(spectrum1.wave), len(spectrum1.intens)

# lets get the spectrum - it is type Syntheticspectrum
# so it carries all its features - we can rotate, shrink it
# shift in rv
# we can view properties of the synthetic spectrum
print spectrum1

# we can of course plot_it
spectrum1.plot(savefig=True, figname='spectrum1.png')

# A great feature of the class is that it remembers all
# loaded spectra until the program ends. This means that
# if your nect interpolation requires similar spectra
# from the grid, everything will be much faster
pars = dict(teff=18300, logg=4.2, z=1.1)
spectrum1= sg.get_synthetic_spectrum(pars, [4250, 4500])

# User can also change the resolution of the grid
# by setting keyword step and the number of the
# spectra that are used for interpolation by setting
# keyword order
# step = wavelength step
# order = maximal number of spectra, that should be used for
# interpolation
pars = dict(teff=29300, logg=3.1, z=0.74)
spectrum2 = sg.get_synthetic_spectrum(pars, [4250, 4500], order=4, step=0.05)

# plot comparison of the two spectra
fig = plt.figure()
ax = fig.add_subplot(111)
spectrum1.plot(ax=ax)
spectrum2.plot(ax=ax)
plt.savefig('comparison.png')






