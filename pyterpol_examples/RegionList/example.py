"""
This script demonstrates capabilities of the RegionList class.
"""

import pyterpol
################################## BARE MINIMUM #########################################

# create an empty class
rl = pyterpol.RegionList()

# add a region - the simplest way
rl.add_region(wmin=4300, wmax=4500)

# add a region, define name
rl.add_region(wmin=6200, wmax=6600, identification='red')

# for some reason we may want to fit a region only for
# one component. Then we have to specify it. Of course
# the component has to be among those defined in StarList
# to work well later.
rl.add_region(wmin=7600, wmax=7800, identification='nir', component='primary')

# We can now check, that every group was assigned different relative
# luminosity group lr
print rl

# What if we want to fit the same relative luminosity for two regions?
# When setting groups manually you have to be careful, not to assign
# group that is the same as one of those created automatically.
# Automatically they are created 1, 2, 3 - unless one has been efined by user.
rl.add_region(wmin=6340, wmax=6350, identification='SiA', groups=dict(lr=100))
rl.add_region(wmin=6365, wmax=6375, identification='SiB', groups=dict(lr=100))

# Now there will be only one lr group for silicon lines
print rl

# if the user want to get a list of defined regions
print rl.get_registered_regions()

# or a list of wavelength limits
print rl.get_wavelengths()

############################### END OF THE SAFE ZONE ####################################

# Here will be demonstration of additional methods, which are not needed for
# usage of the class. It may be potentialy be dangerous to use them, if
# your RegionList has been added to a interface.

# We may just want to create a comparison of observed and synbthetic
# spectra and we may be very lazy. Then it is possible to read the
# regions from synthetic data.
# TODO