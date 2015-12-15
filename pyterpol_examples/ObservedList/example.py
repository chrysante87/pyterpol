"""
This example demonstrates how to prepare observations.
"""

import pyterpol
################################## BARE MINIMUM #########################################

# create a blank list
ol = pyterpol.ObservedList()

# now we are ready to attach some data - lets have a look at the data first
# the spectrum is not a KOREL spectrum, so we do not have to pass additional
# information
obs1 = pyterpol.ObservedSpectrum(filename='o.asc')

# lets plot the spectrum
obs1.plot(figname='observed1.png', savefig=True)

# Lets pretend that the second spectrum is a KOREL spectrum,
# because it is not really important, what is what now.
obs2 = pyterpol.ObservedSpectrum(filename='o2.asc', component='primary', korel=True)
obs2.plot(figname='observed2.png', savefig=True)

# Now attach the spectra to the ObservedList one by one
ol.add_one_observation(obs=obs1)
ol.add_one_observation(obs=obs2)

# review the class
print ol

# the name suggests that the spectra can be attached all
# lets clear the list first
ol.clear_all()

# add them all at once
ol.add_observations([obs1, obs2])

# review
print ol

# we saw that pyterpol complains a lot about not having the errors
# of the spectra. We also saw that no groups were assigned. That is
#  because the default groups are set only by the Interface class.
# lets clear the list again
ol.clear_all()

# It is not necessary to wrap the observations into
# the class ObservedSpectrum. ObservedList does that
# for us. We only have to pass the parameters. Also lets
# pass some errors, and some groups
ol.add_one_observation(filename='o.asc', error=0.01, group=dict(rv=1))
ol.add_one_observation(filename='o2.asc', error=0.01, group=dict(rv=2), component='primary', korel=True)

# We can see that groups were set. In this configuration a separate set of radial
# velocities would be fitted for each spectrum. Such configuration is desirable
# if we work with different observed spectra.
print ol

# lets clear the class for the las time
ol.clear_all()

# If our spectra were different regions from one long spectrum,
# we may want to have the velocity same for each spectrum. Lets
# add the observations as a list of dictionaries
obs = [
    dict(filename='o.asc', error=0.01, group=dict(rv=1)),
    dict(filename='o2.asc', error=0.01, group=dict(rv=1), component='primary', korel=True)
]
ol.add_observations(obs)

# in this configuration there will be only one velocity
# for the two spectra. It has to be stresses that although
# two components have the same group for a parameter,
# THE SAME PARAMETER WILL NEVER BE FITTED. EACH COMPONENT
# GETS ALWAYS ITS OWN PARAMETER FOR EACH GROUP.
print ol


############################### END OF THE SAFE ZONE ####################################

# Here will be demonstration of additional methods, which are not needed for
# usage of the class. It may be potentialy be dangerous to use them, if
# your  ObservedList has been added to a interface.