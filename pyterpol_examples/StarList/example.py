"""
This script demonstrates capabilities of the StarList class.
"""

import pyterpol
################################## BARE MINIMUM #########################################

# create an empty class
sl = pyterpol.StarList()

# pyterpol knows a set of parameters, which are given to a component
# these parameters are teff, logg, z, lr, vrot and rv. Therefore
# adding component is possible by just calling:
sl.add_component()

# review
print sl

# in general it is better to name the component, because
# than it is easier to identify what belongs to what:
sl.add_component('secondary')

# review
print sl

# and the best is to pass values of all parameters
sl.add_component('tertiary', teff=15000., logg=4.0, vrot=10., rv=-20., lr=0.1, z=1.6)

# review
print sl

# what if we do not want the component to have some parameters?
# just pass None for the parameter and forbade pyterpol
# froim using defaults
sl.add_component('quaternary', teff=None, logg=None, vrot=10., rv=-20., lr=0.1, z=None, use_defaults=False)

# In a rare case, when we want to define a parameter,
# which is not listed among the default ones, we can
# add a parameter using method add_parameter_to_component.
# First one must pass name of the component to which
# we add the data (See why it is better to set your
# component names:-) and after that just pass attributes
# of a parameter.
sl.add_parameter_to_component(component='secondary', name='Stupid_parameter',
                              value=6, unit='half_a_dodo', fitted=False,)

# Nevertheless if you grid is given by parameters not
# listed among the default ones, we encourage you
# to add the parameter to default ones.

#review
print sl

# What if we want to see a list of all defined physical parameters
print sl.get_physical_parameters()

############################### END OF THE SAFE ZONE ####################################

# Here will be demonstration of additional methods, which are not needed for
# usage of the class. It may potentialy be dangerous to use them, if
# your StarList has been added to a interface.

# TODO

