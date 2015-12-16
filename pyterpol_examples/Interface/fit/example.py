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

# review the loaded interface
print itf


