"""
Make a comparison with some data created with the old version.
KOREL format of the ouput data is used.
"""

import pyterpol

# 1) setup the model
sl = pyterpol.StarList()
sl.add_component('primary', teff=32554, logg=3.77, vrot=64.5, lr=0.449, rv=6.685, z=1.003)
sl.add_component('secondary', teff=31205, logg=3.36, vrot=213.0, lr=0.551, rv=6.685, z=1.003)
