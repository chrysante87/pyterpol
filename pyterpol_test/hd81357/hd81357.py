import numpy as np
import pyterpol

# 1) Generate region
rl = pyterpol.RegionList()
rl.add_region(wmin=6450, wmax=6600)


# 2) Generate components
sl = pyterpol.StarList()

# compute rotational velocity
prot = 33.77409
r_point = 14.68
i = np.radians(62.7)
vsini = 50.57877*r_point*np.sin(i)/prot
print vsini

# add components
sl.add_component(component='secondary', teff=4200., logg=1.86, vrot=vsini)

# construct the interface
itf = pyterpol.Interface(sl=sl, rl=rl)
itf.set_grid_properties(order=2)
itf.setup()

# write the spectra
itf.write_synthetic_spectra()

