import pyterpol
import matplotlib.pyplot as plt

wmin = 3600
wmax = 4100
sygri = pyterpol.SyntheticGrid(flux_type='absolute')
params = dict(teff=9950, logg=3.7, z=1.0)
spec1 = sygri.get_synthetic_spectrum(params, [wmin, wmax], order=4, step=0.1)
params = dict(teff=10000, logg=3.5, z=1.0)
spec2 = sygri.get_synthetic_spectrum(params, [wmin, wmax], order=4, step=0.1)
params = dict(teff=10000, logg=4.0, z=1.0)
spec3 = sygri.get_synthetic_spectrum(params, [wmin, wmax], order=4, step=0.1)

ax = plt.subplot(111)
spec1.plot(ax=ax)
spec2.plot(ax=ax)
spec3.plot(ax=ax)
plt.show()