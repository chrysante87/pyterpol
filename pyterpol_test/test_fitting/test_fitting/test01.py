"""
Test computation of chi2
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=6320, wmax=6380)
rl.add_region(wmin=6500, wmax=6600)

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=18000., logg=4.5, rv=10., z=1.0, vrot=50.0, lr=0.3)
sl.add_component(component='secondary', teff=25000., logg=4.5, rv=10., z=1.0, vrot=150.0, lr=0.7)

obs = [
    dict(filename='a', error=0.001, group=dict(rv=1)),
    dict(filename='b', error=0.001, group=dict(rv=2)),
    dict(filename='c', error=0.001, group=dict(rv=3))
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

itf = pyterpol.Interface(sl=sl, ol=ol, rl=rl)
itf.setup()

# this reduces the list of observed spectra
reduced = itf.get_comparisons(rv=1)

# setup fitted parameterss
print itf

# this computes the models and chi-square
itf.compute_chi2([
    -100., 100.,
    -100., 100.,
    -100., 100.,
                ])
itf.plot_all_comparisons()

print itf.list_comparisons()
# print itf






