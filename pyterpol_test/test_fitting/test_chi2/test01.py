"""
Test computation of chi2
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=5300, wmax=5500)
rl.add_region(wmin=6500, wmax=6600)

sl = pyterpol.StarList()
sl.add_component(teff=18000., logg=4.5, rv=10., z=1.0, vrot=50.0, lr=0.3)
sl.add_component(teff=25000., logg=4.5, rv=10., z=1.0, vrot=150.0, lr=0.7)

obs = [
    dict(filename='a', error=0.001),
    dict(filename='b', error=0.001),
    dict(filename='c', error=0.001)
]
ol = pyterpol.ObservedList()
ol.add_observations(obs)

itf = pyterpol.Interface(sl=sl, ol=ol, rl=rl)
itf.setup()

reduced = itf.get_comparisons(rv=1)
print itf.read_chi2_from_comparisons()
print itf.read_chi2_from_comparisons(reduced)




