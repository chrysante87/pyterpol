"""
Test of setting up the class to start fitting.
"""
import pyterpol

rl = pyterpol.RegionList()
rl.add_region(wmin=5300, wmax=5500)

sl = pyterpol.StarList()
sl.add_component(teff=10000., logg=4.5, rv=10., z=1.0, vrot=20.0)

itf = pyterpol.Interface(sl=sl, rl=rl, debug=False)
print itf

itf.set_parameter(parname='teff', value=20000., vmin=25000., vmax=15000., fitted=True)
itf.set_parameter(parname='logg', value=3.5, vmin=3., vmax=4., fitted=True)
print itf

parlist = itf.get_fitted_parameters()
print pyterpol.parlist_to_list(parlist)

itf.setup()
itf.populate_comparisons()
itf.plot_all_comparisons()

