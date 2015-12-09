import pyterpol

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=17000., logg=4.0, rv=-100.0, z=1.0, vrot=40.0, lr=0.35)
sl.add_component(component='secondary', teff=26000., logg=4.0, rv=100.0, z=1.0, vrot=140.0, lr=0.65)

sl.set_parameter(component='primary', group=0, name='rv', value=-123123., vmin=-1e6, vmax=1e6, fitted=True)
sl.save('save_starlist.txt')
sl.clear()
print sl

sl.load('save_starlist.txt')
print sl
