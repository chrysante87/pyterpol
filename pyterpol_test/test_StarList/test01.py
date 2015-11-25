import pyterpol

sl = pyterpol.StarList()
sl.add_component(component='primary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)
sl.add_component(component='secondary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)
print sl

sl.clear()

# Try not to import one parameter
sl.add_component(component='primary', teff=None, logg=4.5, vrot=0.0, rv=-20.)

# Try not to setup anything
sl.add_component()

# Try to pass non-sencical parameter
sl.add_component(pes='cerny')

print sl