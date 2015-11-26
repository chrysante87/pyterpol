"""
Test of communication between StarList and ObservedList
"""

import pyterpol

# build the class and add some components
sl = pyterpol.StarList()
sl.add_component(component='primary', teff=15000., logg=4.5, vrot=0.0, rv=-20., groups=dict(z=1))
sl.add_component(component='secondary', teff=15000., logg=4.5, vrot=0.0, rv=-20., groups=dict(z=0))

# create some data
ol = pyterpol.ObservedList()

# create some observations
obs = [
    dict(filename='o.asc'),
    dict(filename='o.asc'),
    dict(filename='o.asc', group=dict(teff=1))
]

ol.add_observations(obs)

# now query the groups and pass them to starlist
sl.set_groups(ol.get_data_groups(['primary','secondary', 'ALL']))
print sl

# completely overwrite default groups
sl.set_groups(ol.get_data_groups(['primary','secondary', 'ALL']), overwrite=True)
print sl

# return groups common for all components
print "Common groups:", sl.get_common_groups()