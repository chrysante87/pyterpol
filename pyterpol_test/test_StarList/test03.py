"""
Test of communication between StarList and ObservedList
"""

import pyterpol

# build the class and add some components
sl = pyterpol.StarList()
sl.add_component(component='primary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)
sl.add_component(component='secondary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)

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
sl.set_groups(ol.get_groups_for_components(['primary','secondary', 'ALL']))
sl.delete_hollow_groups()
print sl