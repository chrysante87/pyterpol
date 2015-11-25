"""
Test of group assignment,
"""

import pyterpol

# build the class and add some components
sl = pyterpol.StarList()
sl.add_component(component='primary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)
sl.add_component(component='secondary', teff=15000., logg=4.5, vrot=0.0, rv=-20.)
# print sl

# add two groups for rv - cloning works
sl.clone_parameter('primary', 'rv', group=1)
sl.clone_parameter('primary', 'rv', group=2)

# define some groups
groups = dict(primary=dict(rv=[0, 1]), all=dict(rv=[10, 11], teff=[0, 1]))
sl.set_groups(groups)
print sl