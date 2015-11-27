"""
Basic testing (creation and attachement of regions) of the RegionList class
"""

import pyterpol

# type
rl = pyterpol.RegionList(debug=True)

# add regions
# the simplest case
rl.add_region(identification='myfirst', wmin=0., wmax=1.) #- correct
# rl.add_region() #- correctly incorrect

# try to add the same region by name
rl.add_region(identification='myfirst') # - correctly unregistered, because it has the same name
rl.add_region(wmin=0., wmax=1.) # - correctly unregistered, because it is the same region
rl.add_region(component='primary', wmin=0., wmax=1.) # - correctly unregistered, because this region was defined for all

# try to add one component - looks good
rl.add_region(identification='mysecond', component='primary', wmin=10., wmax=20.)
rl.add_region(identification='mysecond', component='secondary', wmin=10., wmax=20.)

# try to pass some groups along with trhe rest
rl.add_region(identification='mythird', wmin=100., wmax=200., groups=dict(teff=1))
rl.add_region(identification='myfourth',component='primary', wmin=100., wmax=200., groups=dict(teff=0))
rl.add_region(identification='myfourth',component='secondary', wmin=100., wmax=200., groups=dict(teff=1))
print rl.get_defined_groups()



