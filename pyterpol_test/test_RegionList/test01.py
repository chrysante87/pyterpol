"""
Basic testing of the RegionList class
"""

import pyterpol

# type
rl = pyterpol.RegionList()

# try to add some regions - empty
rl.add_region()
rl.add_region(wmin=4300., wmax=4500., group=0)
print rl

# clear it and star over
rl.clear_all()
print rl

rl.add_region(wmin=4200, wmax=4600)
rl.add_region(wmin=4200, wmax=4500, group=10)
print rl
