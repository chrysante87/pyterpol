"""
Test saving of the class.
"""

import pyterpol
rl = pyterpol.RegionList()
rl.add_region(component='primary', identification='r1', wmin=4500., wmax=4600.)
rl.add_region(component='primary', identification='r2', wmin=4500., wmax=4700.)
rl.add_region(component='primary', identification='r3', wmin=4500., wmax=4800.)
rl.add_region(component='secondary', identification='r4', wmin=4500., wmax=4600.)
rl.add_region(component='secondary', identification='r5', wmin=4500., wmax=4800.)
print rl

rl.save('region_list.txt')
rl.clear_all()
print rl

rl.load('region_list.txt')
print rl










