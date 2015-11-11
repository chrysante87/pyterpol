import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../../')
import pyterpol



# create the grid - custom does nothing
sygri = pyterpol.SyntheticGrid()
print sygri

## create the grid - BSTAR does nothing
sygri = pyterpol.SyntheticGrid('BSTAR')
print sygri

## create the grid - OSTAR does nothing
sygri = pyterpol.SyntheticGrid('OSTAR')
print sygri

## create the grid - POLLUX does nothing
sygri = pyterpol.SyntheticGrid('POLLUX')
print sygri

## create the grid - POLLUX does nothing
sygri = pyterpol.SyntheticGrid('AMBRE')
print sygri

## create the grid - DEFAULT does nothing
sygri = pyterpol.SyntheticGrid('DEFAULT')
print sygri
