"""
Test of basic properties of the Parameter class.
"""

from pyterpol.fitting.parameter import Parameter
from pyterpol.fitting.parameter import parameter_definitions as pd

# test that all parameter defined in parameter_definitions result into
# Parameter - passed
for key in pd.keys():
    p = Parameter(**pd[key])

# lets try effective temperature

# Try to do some changes
p['value'] = 12000.
p['vmin'] = 8000.
p['vmax'] = 20000.
p['fitted'] = True

# try to do incorrect change
p['value'] = 'dog'





