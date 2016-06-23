# defaults settings - for more utility, this was transfered
# to init
import os, inspect

curdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

# DEFINITIONS OF GRIDS OF RELATIVE SPECTRA
# ----------------------------------------------------------------------------------------------------------------------
gridDirectory = os.path.join("/".join(curdir.split('/')[:-1]), 'grids')
# name of the file containing records on synthetic spectra
gridListFile = 'gridlist'

grid_files = dict(
    identification=['DEFAULT', 'OSTAR', 'BSTAR', 'POLLUX', 'AMBRE'],
    directories=[
        ['OSTAR_Z_0.5', 'OSTAR_Z_1.0', 'OSTAR_Z_2.0', 'BSTAR_Z_0.5', 'BSTAR_Z_1.0', 'BSTAR_Z_2.0', 'POLLUX_Z_1.0',
         'AMBRE_Z_1.0'],
        ['OSTAR_Z_0.5', 'OSTAR_Z_1.0', 'OSTAR_Z_2.0'],
        ['BSTAR_Z_0.5', 'BSTAR_Z_1.0', 'BSTAR_Z_2.0'],
        ['POLLUX_Z_1.0'],
        ['AMBRE_Z_1.0']
        ],
    columns=[['FILENAME', 'TEFF', 'LOGG', 'Z'],
             ['FILENAME', 'TEFF', 'LOGG', 'Z'],
             ['FILENAME', 'TEFF', 'LOGG', 'Z'],
             ['FILENAME', 'TEFF', 'LOGG', 'Z'],
             ['FILENAME', 'TEFF', 'LOGG', 'Z']
             ],
    families=[['OSTAR', 'OSTAR', 'OSTAR', 'BSTAR', 'BSTAR', 'BSTAR', 'POLLUX', 'AMBRE'],
              ['OSTAR', 'OSTAR', 'OSTAR'],
              ['BSTAR', 'BSTAR', 'BSTAR'],
              ['POLLUX'],
              ['AMBRE']
              ]
)

# stores default grid order
default_grid_order = ['BSTAR', 'OSTAR', 'AMBRE', 'POLLUX']

# DEFINITIONS OF GRIDS OF ABSOLUTE SPECTRA
# ----------------------------------------------------------------------------------------------------------------------
ABS_gridDirectory = os.path.join("/".join(curdir.split('/')[:-1]), 'grids_ABS')
# name of the file containing records on synthetic spectra
ABS_gridListFile = 'gridlist'

ABS_grid_files = dict(
    identification=['DEFAULT'],
    directories=[
        [],
 ],
    columns=[
        [],
             ],
    families=[
        [],
              ]
)

# stores default grid order
ABS_default_grid_order = []