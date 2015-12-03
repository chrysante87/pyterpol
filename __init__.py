# reads basic classes
from .synthetic.makespectrum import SyntheticSpectrum
from .synthetic.makespectrum import SyntheticGrid
from .observed.observations import ObservedSpectrum
from .fitting.interface import ObservedList
from .fitting.interface import StarList
from .fitting.interface import RegionList
from .fitting.interface import Interface
from .fitting.parameter import Parameter
from .fitting.fitter import Fitter

from .synthetic.auxiliary import parlist_to_list
# setup default directories of the grid