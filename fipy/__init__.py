from pkg_resources import Requirement, resource_string, get_distribution

FiPy = Requirement.parse("FiPy")

__version__ = get_distribution(FiPy).version
__doc__ = resource_string(FiPy,"README.txt")
__docformat__ = 'restructuredtext'

from boundaryConditions import *
from meshes import *
from solvers import *
from steppers import *
from terms import *
from tools import *
from variables import *
from viewers import *
from models import *
