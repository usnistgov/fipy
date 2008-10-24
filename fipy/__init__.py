from pkg_resources import get_distribution

FiPy = get_distribution(__name__)

__version__ = FiPy.version
__doc__ = FiPy.get_metadata("PKG-INFO")
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
