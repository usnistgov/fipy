import os
execfile(os.path.join(__path__[0], '__version__.py'))

f = file(os.path.join(os.path.split(__path__[0])[0], 'README.txt'))
__doc__ = f.read()
f.close()
del f
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
from preconditioners import *
