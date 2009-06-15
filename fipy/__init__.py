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

try:
    from PyTrilinos import Epetra
    
    if Epetra.PyComm().NumProc() > 1:
        raw_input_original = raw_input
        def mpi_raw_input(prompt):
            import sys
            Epetra.PyComm().Barrier()
            sys.stdout.flush()
            if Epetra.PyComm().MyPID() == 0:
                sys.stdout.write(prompt)
                sys.stdout.flush()
                return sys.stdin.readline()
            else:
                return ""
        raw_input = mpi_raw_input

except ImportError:
    pass
