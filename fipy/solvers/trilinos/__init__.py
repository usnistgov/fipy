def solverSuite():
    return 'Trilinos'

from PyTrilinos import ML # Gets around strange Trilinos import-order bugs. 

from preconditioners import *

from linearCGSSolver import LinearCGSSolver
from linearPCGSolver import LinearPCGSolver
from linearGMRESSolver import LinearGMRESSolver
from linearLUSolver import LinearLUSolver
from linearBicgstabSolver import LinearBicgstabSolver
