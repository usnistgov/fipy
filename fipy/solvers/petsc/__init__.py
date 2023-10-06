import logging

_log = logging.getLogger(__name__)

from fipy.solvers.petsc.linearLUSolver import *
from fipy.solvers.petsc.linearPCGSolver import *
from fipy.solvers.petsc.linearGMRESSolver import *
from fipy.solvers.petsc.linearBicgSolver import *
from fipy.solvers.petsc.linearCGSSolver import *
from fipy.solvers.petsc.dummySolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DummySolver
GeneralSolver = DefaultSolver
