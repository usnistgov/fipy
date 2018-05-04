import atexit

import pyamgx

from fipy.solvers.pyamgx.pyAMGXSolver import *
from fipy.solvers.pyamgx.linearLUSolver import *
from fipy.solvers.pyamgx.linearPCGSolver import *
from fipy.solvers.pyamgx.linearGMRESSolver import *

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearGMRESSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"
          ]

__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearLUSolver.__all__)
