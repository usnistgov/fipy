import atexit

import pyamgx

from fipy.solvers.pyamgx.pyAMGXSolver import *
from fipy.solvers.pyamgx.linearCGSolver import *
from fipy.solvers.pyamgx.linearGMRESSolver import *
from fipy.solvers.pyamgx.linearFGMRESSolver import *
from fipy.solvers.pyamgx.linearBiCGStabSolver import *
from fipy.solvers.pyamgx.linearLUSolver import *
from fipy.solvers.pyamgx.aggregationAMGSolver import *
from fipy.solvers.pyamgx.classicalAMGSolver import *

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DefaultSolver
GeneralSolver = DefaultSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"
          ]

__all__.extend(linearCGSolver.__all__)
__all__.extend(linearFGMRESSolver.__all__)
__all__.extend(linearBiCGStabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(aggregationAMGSolver.__all__)
