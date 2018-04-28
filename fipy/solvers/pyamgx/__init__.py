import atexit

import pyamgx

from fipy.solvers.pyamgx.pyAMGXSolver import *
from fipy.solvers.pyamgx.linearCGSSolver import *
from fipy.solvers.pyamgx.linearGMRESSolver import *

pyamgx.initialize()
atexit.register(pyamgx.finalize)

DefaultSolver = LinearGMRESSolver

__all__ = ["DefaultSolver",
          ]

__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
