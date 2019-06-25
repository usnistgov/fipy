from __future__ import unicode_literals
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

DefaultSolver = LinearCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver = LinearLUSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"
          ]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

__all__.extend(linearCGSolver.__all__)
__all__.extend(linearFGMRESSolver.__all__)
__all__.extend(linearBiCGStabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(aggregationAMGSolver.__all__)
