from __future__ import unicode_literals
from fipy.solvers.scipy.linearCGSSolver import *
from fipy.solvers.scipy.linearGMRESSolver import *
from fipy.solvers.scipy.linearBicgstabSolver import *
from fipy.solvers.scipy.linearLUSolver import *
from fipy.solvers.scipy.linearPCGSolver import *

DefaultSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
GeneralSolver = LinearLUSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearBicgstabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
