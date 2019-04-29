from __future__ import unicode_literals
from fipy.solvers.pyAMG.linearGMRESSolver import *
from fipy.solvers.pyAMG.linearCGSSolver import *
from fipy.solvers.pyAMG.linearPCGSolver import *
from fipy.solvers.pyAMG.linearLUSolver import *
from fipy.solvers.pyAMG.linearGeneralSolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearGeneralSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearGeneralSolver.__all__)
