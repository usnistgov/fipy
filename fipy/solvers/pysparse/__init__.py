from fipy.solvers.pysparse.linearCGSSolver import *
from fipy.solvers.pysparse.linearPCGSolver import *
from fipy.solvers.pysparse.linearGMRESSolver import *
from fipy.solvers.pysparse.linearLUSolver import *
from fipy.solvers.pysparse.linearJORSolver import *

from fipy.solvers.pysparse.preconditioners import *

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearLUSolver
DummySolver = DefaultSolver
GeneralSolver =  LinearLUSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearJORSolver.__all__)
__all__.extend(preconditioners.__all__)
