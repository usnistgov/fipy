from linearCGSSolver import *
from linearPCGSolver import *
from linearGMRESSolver import *
from linearLUSolver import *
from linearJORSolver import *

from preconditioners import *

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
