from linearGMRESSolver import *
from linearCGSSolver import *
from linearPCGSolver import *
from linearLUSolver import *
from linearGeneralSolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearGeneralSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearGeneralSolver.__all__)
