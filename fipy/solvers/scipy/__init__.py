from linearCGSSolver import *
from linearGMRESSolver import *
from linearBicgstabSolver import *
from linearLUSolver import *
from linearPCGSolver import *

DefaultSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
GeneralSolver = LinearLUSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearCGSSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearBicgstabSolver.__all__)
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
