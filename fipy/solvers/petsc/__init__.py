from fipy.solvers.petsc.linearLUSolver import *
from fipy.solvers.petsc.linearPCGSolver import *
from fipy.solvers.petsc.linearGMRESSolver import *
from fipy.solvers.petsc.linearBicgSolver import *
from fipy.solvers.petsc.linearCGSSolver import *
from fipy.solvers.petsc.dummySolver import *

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = DummySolver
GeneralSolver = DefaultSolver

__all__ = ["DefaultSolver",
           "DummySolver",
           "DefaultAsymmetricSolver",
           "GeneralSolver"]
           
__all__.extend(linearLUSolver.__all__)
__all__.extend(linearPCGSolver.__all__)
__all__.extend(linearGMRESSolver.__all__)
__all__.extend(linearBicgSolver.__all__)
__all__.extend(linearCGSSolver.__all__)
