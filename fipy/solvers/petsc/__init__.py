def _dealWithPETScImportPathologies():
    ## The scipy import statement is added to allow PETSc to run
    ## without having sporadic deadlocks. This is caused by weird
    ## behavior in scipy and petsc4py depending on the order in which
    ## modules are imported

    try:
        from scipy import stats
    except:
        pass

# _dealWithPETScImportPathologies()

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
