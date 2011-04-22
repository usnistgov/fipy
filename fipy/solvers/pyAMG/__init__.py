from linearGMRESSolver import LinearGMRESSolver
from linearCGSSolver import LinearCGSSolver
from linearPCGSolver import LinearPCGSolver
from fipy.solvers.scipy.linearLUSolver import LinearLUSolver
from linearGeneralSolver import LinearGeneralSolver

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearGMRESSolver
DummySolver = LinearGMRESSolver
GeneralSolver = LinearGeneralSolver
