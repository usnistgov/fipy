from linearCGSSolver import LinearCGSSolver
from linearGMRESSolver import LinearGMRESSolver
from linearBicgstabSolver import LinearBicgstabSolver
from linearLUSolver import LinearLUSolver
from linearPCGSolver import LinearPCGSolver

DefaultSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
GeneralSolver = LinearLUSolver


