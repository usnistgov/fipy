from linearCGSSolver import LinearCGSSolver
from linearGMRESSolver import LinearGMRESSolver
from linearBicgstabSolver import LinearBicgstabSolver
from linearLUSolver import LinearLUSolver

DefaultSolver = LinearLUSolver
DummySolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
GeneralSolver = LinearLUSolver
LinearPCGSolver = LinearLUSolver
