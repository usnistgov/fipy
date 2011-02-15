from linearCGSSolver import LinearCGSSolver
from linearGMRESSolver import LinearGMRESSolver
from linearBicgstabSolver import LinearBicgstabSolver
from linearLUSolver import LinearLUSolver

from fipy.solvers.pysparse.linearPCGSolver import LinearPCGSolver

DefaultSolver = LinearGMRESSolver
DefaultAsymmetricSolver = LinearLUSolver
 
