from linearCGSSolver import LinearCGSSolver
from linearPCGSolver import LinearPCGSolver
from linearGMRESSolver import LinearGMRESSolver
from linearLUSolver import LinearLUSolver
from linearJORSolver import LinearJORSolver

from preconditioners import *

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearLUSolver
