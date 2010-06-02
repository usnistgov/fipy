from fipy.tools import parallel

if parallel.Nproc > 1:
    raise Exception("PySparse solvers cannot be used with multiple processors")

from linearCGSSolver import LinearCGSSolver
from linearPCGSolver import LinearPCGSolver
from linearGMRESSolver import LinearGMRESSolver
from linearLUSolver import LinearLUSolver
from linearJORSolver import LinearJORSolver

from preconditioners import *

DefaultSolver = LinearPCGSolver
DefaultAsymmetricSolver = LinearLUSolver
