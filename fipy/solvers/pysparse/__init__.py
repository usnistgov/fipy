def solverSuite():
    """
    Returns the identifier of the current solver suite - currently, 'Trilinos' or 'Pysparse'.
    """
    return 'Pysparse'

from linearCGSSolver import LinearCGSSolver
from linearPCGSolver import LinearPCGSolver
from linearGMRESSolver import LinearGMRESSolver
from linearLUSolver import LinearLUSolver
from linearJORSolver import LinearJORSolver
