from ..convergence import (Convergence, IterationDivergence,
                           BreakdownDivergence)

class SciPy_Convergence(Convergence):
    """
    """
    status_code = 0
    status_name = "SCIPY_SUCCESS"
    suite = "scipy"

class SciPy_BreakdownDivergence(BreakdownDivergence):
    """
    """
    status_code = -1
    status_name = "SCIPY_ILLEGAL/BREAKDOWN"
    suite = "scipy"

class SciPy_IterationDivergence(IterationDivergence):
    """
    """
    status_code = +1
    status_name = "SCIPY_MAXIT"
    suite = "scipy"
