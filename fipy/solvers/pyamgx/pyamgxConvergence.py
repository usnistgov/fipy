from ..convergence import (Convergence, IterationDivergence,
                           BreakdownDivergence)

class pyamgx_Convergence(Convergence):
    """
    """
    status_code = "success"
    status_name = "success"
    suite = "pyamgx"

class pyamgx_BreakdownDivergence(BreakdownDivergence):
    """
    """
    status_code = "failed"
    status_name = "failed"
    suite = "pyamgx"

class pyamgx_IterationDivergence(IterationDivergence):
    """
    """
    status_code = "diverged"
    status_name = "diverged"
    suite = "pyamgx"
