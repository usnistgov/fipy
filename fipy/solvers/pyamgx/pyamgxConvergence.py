from ..convergence import (
    BreakdownDivergence,
    Convergence,
    Divergence,
    IterationDivergence,
)


class pyamgx_Convergence(Convergence):
    """
    """
    status_code = "success"
    status_name = "AMGX_SOLVE_SUCCESS"
    suite = "pyamgx"

class pyamgx_BreakdownDivergence(BreakdownDivergence):
    """
    """
    status_code = "failed"
    status_name = "AMGX_SOLVE_FAILED"
    suite = "pyamgx"

class pyamgx_IterationDivergence(IterationDivergence):
    """
    """
    status_code = "diverged"
    status_name = "AMGX_SOLVE_DIVERGED"
    suite = "pyamgx"

class pyamgx_Divergence(Divergence):
    """
    """
    status_code = "not_converged"
    status_name = "AMGX_SOLVE_NOT_CONVERGED"
    suite = "pyamgx"
