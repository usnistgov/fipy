from pyamgx import AMGX_SOLVE_SUCCESS, AMGX_SOLVE_FAILED, AMGX_SOLVE_DIVERGED

from ..convergence import (Convergence, AbsoluteToleranceConvergence,
                           RelativeToleranceConvergence, IterationDivergence,
                           Divergence, BreakdownDivergence,
                           IllConditionedDivergence,
                           PreconditioningDivergence,
                           IllConditionedPreconditionerDivergence,
                           OutOfRangeDivergence)
