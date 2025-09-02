from ..convergence import (Convergence, AbsoluteToleranceConvergence,
                           RelativeToleranceConvergence, RHSZeroConvergence,
                           Divergence, IterationDivergence,
                           PreconditioningDivergence, StagnatedDivergence,
                           IllConditionedDivergence, OutOfRangeDivergence)

class Pysparse_AbsoluteToleranceConvergence(AbsoluteToleranceConvergence):
    """Residual 2-norm less than abstol"""
    status_code = 2
    status_name = "Pysparse_CONVERGED_ATOL"
    suite = "pysparse"

class Pysparse_RHSZeroConvergence(RHSZeroConvergence):
    r""":math:`\vec{b} = 0`, so exact solution is :math:`\vec{x} = 0`.
    """
    status_code = 1
    status_name = "Pysparse_CONVERGED_BZERO"
    suite = "pysparse"

class Pysparse_RelativeToleranceConvergence(RelativeToleranceConvergence):
    """Residual 2-norm decreased by a factor of `rtol`, from 2-norm of right
    hand side.
    """
    status_code = 0
    status_name = "Pysparse_CONVERGED_RTOL"
    suite = "pysparse"

class Pysparse_IterationDivergence(IterationDivergence):
    """Ran out of iterations before any convergence criteria was reached"""
    status_code = -1
    status_name = "Pysparse_DIVERGED_MAXITS"
    suite = "pysparse"

class Pysparse_IllConditionedPreconditionerDivergence(PreconditioningDivergence):
    """The system involving the preconditioner was ill-conditioned.
    """
    status_code = -2
    status_name = "Pysparse_DIVERGED_PC_ILL"
    suite = "pysparse"

class Pysparse_NonPosDefPreconditioningDivergence(PreconditioningDivergence):
    r"""An inner product of the form
    :math:`\mathbf{x}^T \mathbf{K}^{-1} \mathbf{x}` was not positive,
    so the preconditioning matrix :math:`\mathbf{K}` does not appear to be
    positive definite.
    """
    status_code = -3
    status_name = "Pysparse_DIVERGED_PC_NONPOSDEF"
    suite = "pysparse"

class Pysparse_IllConditionedDivergence(IllConditionedDivergence):
    """The matrix appears to be very ill-conditioned.
    """
    status_code = -4
    status_name = "Pysparse_DIVERGED_MAT_ILL"
    suite = "pysparse"

class Pysparse_StagnatedDivergence(StagnatedDivergence):
    """The method stagnated.
    """
    status_code = -5
    status_name = "Pysparse_DIVERGED_STAG"
    suite = "pysparse"

class Pysparse_OutOfRangeDivergence(OutOfRangeDivergence):
    """A scalar quantity became too small or too large to continue computing.
    """
    status_code = -6
    status_name = "Pysparse_DIVERGED_RANGE"
    suite = "pysparse"
