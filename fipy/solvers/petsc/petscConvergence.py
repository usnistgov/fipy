from petsc4py import PETSc

from ..convergence import (Convergence, AbsoluteToleranceConvergence,
                           RelativeToleranceConvergence, IterationDivergence,
                           Divergence, BreakdownDivergence,
                           IllConditionedDivergence,
                           PreconditioningDivergence,
                           IllConditionedPreconditionerDivergence,
                           OutOfRangeDivergence)

# "The values KSP_CONVERGED_CG_NEG_CURVE, KSP_CONVERGED_CG_CONSTRAINED, and
# KSP_CONVERGED_STEP_LENGTH are returned only by the special KSPNASH,
# KSPSTCG, and KSPGLTR solvers which are used by the SNESNEWTONTR (trust
# region) solver."

class KSP_AbsoluteToleranceConvergence(AbsoluteToleranceConvergence):
    """Residual 2-norm less than abstol"""
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_ATOL
    status_name = "KSP_CONVERGED_ATOL"
    suite = "petsc"

class KSP_RelativeToleranceConvergence(RelativeToleranceConvergence):
    """Residual 2-norm decreased by a factor of rtol, from 2-norm of right
    hand side.
    """
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_RTOL
    status_name = "KSP_CONVERGED_RTOL"
    suite = "petsc"

class KSP_NormalAbsoluteToleranceConvergence(KSP_AbsoluteToleranceConvergence):
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_ATOL_NORMAL
    status_name = "KSP_CONVERGED_ATOL_NORMAL"

class KSP_NormalRelativeToleranceConvergence(KSP_RelativeToleranceConvergence):
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_RTOL_NORMAL
    status_name = "KSP_CONVERGED_RTOL_NORMAL"

class KSP_IterationConvergence(Convergence):
    """Used by the KSPPREONLY solver after the single iteration of the
    preconditioner is applied.  Also used when the KSPConvergedSkip()
    convergence test routine is set in KSP.
    """
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_ITS
    status_name = "KSP_CONVERGED_ITS"
    suite = "petsc"

class KSP_HappyBreakdownConvergence(Convergence):
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_HAPPY_BREAKDOWN
    status_name = "KSP_CONVERGED_HAPPY_BREAKDOWN"
    suite = "petsc"

class KSP_IteratingConvergence(Convergence):
    """This flag is returned if you call KSPGetConvergedReason() while the
    KSPSolve() is still running.
    """
    status_code = PETSc.KSP.ConvergedReason.CONVERGED_ITERATING
    status_name = "KSP_CONVERGED_ITERATING"
    suite = "petsc"


class KSP_IterationDivergence(IterationDivergence):
    """Ran out of iterations before any convergence criteria was reached"""
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_MAX_IT
    status_name = "KSP_DIVERGED_ITS"
    suite = "petsc"

class KSP_NullDivergence(Divergence):
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_NULL
    status_name = "KSP_DIVERGED_NULL"
    suite = "petsc"

class KSP_ToleranceDivergence(Divergence):
    """Residual norm increased by a factor of divtol.
    """
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_DTOL
    status_name = "KSP_DIVERGED_DTOL"
    suite = "petsc"

class KSP_BreakdownDivergence(BreakdownDivergence):
    """Generic breakdown in method."""
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_BREAKDOWN
    status_name = "KSP_DIVERGED_BREAKDOWN"
    suite = "petsc"

class KSP_BreakdownBICGDivergence(KSP_BreakdownDivergence):
    """Initial residual is orthogonal to preconditioned initial residual.
    Try a different preconditioner, or a different initial Level.)
    """
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_BREAKDOWN_BICG
    status_name = "KSP_DIVERGED_BREAKDOWN_BICG"

class KSP_IndefiniteMatrixDivergence(IllConditionedDivergence):
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_INDEFINITE_MAT
    status_name = "KSP_DIVERGED_INDEFINITE_MAT"
    suite = "petsc"

class KSP_NonSymmetricDivergence(IllConditionedDivergence):
    """It appears the operator or preconditioner is not symmetric and this
    Krylov method (KSPCG, KSPMINRES, KSPCR) requires symmetry
    """
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_NONSYMMETRIC
    status_name = "KSP_DIVERGED_NONSYMMETRIC"
    suite = "petsc"

class KSP_PreconditioningDivergence(PreconditioningDivergence):
    """It was not possible to build or use the requested preconditioner.
    This is usually due to a zero pivot in a factorization.  It can also
    result from a failure in a subpreconditioner inside a nested
    preconditioner such as PCFIELDSPLIT.
    """
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_PCSETUP_FAILED
    status_name = "KSP_DIVERGED_PC_FAILED"
    suite = "petsc"

class KSP_IndefinitePreconditionerDivergence(IllConditionedPreconditionerDivergence):
    """It appears the preconditioner is indefinite (has both positive and
    negative eigenvalues) and this Krylov method (KSPCG) requires it to be
    positive definite.
    """
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_INDEFINITE_PC
    status_name = "KSP_DIVERGED_INDEFINITE_PC"
    suite = "petsc"

class KSP_NanOrInfDivergence(OutOfRangeDivergence):
    """Residual norm became Not-a-number or Inf likely due to 0/0."""
    status_code = PETSc.KSP.ConvergedReason.DIVERGED_NANORINF
    status_name = "KSP_DIVERGED_NANORINF"
    suite = "petsc"
