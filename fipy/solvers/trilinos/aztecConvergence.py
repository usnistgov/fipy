from PyTrilinos import AztecOO

from ..convergence import (Convergence, Divergence,
                           BreakdownDivergence,
                           IllConditionedDivergence,
                           IterationDivergence)

class AZ_NormalConvergence(Convergence):
    """
    """
    status_code = AztecOO.AZ_normal
    status_name = "AZ_normal"
    suite = "trilinos"

class AZ_LossOfAccuracyConvergence(Convergence):
    """
    """
    status_code = AztecOO.AZ_loss
    status_name = "AZ_loss"
    suite = "trilinos"

class AZ_ParameterDivergence(Divergence):
    """
    """
    status_code = AztecOO.AZ_param
    status_name = "AZ_param"
    suite = "trilinos"

class AZ_BreakdownDivergence(BreakdownDivergence):
    """
    """
    status_code = AztecOO.AZ_breakdown
    status_name = "AZ_breakdown"
    suite = "trilinos"

class AZ_IllConditionedDivergence(IllConditionedDivergence):
    """
    """
    status_code = AztecOO.AZ_ill_cond
    status_name = "AZ_ill_cond"
    suite = "trilinos"

class AZ_IterationDivergence(IterationDivergence):
    """
    """
    status_code = AztecOO.AZ_maxits
    status_name = "AZ_maxits"
    suite = "trilinos"
