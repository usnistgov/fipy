from __future__ import division
from __future__ import unicode_literals
from builtins import object
from builtins import str
__docformat__ = 'restructuredtext'

from future.utils import with_metaclass

import json
import logging
import warnings

__all__ = ["ConvergenceBase", "Convergence", "AbsoluteToleranceConvergence",
           "RelativeToleranceConvergence", "RHSZeroConvergence",
           "IterationConvergence", "HappyBreakdownConvergence",
           "IteratingConvergence", "LossOfAccuracyConvergence",
           "Divergence", "IterationDivergence", "BreakdownDivergence",
           "IllConditionedDivergence", "StagnatedDivergence",
           "OutOfRangeDivergence", "PreconditioningDivergence",
           "IllConditionedPreconditionerDivergence", "NullDivergence",
           "ToleranceDivergence"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class _ConvergenceMeta(type):
    # We use __init__ rather than __new__ here because we want
    # to modify attributes of the class *after* they have been
    # created
    # See http://jakevdp.github.io/blog/2012/12/01/a-primer-on-python-metaclasses/
    def __init__(cls, name, bases, dct):
        if not hasattr(cls, "code_registry"):
            # this is the base class, create an empty registry
            # one dictionary for all element types
            cls.code_registry = {}

        if not hasattr(cls, "name_registry"):
            # this is the base class, create an empty registry
            # one dictionary for all element types
            cls.name_registry = {}

        super(_ConvergenceMeta, cls).__init__(name, bases, dct)

        if hasattr(cls, "status_code") and hasattr(cls, "suite"):
            cls.code_registry[(cls.suite, cls.status_code)] = cls

        if hasattr(cls, "status_name") and hasattr(cls, "suite"):
            cls.name_registry[(cls.suite, cls.status_name)] = cls

class DivergenceWarning(UserWarning):
    def __init__(self, divergence):
        msg = "msg={status_name}, code={status_code}, residual={residual}".format(**divergence.info)
        super(DivergenceWarning, self).__init__(msg)

class ConvergenceBase(with_metaclass(_ConvergenceMeta, object)):
    """Information about whether and why a solver converged.

    Attributes
    ----------
    solver : ~fipy.solvers.solver.Solver
        The linear solver that was invoked.
    iterations : int
        The number of linear iterations the solver performed.
    criterion : str
        The :ref:`CONVERGENCE` test used by the solver.
    tolerance_scale : float
        The multiplier applied to the tolerance in order for this solver to
        satisfy `criterion`.
    residual : float
        The unscaled norm of the residual achieved by the solver.
    status_code : int or str
        The canonical return value for this type of convergence.
    status_name : str
        The text representation of `status_code`.
    actual_code : int or str
        The status value actually returned by the solver.
    """

    def __init__(self, solver, iterations, residual, criterion, actual_code=None, **kwargs):
        self.solver = solver
        self.iterations = iterations
        self.criterion = criterion
        self.residual = residual
        if actual_code is None:
            self.actual_code = self.status_code
        else:
            self.actual_code = actual_code

        vars(self).update(kwargs)

        self.log()

    @property
    def info(self):
        info = vars(self).copy()
        info["solver"] = str(info["solver"])
        info.update(vars(self.__class__))
        return {k: v for k, v in info.items() if not k.startswith("_")}

    def log(self, level=logging.DEBUG):
        logger = logging.getLogger(self.__class__.__module__
                                   + "." + self.__class__.__name__)

        logger.log(level, json.dumps(self.info))

    def warn(self):
        pass

    def __str__(self):
        return str(self.info)

class Convergence(ConvergenceBase):
    """Information about why a solver converged.
    """

    message = "User requested convergence criteria is satisfied. Iterations: {0}. Relative error: {1}"

class AbsoluteToleranceConvergence(Convergence):
    """Absolute tolerance satisfied::

       residual < atol * scale
    """
    pass

class RelativeToleranceConvergence(Convergence):
    """Relative tolerance satisfied::

       residual < rtol * scale
    """
    pass

class RHSZeroConvergence(Convergence):
    r""":math:`\vec{b} = 0`, so exact solution is :math:`\vec{x} = 0`.
    """
    pass

class IterationConvergence(Convergence):
    """Requested iterations complete (and no residual calculated).
    """
    pass

class HappyBreakdownConvergence(Convergence):
    '''"Exact" solution found and more iterations will just make things worse.
    '''
    pass

class IteratingConvergence(Convergence):
    """Solve still in progress.
    """
    pass

class LossOfAccuracyConvergence(Convergence):
    """Numerical loss of precision occurred.
    """
    pass

class Divergence(ConvergenceBase):
    """Information about why a solver diverged.
    """

    def warn(self):
        warnings.warn(DivergenceWarning(self), stacklevel=5)

class IterationDivergence(Divergence):
    """Exceeded maximum iterations.
    """
    pass

class BreakdownDivergence(Divergence):
    """Method broke down.
    """
    pass

class IllConditionedDivergence(Divergence):
    """Matrix was ill-conditioned.
    """
    pass

class StagnatedDivergence(Divergence):
    """The method stagnated.
    """
    pass

class OutOfRangeDivergence(Divergence):
    """A value became too small, too large, or invalid.
    """
    pass

class PreconditioningDivergence(Divergence):
    """A problem with the preconditioner.
    """
    pass

class IllConditionedPreconditionerDivergence(PreconditioningDivergence):
    """Preconditioner is ill-conditioned.
    """
    pass

class NullDivergence(Divergence):
    """Breakdown when solving the Hessenberg system within GMRES.
    """
    pass

class ToleranceDivergence(Divergence):
    """Residual norm increased too much.
    """
    pass
