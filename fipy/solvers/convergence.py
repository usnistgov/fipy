from __future__ import division
from __future__ import unicode_literals
from builtins import object
from builtins import str
__docformat__ = 'restructuredtext'

from future.utils import with_metaclass

import logging
import warnings

__all__ = ["ConvergenceBase", "Convergence", "AbsoluteToleranceConvergence",
           "RelativeToleranceConvergence", "Divergence",
           "IterationDivergence", "BreakdownDivergence",
           "IllConditionedDivergence", "OutOfRangeDivergence",
           "PreconditioningDivergence",
           "IllConditionedPreconditionerDivergence"]
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

class ConvergenceBase(with_metaclass(_ConvergenceMeta, object)):
    def __init__(self, solver, iterations, residual, criterion, actual_code=None, **kwargs):
        self.solver = solver
        self.iterations = iterations
        self.residual = residual
        self.criterion = criterion
        if actual_code is None:
            self.actual_code = self.status_code
        else:
            self.actual_code = actual_code

        self.info = vars(self).copy()
        self.info["status_name"] = self.status_name
        self.info["status_code"] = self.status_code
        self.info["max_iterations"] = self.solver.iterations
        self.info.update(kwargs)

        self._log = logging.getLogger(solver.__class__.__module__
                                      + "." + solver.__class__.__name__)

    def warn(self):
        pass

    def __str__(self):
        return str(self.info)

class Convergence(ConvergenceBase):
     message = "User requested convergence criteria is satisfied. Iterations: {0}. Relative error: {1}"

class AbsoluteToleranceConvergence(Convergence):
    pass

class RelativeToleranceConvergence(Convergence):
    pass

class Divergence(ConvergenceBase):

    def warn(self):
        warnings.warn("({status_code}, {status_name}): {residual}".format(**self.info), stacklevel=5)

class IterationDivergence(Divergence):
    pass

class BreakdownDivergence(Divergence):
    pass

class IllConditionedDivergence(Divergence):
    pass

class OutOfRangeDivergence(Divergence):
    pass

class PreconditioningDivergence(Divergence):
    pass

class IllConditionedPreconditionerDivergence(PreconditioningDivergence):
    pass

