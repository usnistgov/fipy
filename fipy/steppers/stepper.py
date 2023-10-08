from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["Stepper"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class Stepper(object):
    """Rudimentary utility class for iterating time steps

    Use `steppyngstounes <https://pages.nist.gov/steppyngstounes/en/latest>`_
    instead.
    """

    def __init__(self, vardata=()):
        self.vardata = vardata

    def sweepFn(vardata, dt, *args, **kwargs):
        residual = 0
        for var, eqn, bcs in vardata:
            residual = max(residual, eqn.sweep(var=var, dt=dt, boundaryConditions=bcs))

        return residual
    sweepFn = staticmethod(sweepFn)

    def successFn(vardata, dt, dtPrev, elapsed, *args, **kwargs):
        pass
    successFn = staticmethod(successFn)

    def failFn(vardata, dt, *args, **kwargs):
        pass
    failFn = staticmethod(failFn)

    def _lowerBound(self, dt):
        dt = max(dt, self.dtMin)
        if self.elapsed + dt == self.elapsed:
            raise FloatingPointError("step size underflow: %g + %g == %g" % (self.elapsed, dt, self.elapsed))

        return dt

    def _step(self, dt, dtPrev, sweepFn, failFn, *args, **kwargs):
        sweepFn(vardata=self.vardata, dt=dt, *args, **kwargs)
        return dt, dt

    def step(self, dt, dtTry=None, dtMin=None, dtPrev=None,
             sweepFn=None, successFn=None, failFn=None, *args, **kwargs):
        sweepFn = sweepFn or self.sweepFn
        successFn = successFn or self.successFn
        failFn = failFn or self.failFn

        dtTry = dtTry or dtMin or dt
        dtPrev = dtPrev or dtMin
        self.dtMin = dtMin or 0.

        self.elapsed = 0.

        while self.elapsed < dt:
            dtMax = dt - self.elapsed
            if dtTry > dtMax:
                dtSave = dtTry
                dtTry = dtMax
            else:
                dtSave = None

            for var, eqn, bcs in self.vardata:
                var.updateOld()

            dtPrev, dtTry = self._step(dt=dtTry, dtPrev=dtPrev,
                                       sweepFn=sweepFn, failFn=failFn,
                                       *args, **kwargs)

            self.elapsed += dtPrev

            successFn(vardata=self.vardata,
                      dtPrev=dtPrev, elapsed=self.elapsed, dt=dt, *args, **kwargs)

            dtTry = max(dtTry, self.dtMin)

        return dtSave or dtPrev, dtSave or dtTry
