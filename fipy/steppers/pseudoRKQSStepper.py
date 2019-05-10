from __future__ import unicode_literals
from fipy.steppers.stepper import Stepper

__all__ = ["PseudoRKQSStepper"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PseudoRKQSStepper(Stepper):
    """
    Adaptive stepper based on the ``rkqs`` (Runge-Kutta
    "quality-controlled" stepper) algorithm of Numerical Recipes in C: 2nd
    Edition, Section 16.2.

    Not really appropriate, since we're not doing Runge-Kutta steps
    in the first place, but works OK.
    """
    def __init__(self, vardata=(), safety=0.9, pgrow=-0.2, pshrink=-0.25, errcon=1.89e-4):
        Stepper.__init__(self, vardata=vardata)
        self.safety = safety
        self.pgrow = pgrow
        self.pshrink = pshrink
        self.errcon = errcon

    def _step(self, dt, dtPrev, sweepFn, failFn, *args, **kwargs):
        residual = 1e100
        while residual > 1.:
            residual = sweepFn(vardata=self.vardata, dt=dt, *args, **kwargs)

            if residual > 1.:
                # step failed
                failFn(vardata=self.vardata, dt=dt, *args, **kwargs)

                # revert
                for var, eqn, bcs in self.vardata:
                    var.setValue(var.old)

                    dt = max(self.safety * dt * residual**self.pgrow, 0.1 * dt)

                dt = self._lowerBound(dt)

        if residual > self.errcon:
            dtNext = dt * self.safety * residual**self.pshrink
        else:
            dtNext = 5 * dt

        return dt, dtNext
