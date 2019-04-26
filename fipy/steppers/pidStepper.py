from __future__ import division
from __future__ import unicode_literals
from fipy.steppers.stepper import Stepper

__all__ = ["PIDStepper"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class PIDStepper(Stepper):
    """
    Adaptive stepper using a PID controller, based on::

        @article{PIDpaper,
           author =  {A. M. P. Valli and G. F. Carey and A. L. G. A. Coutinho},
           title =   {Control strategies for timestep selection in finite element
                      simulation of incompressible flows and coupled
                      reaction-convection-diffusion processes},
           journal = {Int. J. Numer. Meth. Fluids},
           volume =  47,
           year =    2005,
           pages =   {201-231},
        }
    """
    def __init__(self, vardata=(), proportional=0.075, integral=0.175, derivative=0.01):
        Stepper.__init__(self, vardata=vardata)

        self.proportional = proportional
        self.integral = integral
        self.derivative = derivative

        self.error = [1., 1., 1.]
        self.nrej = 0

    def _step(self, dt, dtPrev, sweepFn, failFn, *args, **kwargs):
        while True:
            self.error[2] = sweepFn(vardata=self.vardata, dt=dt, *args, **kwargs)

            # omitting nsa > nsaMax check since it's unclear from
            # the paper what it's supposed to do
            if self.error[2] > 1. and dt > self.dtMin:
                # reject the timestep
                failFn(vardata=self.vardata, dt=dt, *args, **kwargs)

                self.nrej += 1

                for var, eqn, bcs in self.vardata:
                    var.setValue(var.old)

                factor = min(1. / self.error[2], 0.8)

                dt = self._lowerBound(factor * dt)

                dtPrev = dt**2 / dtPrev
            else:
                # step succeeded
                break

        dtNext = dtPrev * ((self.error[1] / self.error[2])**self.proportional
                           * (1. / self.error[2])**self.integral
                           * (self.error[1]**2 / (self.error[2] * self.error[0]))**self.derivative)

        self.error[0] = self.error[1]
        self.error[1] = self.error[2]

        return dt, dtNext
