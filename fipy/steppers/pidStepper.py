## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "pidStepper.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ########################################################################
 ##

from fipy.steppers.stepper import Stepper

__all__ = ["PIDStepper"]

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
        while 1:
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
