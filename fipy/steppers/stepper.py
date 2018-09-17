## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "stepper.py"
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

__docformat__ = 'restructuredtext'

__all__ = ["Stepper"]

class Stepper:
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
