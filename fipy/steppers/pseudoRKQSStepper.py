## -*-Pyth-*-
 # ########################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "pseudoRKQSStepper.py"
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

__all__ = ["PseudoRKQSStepper"]

class PseudoRKQSStepper(Stepper):
    """
    Adaptive stepper based on the ``rkqs`` (Runge-Kutta
    "quality-controlled" stepper) algorithm of numerixal Recipes in C: 2nd
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
