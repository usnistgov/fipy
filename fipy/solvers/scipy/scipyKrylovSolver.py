#!/usr/bin/env python

##
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "scipyKrylovSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

import os
import warnings

from fipy.solvers.scipy.scipySolver import _ScipySolver
from fipy.solvers import (MaximumIterationWarning,
                          IllegalInputOrBreakdownWarning)

class _ScipyKrylovSolver(_ScipySolver):
    """
    The base `ScipyKrylovSolver` class.

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def _countIterations(self, xk):
        self.actualIterations += 1
        
    def _solve_(self, L, x, b):
        A = L.matrix
        if self.preconditioner is None:
            M = None
        else:
            M = self.preconditioner._applyToMatrix(A)

        self.actualIterations = 0
        x, info = self.solveFnc(A, b, x,
                                tol=self.tolerance,
                                maxiter=self.iterations,
                                M=M, 
                                callback=self._countIterations)
        self.status['iterations'] = self.actualIterations
        if info == 0:
            self.status['code'] = "Success"
        elif info < 0:
            self.status['code'] = IllegalInputOrBreakdownWarning.__class__.__name__
        elif info > 0:
            self.status['code'] = MaximumIterationWarning.__class__.__name__
            
        self._raiseWarning(info, self.actualIterations, 0.)

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (self.actualIterations, self.iterations))
            
            if info > 0:
                PRINT('tolerance not achieved in {0} iterations'.format(info))
            if info < 0:
                PRINT('illegal input or breakdown: {0}'.format(info))

        return x

    def _raiseWarning(self, info, iter, relres):
        # 0 : successful exit
        # >0 : convergence to tolerance not achieved, number of iterations
        # <0 : illegal input or breakdown

        if info < 0:
            # is stacklevel=5 always what's needed to get to the user's scope?
            warnings.warn(IllegalInputOrBreakdownWarning(self, iter, relres), stacklevel=5)
        elif info > 0:
            warnings.warn(MaximumIterationWarning(self, iter, relres), stacklevel=5)
