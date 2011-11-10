#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearCGSSolver.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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

import sys

from pyamg import smoothed_aggregation_solver
from fipy.solvers.pyAMG.pyAMGSolver import PyAMGSolver

__all__ = ["SmoothedAggregationSolver"]

class SmoothedAggregationSolver(PyAMGSolver):

    def __init__(self, *args, **kwargs):
        super(SmoothedAggregationSolver, self).__init__(*args, **kwargs)

        self.solveFnc = smoothed_aggregation_solver
        self.setupOptionsDict = {"max_coarse": 500}
        self.solveOptionsDict = {"maxiter": self.iterations,
                                 "tol": self.tolerance}

    def _canSolveAssymetric(self):
        return True

