#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "multilevelNSSAPreconditioner.py"
 #                                    created: 06/25/07
 #                                last update: 11/8/11 {4:09:22 PM}
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Maxsim Gibiansky <maxsim.gibiansky@nist.gov>
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2007-06-25 MLG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelNSSAPreconditioner"]

class MultilevelNSSAPreconditioner(Preconditioner):
    """
    Energy-based minimizing smoothed aggregation suitable for highly
    convective non-symmetric fluid flow problems.
    """
    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return
        
        self.Prec = ML.MultiLevelPreconditioner(matrix, False)

        self.Prec.SetParameterList({"output": 0,
                                    "max levels" : 10,
                                    "prec type" : "MGW",
                                    "increasing or decreasing" : "increasing",
                                    "aggregation: type" : "Uncoupled-MIS",
                                    "energy minimization: enable" : True,
                                    "eigen-analysis: type" : "power-method",
                                    "eigen-analysis: iterations" : 20,
                                    "smoother: sweeps" : 4,
                                    "smoother: damping factor" : 0.67,
                                    "smoother: pre or post" : 'post',
                                    "smoother: type" : "symmetric Gauss-Seidel",
                                    "coarse: type" : 'Amesos-KLU',
                                    "coarse: max size" : 256
                                    })

        self.Prec.ComputePreconditioner()
        
        solver.SetPrecOperator(self.Prec)
        

        
