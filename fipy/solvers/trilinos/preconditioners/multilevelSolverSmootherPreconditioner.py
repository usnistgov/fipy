#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "multilevelSolverSmootherPreconditioner.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from PyTrilinos import ML

from fipy.solvers.trilinos.preconditioners.preconditioner import Preconditioner

__all__ = ["MultilevelSolverSmootherPreconditioner"]

class MultilevelSolverSmootherPreconditioner(Preconditioner):
    """
    Multilevel preconditioner for Trilinos solvers using Aztec solvers
    as smoothers.
    
    """
    def __init__(self, levels=10):
        """
        Initialize the multilevel preconditioner

        - `levels`: Maximum number of levels
        """
        self.levels = levels

    def _applyToSolver(self, solver, matrix):
        if matrix.NumGlobalNonzeros() <= matrix.NumGlobalRows():
            return
        
        self.Prec = ML.MultiLevelPreconditioner(matrix, False)
        self.Prec.SetParameterList({"output": 0, "smoother: type" : "Aztec", "smoother: Aztec as solver" : True})
        self.Prec.ComputePreconditioner()
        solver.SetPrecOperator(self.Prec)
        

        
