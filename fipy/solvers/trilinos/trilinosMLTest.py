#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosMLTest.py"
 #                                    created: 06/07/07 
 #                                last update: 06/11/07 
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
 #  2007-06-04 MLG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import sys

from fipy.solvers.trilinosSolver import TrilinosSolver

try:
    from PyTrilinos import Epetra
    from PyTrilinos import EpetraExt
    from PyTrilinos import Amesos
    from PyTrilinos import AztecOO
    from PyTrilinos import ML
except:
    raise(ImportError, 
          """Epetra, EpetraExt, Amesos, AztecOO, and ML must be available from
             PyTrilinos to use multilevel preconditioners.""")
    
try:
    from PyTrilinos import IFPACK
except:
    pass

from fipy.tools import numerix

class TrilinosMLTest(TrilinosSolver):

    """
    This solver class does not actually solve the system, but outputs 
    information about what ML preconditioner settings will work best.

    Will probably be removed by the time this gets integrated in

    """
    
    def __init__(self, tolerance=1e-10, iterations=5, steps=None, MLOptions={}, ):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterations to perform per test.
          - `steps`: A deprecated name for `iterations`.

          - `MLOptions`: Options to pass to ML. A dictionary of {option:value} 
                         pairs. This will be passed to ML.SetParameterList. 
                         
          For detailed information on the possible parameters for ML, see
          http://trilinos.sandia.gov/packages/ml/documentation.html

          Currently, passing options to Aztec through ML is not supported.
         """
                
        # Former default was MLList = {"max levels" : 4, "smoother: type" : "symmetric Gauss-Seidel", "aggregation: type" : "Uncoupled"} (parameters kept for reference for now)
                
        TrilinosSolver.__init__(self, tolerance=tolerance, 
                                iterations=iterations, steps=steps)

        self.MLOptions = MLOptions
        if not self.MLOptions.has_key("test: max iters"):
            self.MLOptions["test: max iters"] = iterations
        
        if not self.MLOptions.has_key("test: tolerance"):
            self.MLOptions["test: tolerance"] = tolerance

        
    
    def _applyTrilinosSolver(self, A, LHS, RHS):

        Prec = ML.MultiLevelPreconditioner(A, False)
        
        Prec.SetParameterList(self.MLOptions)
        Prec.ComputePreconditioner()
        Prec.AnalyzeSmoothers()
        raw_input("Results of preconditioner test shown above. Press enter to quit.")
        import sys
        sys.exit(0)
