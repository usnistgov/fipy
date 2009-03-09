#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosAztecOOSolver.py"
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

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

from PyTrilinos import AztecOO

class TrilinosAztecOOSolver(TrilinosSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses. It provides the code to call all solvers from the Trilinos AztecOO package.

    """
      
    def __init__(self, tolerance=1e-10, iterations=1000, steps=None, precon=JacobiPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `steps`: A deprecated name for `iterations`.
          - `precon`: Preconditioner object to use. 

        """
        if self.__class__ is TrilinosAztecOOSolver:
            raise NotImplementedError, "can't instantiate abstract base class"
            
        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations, steps=steps, precon=None)
        self.preconditioner = precon

    def _applyTrilinosSolver(self, A, LHS, RHS):

        solver = AztecOO.AztecOO(A, LHS, RHS)
        solver.SetAztecOption(AztecOO.AZ_solver, self.solver)

##        solver.SetAztecOption(AztecOO.AZ_kspace, 100)

        solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)

        if self.preconditioner is not None:
            self.preconditioner._applyToSolver(solver=solver, matrix=A)
        else:
            solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)
        
        output = solver.Iterate(self.iterations, self.tolerance)

##         status = solver.GetAztecStatus()

##         print
##         print 'AztecOO.AZ_its:',status[AztecOO.AZ_its]
##         failure = {AztecOO.AZ_normal : 'AztecOO.AZ_normal',
##                    AztecOO.AZ_param : 'AztecOO.AZ_param',
##                    AztecOO.AZ_breakdown : 'AztecOO.AZ_breakdown',
##                    AztecOO.AZ_loss : 'AztecOO.AZ_loss',
##                    AztecOO.AZ_ill_cond : 'AztecOO.AZ_ill_cond',
##                    AztecOO.AZ_maxits : 'AztecOO.AZ_maxits'}
##         print 'stuff',stuff
##         print 'failure',failure[status[AztecOO.AZ_why]]
                                
##         print 'AztecOO.AZ_r:',status[AztecOO.AZ_r]
##         print 'AztecOO.AZ_scaled_r:',status[AztecOO.AZ_scaled_r]
##         print 'AztecOO.AZ_rec_r:',status[AztecOO.AZ_rec_r]
##         print 'AztecOO.AZ_solve_time:',status[AztecOO.AZ_solve_time]
##         print 'AztecOO.AZ_Aztec_version:',status[AztecOO.AZ_Aztec_version]

        return output
