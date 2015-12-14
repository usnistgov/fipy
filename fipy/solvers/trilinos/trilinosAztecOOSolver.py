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

import os

from PyTrilinos import AztecOO

from fipy.solvers.trilinos.trilinosSolver import TrilinosSolver
from fipy.solvers.trilinos.preconditioners.jacobiPreconditioner import JacobiPreconditioner

__all__ = ["TrilinosAztecOOSolver"]

class TrilinosAztecOOSolver(TrilinosSolver):

    """
    .. attention:: This class is abstract, always create on of its subclasses. It provides the code to call all solvers from the Trilinos AztecOO package.

    """

    def __init__(self, tolerance=1e-10, iterations=1000, precon=JacobiPreconditioner()):
        """
        :Parameters:
          - `tolerance`: The required error tolerance.
          - `iterations`: The maximum number of iterative steps to perform.
          - `precon`: Preconditioner object to use.

        """
        if self.__class__ is TrilinosAztecOOSolver:
            raise NotImplementedError, "can't instantiate abstract base class"

        TrilinosSolver.__init__(self, tolerance=tolerance,
                                iterations=iterations, precon=None)
        self.preconditioner = precon

    def _solve_(self, L, x, b):

        Solver = AztecOO.AztecOO(L, x, b)
        Solver.SetAztecOption(AztecOO.AZ_solver, self.solver)

##        Solver.SetAztecOption(AztecOO.AZ_kspace, 30)

        Solver.SetAztecOption(AztecOO.AZ_output, AztecOO.AZ_none)

        if self.preconditioner is not None:
            self.preconditioner._applyToSolver(solver=Solver, matrix=L)
        else:
            Solver.SetAztecOption(AztecOO.AZ_precond, AztecOO.AZ_none)

        output = Solver.Iterate(self.iterations, self.tolerance)

        if self.preconditioner is not None:
            if hasattr(self.preconditioner, 'Prec'):
                del self.preconditioner.Prec

        if 'FIPY_VERBOSE_SOLVER' in os.environ:
            status = Solver.GetAztecStatus()

            from fipy.tools.debug import PRINT
            PRINT('iterations: %d / %d' % (status[AztecOO.AZ_its], self.iterations))
            failure = {AztecOO.AZ_normal : 'AztecOO.AZ_normal',
                       AztecOO.AZ_param : 'AztecOO.AZ_param',
                       AztecOO.AZ_breakdown : 'AztecOO.AZ_breakdown',
                       AztecOO.AZ_loss : 'AztecOO.AZ_loss',
                       AztecOO.AZ_ill_cond : 'AztecOO.AZ_ill_cond',
                       AztecOO.AZ_maxits : 'AztecOO.AZ_maxits'}

            PRINT('failure',failure[status[AztecOO.AZ_why]])

            PRINT('AztecOO.AZ_r:',status[AztecOO.AZ_r])
            PRINT('AztecOO.AZ_scaled_r:',status[AztecOO.AZ_scaled_r])
            PRINT('AztecOO.AZ_rec_r:',status[AztecOO.AZ_rec_r])
            PRINT('AztecOO.AZ_solve_time:',status[AztecOO.AZ_solve_time])
            PRINT('AztecOO.AZ_Aztec_version:',status[AztecOO.AZ_Aztec_version])

        return output
