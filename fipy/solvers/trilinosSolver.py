#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosSolver.py"
 #                                    created: 06/07/07 
 #                                last update: 06/25/07 
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

from fipy.solvers.solver import Solver
from fipy.tools.trilinosMatrix import _TrilinosMatrix

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.tools import numerix

class TrilinosSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """

    def __init__(self, tolerance=1e-10, iterations=1000, steps=None):
        """
        :Parameters:
        - `tolerance`: The required error tolerance.
        - `iterations`: The maximum number of iterative steps to perform.
        - `steps`: A deprecated name for `iterations`.
    """
        TrilinosSolver.__init__(self, tolerance=tolerance,
                                      iterations=iterations, steps=steps)
        self.preconditioner = None

    
    def _makeTrilinosMatrix(self, L):
        """ 
        Takes in a pySparse matrix and returns an Epetra.CrsMatrix . 
        Slow, but works. Scales linearly.
        """
        Comm = Epetra.PyComm() # For now, no args, Communicator is serial
        m,n = L._getMatrix().shape 

        Map = Epetra.Map(m, 0, Comm)

        #A = Epetra.FECrsMatrix(Epetra.Copy, Map, n)
        #A.InsertGlobalValues(\
        #                Epetra.IntSerialDenseVector(range(0,m)),\
        #                Epetra.IntSerialDenseVector(range(0,n)),\
        #                Epetra.SerialDenseMatrix(L.getNumpyArray()))
        # Replaced with writing to/reading from matrixmarket format temporary file
        
        import tempfile
        import os

        filename = tempfile.mktemp(suffix=".mm")
        L._getMatrix().export_mtx(filename)
        (ierr, A) = EpetraExt.MatrixMarketFileToCrsMatrix(filename, Map)
        
        # File on disk replaced with pipe (NOT REPLACED - NOT WORKING)
        #os.mkfifo(filename)
        #pid = os.fork()     
        #
        #if(pid == 0):
        #    L._getMatrix().export_mtx(filename)
        #    import sys
        #    sys.exit(0)
        #
        #(ierr, A) = EpetraExt.MatrixMarketFileToCrsMatrix(filename, Map)
        #
        #os.waitpid(pid, 0)

        A.FillComplete()
        A.OptimizeStorage()
        return A
    
    def _solve(self, L, x, b):

        if not isinstance(L, _TrilinosMatrix):
            A = self._makeTrilinosMatrix(L)
        else:
            A = L._getMatrix()
            A.FillComplete()
            A.OptimizeStorage()

        RHS = Epetra.Vector(b)
        LHS = Epetra.Vector(x)

        self._applyTrilinosSolver(A, LHS, RHS)
        x = numerix.array(LHS)

    def _getMatrixClass(self):
        return _TrilinosMatrix

    def _applyTrilinosSolver(self):
        pass

    def setPreconditioner(self, preconditioner)
    """
    Define which preconditioner the solver should use.
    Only available for Trilinos solvers.

    :Parameters:
    - `preconditioner`: the preconditioner object to use.
    """
        self.preconditioner = preconditioner
