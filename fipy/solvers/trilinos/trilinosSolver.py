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
    
    def _makeTrilinosMatrix(self, L):
        """ 
        Takes in a pySparse matrix and returns an Epetra.CrsMatrix . 
        Slow, but works. Scales linearly.
        """
        Comm = Epetra.PyComm() 
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
        
        os.remove(filename)

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


        # Here, need to make maps (a root map and one from the matrix)
        
        # Import to the Matrix's map - so that the vectors passed to ApplyTrilinosSolver
        # are ready for TrilinosSolver application
        if isinstance(L, _TrilinosMatrix) and L.parallel:
            # We're in parallel mode
            A.GlobalAssemble()
            DistributedMap = L.map
            if L.comm.MyPID() == 0:
                myElements=A.NumGlobalRows()
            else:
                myElements=0
            RootMap = Epetra.Map(-1, range(0, myElements), 0, L.comm)

            RootToDist = Epetra.Import(DistributedMap, RootMap)

            rLHS = Epetra.Vector(RootMap, x)
            rRHS = Epetra.Vector(RootMap, b)

            LHS = Epetra.Vector(DistributedMap)
            RHS = Epetra.Vector(DistributedMap)
            
            LHS.Import(rLHS, RootToDist, Epetra.Insert)
            RHS.Import(rRHS, RootToDist, Epetra.Insert)

        else:
            # Serially, just cast the vectors from numpy to Trilinos
            RHS = Epetra.Vector(b)
            LHS = Epetra.Vector(x)

        self._applyTrilinosSolver(A, LHS, RHS)


        if isinstance(L, _TrilinosMatrix) and L.parallel:
            # Now, to give each processor a copy of the vector
            
            PersonalMap = Epetra.Map(-1, range(0, A.NumGlobalRows), 0, Comm)
            DistToPers = Epetra.Import(PersonalMap, DistributedMap)

            PersonalLHS = Epetra.Vector(PersonalMap)
            PersonalLHS.Import(LHS, DistToPers, Epetra.Insert)

            x = numerix.array(PersonalLHS)
        else: 
            x = numerix.array(LHS)

    def _getMatrixClass(self):
        return _TrilinosMatrix

    def _applyTrilinosSolver(self):
        pass
