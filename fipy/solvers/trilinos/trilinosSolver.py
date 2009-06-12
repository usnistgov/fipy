#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosSolver.py"
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

from fipy.solvers.solver import Solver
from fipy.tools.trilinosMatrix import _TrilinosMatrix
from fipy.tools.pysparseMatrix import _PysparseMatrix
from fipy.tools.trilinosMatrix import _numpyToTrilinosVector
from fipy.tools.trilinosMatrix import _trilinosToNumpyVector

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.tools import numerix

class TrilinosSolver(Solver):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    """
    def __init__(self, *args, **kwargs):
        if self.__class__ is TrilinosSolver:
            raise NotImplementedError, "can't instantiate abstract base class"
        else:
            Solver.__init__(self, *args, **kwargs)
            
    def _solve(self):
        mesh = self.var.getMesh()
        comm = Epetra.PyComm()
        
        globalNonOverlappingCellIDs = mesh._getGlobalNonOverlappingCellIDs()
        globalOverlappingCellIDs = mesh._getGlobalOverlappingCellIDs()
        localNonOverlappingCellIDs = mesh._getLocalNonOverlappingCellIDs()
        localOverlappingCellIDs = mesh._getLocalOverlappingCellIDs()
        
        nonOverlappingMap = Epetra.Map(-1, list(globalNonOverlappingCellIDs), 0, comm)
        nonOverlappingVector = Epetra.Vector(nonOverlappingMap, self.var[localNonOverlappingCellIDs])
        
        nonOverlappingRHSvector = Epetra.Vector(nonOverlappingMap, self.RHSvector[localNonOverlappingCellIDs])

        globalMatrix = Epetra.CrsMatrix(Epetra.Copy, nonOverlappingMap, -1)
        
        A = self.matrix[localNonOverlappingCellIDs, localOverlappingCellIDs].matrix
        
        values, irow, jcol = A.find()
        globalMatrix.InsertGlobalValues(globalNonOverlappingCellIDs[irow], 
                                        globalOverlappingCellIDs[jcol], 
                                        values)
        
        globalMatrix.FillComplete()
        globalMatrix.OptimizeStorage()
        
        self._solve_(globalMatrix, nonOverlappingVector, nonOverlappingRHSvector)
        
        overlappingMap =  Epetra.Map(-1, list(globalOverlappingCellIDs), 0, comm)
        overlappingVector = Epetra.Vector(overlappingMap, self.var)
        
        overlappingVector.Import(nonOverlappingVector, Epetra.Import(overlappingMap, nonOverlappingMap), Epetra.Insert)
        
        self.var.setValue(overlappingVector)
            
    def _getMatrixClass(self):
        # an ugly expediency (I blame Wheeler)
        return _PysparseMatrix

