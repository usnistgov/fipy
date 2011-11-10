#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "offsetSparseMatrix.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  summationTerm.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from fipy.tools import numerix

__all__ = ["OffsetSparseMatrix"]

def OffsetSparseMatrix(SparseMatrix, numberOfVariables, numberOfEquations):
    """
    Used in binary terms. equationIndex and varIndex need to be set statically before instantiation.
    """
    
    class OffsetSparseMatrixClass(SparseMatrix):
        equationIndex = 0
        varIndex = 0
        
        def __init__(self, mesh, bandwidth=0, sizeHint=None, 
                     numberOfVariables=numberOfVariables, numberOfEquations=numberOfEquations):
            SparseMatrix.__init__(self, mesh=mesh, bandwidth=bandwidth, sizeHint=sizeHint, 
                                  numberOfVariables=numberOfVariables, numberOfEquations=numberOfEquations) 

        def put(self, vector, id1, id2):
            SparseMatrix.put(self, vector, id1 + self.mesh.numberOfCells * self.equationIndex, id2 + self.mesh.numberOfCells * self.varIndex)

        def addAt(self, vector, id1, id2):
            SparseMatrix.addAt(self, vector, id1 + self.mesh.numberOfCells * self.equationIndex, id2 + self.mesh.numberOfCells * self.varIndex)

        def addAtDiagonal(self, vector):
            if type(vector) in [type(1), type(1.)]:
                tmp = numerix.zeros((self.mesh.numberOfCells,), 'd')
                tmp[:] = vector
                SparseMatrix.addAtDiagonal(self, tmp)
            else:
                SparseMatrix.addAtDiagonal(self, vector)

    return OffsetSparseMatrixClass
