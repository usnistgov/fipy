#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "binaryTerm.py"
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

from fipy.terms.binaryTerm import _BinaryTerm
##from fipy.matrices.pysparseMatrix import _CoupledPysparseMeshMatrix
from fipy.variables.coupledCellVariable import _CoupledCellVariable
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

class _CoupledBinaryTerm(_BinaryTerm):
    def __init__(self, term, other):
        _BinaryTerm.__init__(self, term, other)
        if len(self._getVars()) < len(self._getCoupledTerms()):
            raise Exception, 'Different number of solution variables and equations.'
    
    def _getCoupledTerms(self):
        return self.term._getCoupledTerms() + self.other._getCoupledTerms()

    def _verifyVar(self, var):
        if var is not None:
            raise Exception, 'The solution variable should not be specified.'

        if len(self._getVars()) != len(self._getCoupledTerms()):
            raise Exception, 'Different number of solution variables and equations.'

        return _CoupledCellVariable(self._getVars())
    
    def _buildMatrix(self, var, SparseMatrix,  boundaryConditions=(), dt=1.0, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Offset tests

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v0 = CellVariable(mesh=m, value=0.)
        >>> v1 = CellVariable(mesh=m, value=1.)        
        >>> eq0 = TransientTerm(var=v0) - DiffusionTerm(coeff=1., var=v0) - DiffusionTerm(coeff=2., var=v1)
        >>> eq1 = TransientTerm(var=v1) - DiffusionTerm(coeff=3., var=v0) - DiffusionTerm(coeff=4., var=v1) 
        >>> eq = eq0 & eq1
        >>> var = eq._verifyVar(None)
        >>> solver = DefaultSolver()
        >>> var, matrix, RHSvector = eq._buildMatrix(var=var, SparseMatrix=DefaultSolver()._getMatrixClass()) 
        >>> print var.getValue()
        [ 0.  0.  0.  1.  1.  1.]
        >>> print RHSvector.getValue()
        [ 0.  0.  0.  1.  1.  1.]
        >>> print matrix #doctest: +NORMALIZE_WHITESPACE
         2.000000  -1.000000      ---     2.000000  -2.000000      ---
        -1.000000   3.000000  -1.000000  -2.000000   4.000000  -2.000000
            ---    -1.000000   2.000000      ---    -2.000000   2.000000
         3.000000  -3.000000      ---     5.000000  -4.000000      ---
        -3.000000   6.000000  -3.000000  -4.000000   9.000000  -4.000000
            ---    -3.000000   3.000000      ---    -4.000000   5.000000
        
        
        """

        numberOfCells = var.getMesh().getNumberOfCells()
        numberOfVariables = len(self._getVars())
        
        matrix = 0
        RHSvectorsJ = []

        for i, term in enumerate(self._getCoupledTerms()):

            RHSvector = 0
            for j, tmpVar in enumerate(self._getVars()):

                class OffsetSparseMatrix(SparseMatrix):
                    def __init__(self, mesh, bandwidth=0, sizeHint=None, numberOfVariables=numberOfVariables):
                        SparseMatrix.__init__(self, mesh=mesh, bandwidth=bandwidth, sizeHint=sizeHint, numberOfVariables=numberOfVariables)

                    def put(self, vector, id1, id2):
                        SparseMatrix.put(self, vector, id1 + numberOfCells * i, id2 + numberOfCells * j)

                    def addAt(self, vector, id1, id2):
                        SparseMatrix.addAt(self, vector, id1 + numberOfCells * i, id2 + numberOfCells * j)

                    def addAtDiagonal(self, vector):
                        if type(vector) in [type(1), type(1.)]:
                            tmp = numerix.zeros((numberOfCells,), 'd')
                            tmp[:] = vector
                            SparseMatrix.addAtDiagonal(self, tmp)
                        else:
                            SparseMatrix.addAtDiagonal(self, vector)
                            
                tmpVar, tmpMatrix, tmpRHSvector = term._buildMatrix(tmpVar,
                                                                    OffsetSparseMatrix,
                                                                    boundaryConditions=(),
                                                                    dt=dt,
                                                                    transientGeomCoeff=term._getTransientGeomCoeff(tmpVar.getMesh()),
                                                                    diffusionGeomCoeff=term._getDiffusionGeomCoeff(tmpVar.getMesh()))

                RHSvector += tmpRHSvector
                matrix += tmpMatrix

            RHSvectorsJ += [CellVariable(value=RHSvector, mesh=var.getMesh())]

        RHSvector = _CoupledCellVariable(RHSvectorsJ)

	return (var, matrix, RHSvector)

    def __repr__(self):

        return '(' + repr(self.term) + ' & ' + repr(self.other) + ')'

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()

