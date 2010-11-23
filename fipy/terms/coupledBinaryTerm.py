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
from fipy.matrices.pysparseMatrix import _CoupledPysparseMeshMatrix
from fipy.variables.coupledCellVariable import _CoupledCellVariable
from fipy.variables.cellVariable import CellVariable

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

        matricesIJ = []
        RHSvectorsJ = []

        for term in self._getCoupledTerms():

            matricesI = []
            RHSvector = 0
            for tmpVar in self._getVars():
                tmpMatrix = SparseMatrix.__class__(mesh=var.getMesh())
                tmpVar, tmpMatrix, tmpRHSvector = term._buildMatrix(tmpVar,
                                                                    tmpMatrix,
                                                                    boundaryConditions=(),
                                                                    dt=dt,
                                                                    transientGeomCoeff=term._getTransientGeomCoeff(tmpVar.getMesh()),
                                                                    diffusionGeomCoeff=term._getDiffusionGeomCoeff(tmpVar.getMesh()))

                RHSvector += tmpRHSvector
                matricesI += [tmpMatrix]

            RHSvectorsJ += [CellVariable(value=RHSvector, mesh=var.getMesh())]
            matricesIJ += [matricesI]

        matrix = tmpMatrix.getCoupledClass()(mesh=var.getMesh(), matrices=matricesIJ)
        RHSvector = _CoupledCellVariable(RHSvectorsJ)

	return (var, matrix, RHSvector)

    def __repr__(self):

        return '(' + repr(self.term) + ' & ' + repr(self.other) + ')'

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()

