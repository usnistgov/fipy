#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellTerm.py"
 #                                    created: 11/12/03 {11:00:54 AM} 
 #                                last update: 2/25/05 {5:23:49 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.term import Term
from fipy.tools.inline import inline

from fipy.tools.sparseMatrix import SparseMatrix

class _CellTerm(Term):
    def __init__(self, coeff = 1.):
	Term.__init__(self, coeff = coeff)
        self.coeffVectors = None

    def _calcCoeffVectors(self, mesh):
	coeff = self._getGeomCoeff(mesh)
	weight = self._getWeight(mesh)
	
	self.coeffVectors = {
	    'diagonal': coeff * weight['diagonal'],
	    'old value': coeff * weight['old value'],
	    'b vector': coeff * weight['b vector'],
	    'new value': coeff * weight['new value']
	}

    def _getCoeffVectors(self, mesh):
	if self.coeffVectors is None:
	    self._calcCoeffVectors(mesh)
	return self.coeffVectors
	
    def _buildMatrixPy(self, L, oldArray, b, dt, coeffVectors):
        N = len(oldArray)

	b += Numeric.array(oldArray) * coeffVectors['old value'][:] / dt
	b += Numeric.ones([N]) * coeffVectors['b vector'][:]
	L.addAtDiagonal(Numeric.ones([N]) * coeffVectors['new value'][:] / dt)
        L.addAtDiagonal(Numeric.ones([N]) * coeffVectors['diagonal'][:])

    def _buildMatrixIn(self, L, oldArray, b, dt, coeffVectors):
        updatePyArray = Numeric.zeros((oldArray.getMesh().getNumberOfCells()),'d')

        inline._runInlineLoop1("""
            b(i) += oldArray(i) * oldCoeff(i) / dt;
            b(i) += bCoeff(i);
            updatePyArray(i) += newCoeff(i) / dt;
            updatePyArray(i) += diagCoeff(i);
        """,b = b[:],
            oldArray = oldArray.getNumericValue(),
            oldCoeff = Numeric.array(coeffVectors['old value'][:]),
            bCoeff = Numeric.array(coeffVectors['b vector'][:]),
            newCoeff = Numeric.array(coeffVectors['new value'][:]),
            diagCoeff = Numeric.array(coeffVectors['diagonal'][:]),
            updatePyArray = updatePyArray[:],
            ni = len(updatePyArray[:]),
            dt = dt)

	L.addAtDiagonal(updatePyArray)

    def _buildMatrix(self, var, boundaryConditions = (), dt = 1.):
	N = len(var)
	b = Numeric.zeros((N),'d')
	L = SparseMatrix(size = N)
	
	coeffVectors = self._getCoeffVectors(var.getMesh())

	inline._optionalInline(self._buildMatrixIn, self._buildMatrixPy, L, var.getOld(), b, dt, coeffVectors)
	
	return (L, b)


