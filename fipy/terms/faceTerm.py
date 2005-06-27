#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceTerm.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 4/4/05 {2:27:08 PM} 
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
 # system.  NIST assumes no responsibility whatsever for its use by
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
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##
 
import Numeric

from fipy.terms.term import Term
import fipy.tools.vector
import fipy.tools.array as array
from fipy.tools.inline import inline
from fipy.tools.sparseMatrix import _SparseMatrix

class _FaceTerm(Term):
    def __init__(self, coeff = 1.):
	Term.__init__(self, coeff = coeff)
        self.coeffMatrix = None
            
    def _getCoeffMatrix(self, mesh, weight):
	coeff = self._getGeomCoeff(mesh)
        if self.coeffMatrix is None:
            self.coeffMatrix = {'cell 1 diag' : coeff * weight['cell 1 diag'],
                                'cell 1 offdiag': coeff * weight['cell 1 offdiag'],
                                'cell 2 diag': coeff * weight['cell 2 diag'],
                                'cell 2 offdiag': coeff * weight['cell 2 offdiag']}
        return self.coeffMatrix

    def _implicitBuildMatrix(self, L, id1, id2, b, weight, mesh, boundaryConditions):
	coeffMatrix = self._getCoeffMatrix(mesh, weight)
	
	interiorFaceIDs = mesh.getInteriorFaceIDs()
	
	L.addAt(array.take(coeffMatrix['cell 1 diag'], interiorFaceIDs),    id1, id1)
	L.addAt(array.take(coeffMatrix['cell 1 offdiag'], interiorFaceIDs), id1, id2)
	L.addAt(array.take(coeffMatrix['cell 2 offdiag'], interiorFaceIDs), id2, id1)
	L.addAt(array.take(coeffMatrix['cell 2 diag'], interiorFaceIDs),    id2, id2)

        N = mesh.getNumberOfCells()
	M = mesh._getMaxFacesPerCell()

        for boundaryCondition in boundaryConditions:
            LL, bb = boundaryCondition._buildMatrix(N, M, coeffMatrix)
            L += LL
            b += bb

    def _explicitBuildMatrix(self, oldArray, id1, id2, b, weight, mesh, boundaryConditions, dt):

	coeffMatrix = self._getCoeffMatrix(mesh, weight)

        inline._optionalInline(self._explicitBuildMatrixIn, self._explicitBuildMatrixPy, oldArray, id1, id2, b, coeffMatrix, mesh, dt)

        N = mesh.getNumberOfCells()
	M = mesh._getMaxFacesPerCell()

        for boundaryCondition in boundaryConditions:

            LL,bb = boundaryCondition._buildMatrix(N, M, coeffMatrix)
            if LL != 0:
##		b -= LL.takeDiagonal() * Numeric.array(oldArray)
                b -= LL * Numeric.array(oldArray)
	    b += bb

    def _explicitBuildMatrixIn(self, oldArray, id1, id2, b, weightedStencilCoeff, mesh, dt):

        oldArrayId1, oldArrayId2 = self._getOldAdjacentValues(oldArray, id1, id2, dt)
	weight = self._getWeight(mesh)['explicit']
        coeff = Numeric.array(self._getGeomCoeff(mesh))
        Nfac = mesh._getNumberOfFaces()

        cell1Diag = Numeric.zeros((Nfac,),'d')
        cell1Diag[:] = weight['cell 1 diag']
        cell1OffDiag = Numeric.zeros((Nfac,),'d')
        cell1OffDiag[:] = weight['cell 1 offdiag']
        cell2Diag = Numeric.zeros((Nfac,),'d')
        cell2Diag[:] = weight['cell 2 diag']
        cell2OffDiag = Numeric.zeros((Nfac,),'d')
        cell2OffDiag[:] = weight['cell 2 offdiag']

##         cell1Diag = Numeric.resize(Numeric.array(weight['cell 1 diag']), (Nfac,))
##         cell1OffDiag = Numeric.resize(Numeric.array(weight['cell 1 offdiag']), (Nfac,))
##         cell2Diag = Numeric.resize(Numeric.array(weight['cell 2 diag']), (Nfac,))
##         cell2OffDiag = Numeric.resize(Numeric.array(weight['cell 2 offdiag']), (Nfac,))

	inline._runInlineLoop1("""
	    long int faceID = faceIDs(i);
	    long int cellID1 = id1(i);
	    long int cellID2 = id2(i);
            
	    b(cellID1) += -coeff(faceID) * (cell1Diag(faceID) * oldArrayId1(i) + cell1OffDiag(faceID) * oldArrayId2(i));
	    b(cellID2) += -coeff(faceID) * (cell2Diag(faceID) * oldArrayId2(i) + cell2OffDiag(faceID) * oldArrayId1(i));
	""",oldArrayId1 = Numeric.array(oldArrayId1),
            oldArrayId2 = Numeric.array(oldArrayId2),
	    id1 = id1,
	    id2 = id2,
	    b = b,
	    cell1Diag = cell1Diag,
	    cell1OffDiag = cell1OffDiag,
	    cell2Diag = cell2Diag,
	    cell2OffDiag = cell2OffDiag,
	    coeff = coeff,
	    faceIDs = mesh.getInteriorFaceIDs(),
	    ni = len(mesh.getInteriorFaceIDs()))

    def _explicitBuildMatrixPy(self, oldArray, id1, id2, b, coeffMatrix, mesh, dt):
        oldArrayId1, oldArrayId2 = self._getOldAdjacentValues(oldArray, id1, id2, dt)

	interiorFaceIDs = mesh.getInteriorFaceIDs()

	cell1diag = array.take(coeffMatrix['cell 1 diag'], interiorFaceIDs)
	cell1offdiag = array.take(coeffMatrix['cell 1 offdiag'], interiorFaceIDs)
	cell2diag = array.take(coeffMatrix['cell 2 diag'], interiorFaceIDs)
	cell2offdiag = array.take(coeffMatrix['cell 2 offdiag'], interiorFaceIDs)

	fipy.tools.vector.putAdd(b, id1, -(cell1diag * oldArrayId1[:] + cell1offdiag * oldArrayId2[:]))
	fipy.tools.vector.putAdd(b, id2, -(cell2diag * oldArrayId2[:] + cell2offdiag * oldArrayId1[:]))

    def _getOldAdjacentValues(self, oldArray, id1, id2, dt):
	return array.take(oldArray, id1), array.take(oldArray, id2)

    def _buildMatrix(self, var, boundaryConditions = (), dt = 1.):
	"""Implicit portion considers
	"""

	mesh = var.getMesh()
	
	id1, id2 = mesh._getAdjacentCellIDs()
	id1 = array.take(id1, mesh.getInteriorFaceIDs())
	id2 = array.take(id2, mesh.getInteriorFaceIDs())
	
        N = len(var)
        b = Numeric.zeros((N),'d')
        L = _SparseMatrix(size = N)

	weight = self._getWeight(mesh)

        if weight.has_key('implicit'):
	    self._implicitBuildMatrix(L, id1, id2, b, weight['implicit'], mesh, boundaryConditions)

        if weight.has_key('explicit'):
            self._explicitBuildMatrix(var.getOld(), id1, id2, b, weight['explicit'], mesh, boundaryConditions, dt)

        return (L, b)

