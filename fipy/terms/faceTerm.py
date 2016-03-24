#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "faceTerm.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import os

from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm
from fipy.tools import vector
from fipy.tools import numerix
from fipy.tools import inline

__all__ = ["FaceTerm"]

class FaceTerm(_NonDiffusionTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1., var=None):
        if self.__class__ is FaceTerm:
            raise NotImplementedError("can't instantiate abstract base class")

        _NonDiffusionTerm.__init__(self, coeff=coeff, var=var)
        self.coeffMatrix = None

    def _getCoeffMatrix_(self, var, weight):
        coeff = self._getGeomCoeff(var)

        if self.coeffMatrix is None:
            self.coeffMatrix = {'cell 1 diag' : coeff * weight['cell 1 diag'],
                                'cell 1 offdiag': coeff * weight['cell 1 offdiag'],
                                'cell 2 diag': coeff * weight['cell 2 diag'],
                                'cell 2 offdiag': coeff * weight['cell 2 offdiag']}
        return self.coeffMatrix

    def _implicitBuildMatrix_(self, SparseMatrix, L, id1, id2, b, weight, var, boundaryConditions, interiorFaces, dt):
        mesh = var.mesh
        coeffMatrix = self._getCoeffMatrix_(var, weight)

        id1 = self._reshapeIDs(var, id1)
        id2 = self._reshapeIDs(var, id2)

        L.addAt(numerix.take(coeffMatrix['cell 1 diag'], interiorFaces, axis=-1).ravel(), id1.ravel(), id1.swapaxes(0,1).ravel())
        L.addAt(numerix.take(coeffMatrix['cell 1 offdiag'], interiorFaces, axis=-1).ravel(), id1.ravel(), id2.swapaxes(0,1).ravel())
        L.addAt(numerix.take(coeffMatrix['cell 2 offdiag'], interiorFaces, axis=-1).ravel(), id2.ravel(), id1.swapaxes(0,1).ravel())
        L.addAt(numerix.take(coeffMatrix['cell 2 diag'], interiorFaces, axis=-1).ravel(), id2.ravel(), id2.swapaxes(0,1).ravel())

        N = mesh.numberOfCells
        M = mesh._maxFacesPerCell

        for boundaryCondition in boundaryConditions:
            LL, bb = boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffMatrix)

            if 'FIPY_DISPLAY_MATRIX' in os.environ:
                self._viewer.title = r"%s %s" % (boundaryCondition.__class__.__name__, self.__class__.__name__)
                self._viewer.plot(matrix=LL, RHSvector=bb)
                from fipy import raw_input
                input()

            L += LL
            b += bb

    def _explicitBuildMatrix_(self, SparseMatrix, oldArray, id1, id2, b, weight, var, boundaryConditions, interiorFaces, dt):

        mesh = var.mesh
        coeffMatrix = self._getCoeffMatrix_(var, weight)

        self._explicitBuildMatrixInline_(oldArray=oldArray, id1=id1, id2=id2, b=b, coeffMatrix=coeffMatrix,
                                         mesh=var.mesh, interiorFaces=interiorFaces, dt=dt, weight=weight)

        N = mesh.numberOfCells
        M = mesh._maxFacesPerCell

        for boundaryCondition in boundaryConditions:

            LL,bb = boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffMatrix)
            if LL != 0:
##              b -= LL.takeDiagonal() * numerix.array(oldArray)
                b -= LL * numerix.array(oldArray)
            b += bb

    if inline.doInline:
        def _explicitBuildMatrixInline_(self, oldArray, id1, id2, b, coeffMatrix, mesh, interiorFaces, dt, weight):

            oldArrayId1, oldArrayId2 = self._getOldAdjacentValues(oldArray, id1, id2, dt)
            coeff = numerix.array(self._getGeomCoeff(oldArray))
            Nfac = mesh.numberOfFaces

            cell1Diag = numerix.zeros((Nfac,),'d')
            cell1Diag[:] = weight['cell 1 diag']
            cell1OffDiag = numerix.zeros((Nfac,),'d')
            cell1OffDiag[:] = weight['cell 1 offdiag']
            cell2Diag = numerix.zeros((Nfac,),'d')
            cell2Diag[:] = weight['cell 2 diag']
            cell2OffDiag = numerix.zeros((Nfac,),'d')
            cell2OffDiag[:] = weight['cell 2 offdiag']

            inline._runInline("""
                long int faceID = faceIDs[i];
                long int cellID1 = id1[i];
                long int cellID2 = id2[i];

                b[cellID1] += -coeff[faceID] * (cell1Diag[faceID] * oldArrayId1[i] + cell1OffDiag[faceID] * oldArrayId2[i]);
                b[cellID2] += -coeff[faceID] * (cell2Diag[faceID] * oldArrayId2[i] + cell2OffDiag[faceID] * oldArrayId1[i]);
            """,oldArrayId1 = numerix.array(oldArrayId1),
                oldArrayId2 = numerix.array(oldArrayId2),
                id1 = id1,
                id2 = id2,
                b = b,
                cell1Diag = cell1Diag,
                cell1OffDiag = cell1OffDiag,
                cell2Diag = cell2Diag,
                cell2OffDiag = cell2OffDiag,
                coeff = coeff,
                faceIDs = interiorFaces,
                ni = len(interiorFaces))
    else:
        def _explicitBuildMatrixInline_(self, oldArray, id1, id2, b, coeffMatrix, mesh, interiorFaces, dt, weight):

            oldArrayId1, oldArrayId2 = self._getOldAdjacentValues(oldArray, id1, id2, dt=dt)


            cell1diag = numerix.take(coeffMatrix['cell 1 diag'], interiorFaces)
            cell1offdiag = numerix.take(coeffMatrix['cell 1 offdiag'], interiorFaces)
            cell2diag = numerix.take(coeffMatrix['cell 2 diag'], interiorFaces)
            cell2offdiag = numerix.take(coeffMatrix['cell 2 offdiag'], interiorFaces)

            vector.putAdd(b, id1, -(cell1diag * oldArrayId1 + cell1offdiag * oldArrayId2))
            vector.putAdd(b, id2, -(cell2diag * oldArrayId2 + cell2offdiag * oldArrayId1))

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """Implicit portion considers
        """
        mesh = var.mesh
        id1, id2 = mesh._adjacentCellIDs
        interiorFaces = numerix.nonzero(mesh.interiorFaces)[0]

        id1 = numerix.take(id1, interiorFaces)
        id2 = numerix.take(id2, interiorFaces)

        b = numerix.zeros(var.shape,'d').ravel()
        L = SparseMatrix(mesh=mesh)

        weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)

        if 'implicit' in weight:
            self._implicitBuildMatrix_(SparseMatrix, L, id1, id2, b, weight['implicit'], var, boundaryConditions, interiorFaces, dt)

        if 'explicit' in weight:
            self._explicitBuildMatrix_(SparseMatrix, var.old, id1, id2, b, weight['explicit'], var, boundaryConditions, interiorFaces, dt)

        return (var, L, b)
