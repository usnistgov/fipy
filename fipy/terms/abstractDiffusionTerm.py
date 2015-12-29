#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "baseDiffusionTerm.py"
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

__all__ = []

import os

from fipy.terms.unaryTerm import _UnaryTerm
from fipy.tools import numerix
from fipy.terms import TermMultiplyError
from fipy.terms import AbstractBaseClassError
from fipy.variables.faceVariable import FaceVariable

class _AbstractDiffusionTerm(_UnaryTerm):

    def __init__(self, coeff = (1.,), var=None):
        if self.__class__ is _AbstractDiffusionTerm:
            raise AbstractBaseClassError

        if type(coeff) not in (type(()), type([])):
            coeff = [coeff]

        self.order = len(coeff) * 2


        if len(coeff) > 0:
            self.nthCoeff = coeff[0]

            from fipy.variables.variable import Variable
            if not isinstance(self.nthCoeff, Variable):
                self.nthCoeff = Variable(value = self.nthCoeff)

            from fipy.variables.cellVariable import CellVariable
            if isinstance(self.nthCoeff, CellVariable):
                self.nthCoeff = self.nthCoeff.arithmeticFaceValue

        else:
            self.nthCoeff = None

        _UnaryTerm.__init__(self, coeff=coeff, var=var)

        if self.order > 0:
            self.lowerOrderDiffusionTerm = self.__class__(coeff = coeff[1:])

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            self.coeff[0] = other * self.coeff[0]
            return self.__class__(coeff=self.coeff, var=self.var)
        else:
            raise TermMultiplyError

    __rmul__ = __mul__

    def __neg__(self):
        negatedCoeff = list(self.coeff)

        from fipy.variables.variable import Variable
        if isinstance(negatedCoeff[0], (list, tuple)):
            negatedCoeff[0] = -numerix.array(negatedCoeff[0])
        else:
            negatedCoeff[0] = -negatedCoeff[0]

        return self.__class__(coeff=negatedCoeff, var=self.var)

    def __getBoundaryConditions(self, boundaryConditions):
        higherOrderBCs = []
        lowerOrderBCs = []

        for bc in boundaryConditions:
            bcDeriv = bc._getDerivative(self.order - 2)
            if bcDeriv:
                higherOrderBCs.append(bcDeriv)
            else:
                lowerOrderBCs.append(bc)

        return higherOrderBCs, lowerOrderBCs

    def __getRotationTensor(self, mesh):
        if not hasattr(self, 'rotationTensor'):

            rotationTensor = FaceVariable(mesh=mesh, rank=2)

            rotationTensor[:, 0] = self._getNormals(mesh)

            if mesh.dim == 2:
                rotationTensor[:,1] = rotationTensor[:,0].dot((((0, 1), (-1, 0))))
            elif mesh.dim ==3:
                epsilon = 1e-20

                div = numerix.sqrt(1 - rotationTensor[2,0]**2)
                flag = numerix.resize(div > epsilon, (mesh.dim, mesh.numberOfFaces))

                rotationTensor[0, 1] = 1
                rotationTensor[:, 1] = numerix.where(flag,
                                                     rotationTensor[:,0].dot((((0, 1, 0), (-1, 0, 0), (0, 0, 0)))) / div,
                                                     rotationTensor[:, 1])


                rotationTensor[1, 2] = 1
                rotationTensor[:, 2] = numerix.where(flag,
                                                     rotationTensor[:,0] * rotationTensor[2,0] / div,
                                                     rotationTensor[:, 2])
                rotationTensor[2, 2] = -div

            self.rotationTensor = rotationTensor

        return self.rotationTensor

    def __calcAnisotropySource(self, coeff, mesh, var):

        if not hasattr(self, 'anisotropySource'):
            if len(coeff) > 1:
                unconstrainedVar = var + 0
                gradients = unconstrainedVar.grad.harmonicFaceValue.dot(self.__getRotationTensor(mesh))
                from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
                self.anisotropySource = _AddOverFacesVariable(gradients[1:].dot(coeff[1:])) * mesh.cellVolumes

    def _calcGeomCoeff(self, var):

        mesh = var.mesh
        if self.nthCoeff is not None:

            coeff = self.nthCoeff

            shape = numerix.getShape(coeff)

            if isinstance(coeff, FaceVariable):
                rank = coeff.rank
            else:
                rank = len(shape)

            if var.rank == 0:
                anisotropicRank = rank
            elif var.rank == 1:
                anisotropicRank = rank - 2
            else:
                raise IndexError, 'the solution variable has the wrong rank'

            if anisotropicRank == 0 and self._treatMeshAsOrthogonal(mesh):

                if coeff.shape != () and not isinstance(coeff, FaceVariable):
                    coeff = coeff[...,numerix.newaxis]

                tmpBop = (coeff * FaceVariable(mesh=mesh, value=mesh._faceAreas) / mesh._cellDistances)[numerix.newaxis, :]

            else:

                if anisotropicRank == 1 or anisotropicRank == 0:
                    coeff = coeff * numerix.identity(mesh.dim)

                if anisotropicRank > 0:
                    shape = numerix.getShape(coeff)
                    if mesh.dim != shape[0] or mesh.dim != shape[1]:
                        raise IndexError, 'diffusion coefficent tensor is not an appropriate shape for this mesh'

                faceNormals = FaceVariable(mesh=mesh, rank=1, value=mesh.faceNormals)
                rotationTensor = self.__getRotationTensor(mesh)
                rotationTensor[:,0] = rotationTensor[:,0] / mesh._cellDistances

                tmpBop = faceNormals.dot(coeff).dot(rotationTensor) * mesh._faceAreas

            return tmpBop

        else:

            return None

    def _getCoefficientMatrixForTests(self, SparseMatrix, var, coeff):
        """
        This method was introduced because __getCoefficientMatrix is private, but
        the tests in DiffusionTerm need to call it.
        """
        return self.__getCoefficientMatrix(SparseMatrix, var, coeff)

    def __getCoefficientMatrix(self, SparseMatrix, var, coeff):
        mesh = var.mesh

        id1, id2 = mesh._adjacentCellIDs
        interiorFaces = numerix.nonzero(mesh.interiorFaces)[0]

        id1 = numerix.take(id1, interiorFaces)
        id2 = numerix.take(id2, interiorFaces)

        id1 = self._reshapeIDs(var, id1)
        id2 = self._reshapeIDs(var, id2)

##         print 'id1',id1
##         print 'id2',id2

        coefficientMatrix = SparseMatrix(mesh=mesh, bandwidth = mesh._maxFacesPerCell + 1)
        interiorCoeff = numerix.take(coeff, interiorFaces, axis=-1).ravel()
        coefficientMatrix.addAt(interiorCoeff, id1.ravel(), id1.swapaxes(0,1).ravel())
        coefficientMatrix.addAt(-interiorCoeff, id1.ravel(), id2.swapaxes(0,1).ravel())
        coefficientMatrix.addAt(-interiorCoeff, id2.ravel(), id1.swapaxes(0,1).ravel())
        coefficientMatrix.addAt(interiorCoeff, id2.ravel(), id2.swapaxes(0,1).ravel())

##         print 'coefficientMatrix',coefficientMatrix
##         raw_input('stopped')

##         interiorCoeff = numerix.array(coeff)

##         interiorCoeff[...,mesh.exteriorFaces.value] = 0
##         print 'interiorCoeff',interiorCoeff
##         interiorCoeff = numerix.take(interiorCoeff, mesh.cellFaceIDs, axis=-1)
## ##        print interiorCoeff.shape
## ##        print interiorCoeff[:,:,0]
## ##        print interiorCoeff[:,:,1]

##         coefficientMatrix = SparseMatrix(mesh=mesh, bandwidth = mesh._maxFacesPerCell + 1)

## ##        print 'numerix.sum(interiorCoeff, -2)',numerix.sum(interiorCoeff, -2)
## ##        print numerix.sum(interiorCoeff, -2).ravel()
## ##        raw_input('stopped')
## ##        coefficientMatrix.addAtDiagonal(numerix.sum(interiorCoeff, -2).ravel())
## ##        print 'coefficientMatrix',coefficientMatrix

##         del interiorCoeff

##         interiorFaces = mesh.interiorFaceIDs
##         interiorFaceCellIDs = mesh.interiorFaceCellIDs

##         interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)

##         print 'interiorCoeff',interiorCoeff
##         raw_input('stopped')
##         coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[0], interiorFaceCellIDs[1])
##         interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)
##         coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[1], interiorFaceCellIDs[0])

        return coefficientMatrix

    def __bcAdd(self, coefficientMatrix, boundaryB, LL, bb):
        coefficientMatrix += LL
        boundaryB += bb

    def __doBCs(self, SparseMatrix, higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB):
        for boundaryCondition in higherOrderBCs:
            LL, bb = boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffs)
            if 'FIPY_DISPLAY_MATRIX' in os.environ:
                self._viewer.title = r"%s %s" % (boundaryCondition.__class__.__name__, self.__class__.__name__)
                self._viewer.plot(matrix=LL, RHSvector=bb)
                from fipy import raw_input
                raw_input()
            self.__bcAdd(coefficientMatrix, boundaryB, LL, bb)

        return coefficientMatrix, boundaryB

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Test to ensure that a changing coefficient influences the boundary conditions.

        >>> from fipy import *
        >>> m = Grid2D(nx=2, ny=2)
        >>> v = CellVariable(mesh=m)
        >>> c0 = Variable(1.)
        >>> v.constrain(c0, where=m.facesLeft)

        Diffusion will only be in the y-direction

        >>> coeff = Variable([[0. , 0.], [0. , 1.]])
        >>> eq = DiffusionTerm(coeff)
        >>> eq.solve(v, solver=DummySolver())
        >>> print v
        [ 0.  0.  0.  0.]

        Change the coefficient.

        >>> coeff[0, 0] = 1.
        >>> eq.solve(v)
        >>> print v
        [ 1.  1.  1.  1.]

        Change the constraints.

        >>> c0.setValue(2.)
        >>> v.constrain(3., where=m.facesRight)
        >>> print v.faceValue.constraintMask
        [False False False False False False  True False  True  True False  True]
        >>> eq.solve(v)
        >>> print v
        [ 2.25  2.75  2.25  2.75]

        """

        var, L, b = self.__higherOrderbuildMatrix(var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)
        mesh = var.mesh

        if self.order == 2:

            if (not hasattr(self, 'constraintL')) or (not hasattr(self, 'constraintB')):

                normals = FaceVariable(mesh=mesh, rank=1, value=mesh._orientedFaceNormals)

                if len(var.shape) == 1 and len(self.nthCoeff.shape) > 1:
                    nthCoeffFaceGrad = var.faceGrad.dot(self.nthCoeff)
                    normalsNthCoeff =  normals.dot(self.nthCoeff)
                else:

                    if self.nthCoeff.shape != () and not isinstance(self.nthCoeff, FaceVariable):
                        coeff = self.nthCoeff[...,numerix.newaxis]
                    else:
                        coeff = self.nthCoeff

                    nthCoeffFaceGrad = coeff[numerix.newaxis] * var.faceGrad[:,numerix.newaxis]
                    s = (slice(0,None,None),) + (numerix.newaxis,) * (len(coeff.shape) - 1) + (slice(0,None,None),)
                    normalsNthCoeff = coeff[numerix.newaxis] * normals[s]

                self.constraintB = -(var.faceGrad.constraintMask * nthCoeffFaceGrad).divergence * mesh.cellVolumes

                constrainedNormalsDotCoeffOverdAP = var.arithmeticFaceValue.constraintMask * \
                                                    normalsNthCoeff / mesh._cellDistances

                self.constraintB -= (constrainedNormalsDotCoeffOverdAP * var.arithmeticFaceValue).divergence * mesh.cellVolumes

                ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))

                self.constraintL = -constrainedNormalsDotCoeffOverdAP.divergence * mesh.cellVolumes

            ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
            L.addAt(self.constraintL.ravel(), ids.ravel(), ids.swapaxes(0,1).ravel())
            b += numerix.reshape(self.constraintB.ravel(), ids.shape).sum(-2).ravel()

        return (var, L, b)

    def __higherOrderbuildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        mesh = var.mesh

        N = mesh.numberOfCells
        M = mesh._maxFacesPerCell

        if self.order > 2:

            higherOrderBCs, lowerOrderBCs = self.__getBoundaryConditions(boundaryConditions)

            var, lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm._buildMatrix(var = var, SparseMatrix=SparseMatrix,
                                                                                      boundaryConditions = lowerOrderBCs,
                                                                                      dt = dt, transientGeomCoeff=transientGeomCoeff,
                                                                                      diffusionGeomCoeff=diffusionGeomCoeff)
            del lowerOrderBCs

            lowerOrderb = lowerOrderb / mesh.cellVolumes
            volMatrix = SparseMatrix(mesh=var.mesh, bandwidth = 1)

            volMatrix.addAtDiagonal(1. / mesh.cellVolumes)
            lowerOrderL = volMatrix * lowerOrderL
            del volMatrix

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(var)[0]
                minusCoeff = -coeff

                coeff.dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':     minusCoeff,
                    'cell 1 offdiag':  coeff
                    }
                del coeff
                del minusCoeff

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']


            mm = self.__getCoefficientMatrix(SparseMatrix, var, self.coeffDict['cell 1 diag'])
            L, b = self.__doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict,
                               mm, numerix.zeros(len(var.ravel()),'d'))

            del higherOrderBCs
            del mm

            b = L * lowerOrderb + b
            del lowerOrderb

            L = L * lowerOrderL
            del lowerOrderL

        elif self.order == 2:

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(var)
                minusCoeff = -coeff[0]

                coeff[0].dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':    minusCoeff,
                    'cell 1 offdiag':  coeff[0]
                    }

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']

                self.__calcAnisotropySource(coeff, mesh, var)

                del coeff
                del minusCoeff

            higherOrderBCs, lowerOrderBCs = self.__getBoundaryConditions(boundaryConditions)
            del lowerOrderBCs

            L, b = self.__doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict,
                               self.__getCoefficientMatrix(SparseMatrix, var, self.coeffDict['cell 1 diag']), numerix.zeros(len(var.ravel()),'d'))

            if hasattr(self, 'anisotropySource'):
                b -= self.anisotropySource

            del higherOrderBCs


        else:

            L = SparseMatrix(mesh=mesh)
            L.addAtDiagonal(mesh.cellVolumes)
            b = numerix.zeros(len(var.ravel()),'d')

        return (var, L, b)

    def _getDiffusionGeomCoeff(self, var):
        if var is self.var or self.var is None:
            return self._getGeomCoeff(var)
        else:
            return None

    @property
    def _diffusionVars(self):
        return self._vars

    def _treatMeshAsOrthogonal(self, mesh):
        raise NotImplementedError

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
