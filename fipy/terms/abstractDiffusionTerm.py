from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import os

from fipy import input
from fipy.solvers import INDEX_TYPE
from fipy.terms.unaryTerm import _UnaryTerm
from fipy.tools import numerix
from fipy.terms import TermMultiplyError
from fipy.terms import AbstractBaseClassError
from fipy.variables import Variable, FaceVariable, CellVariable

class _AbstractDiffusionTerm(_UnaryTerm):

    def __init__(self, coeff = (1.,), var=None):
        if self.__class__ is _AbstractDiffusionTerm:
            raise AbstractBaseClassError

        self.constraintL = {}
        self.constraintB = {}

        if type(coeff) not in (type(()), type([])):
            coeff = [coeff]

        self.order = len(coeff) * 2


        if len(coeff) > 0:
            self.nthCoeff = coeff[0]

            if not isinstance(self.nthCoeff, Variable):
                self.nthCoeff = Variable(value=self.nthCoeff)

            if isinstance(self.nthCoeff, CellVariable):
                self.nthCoeff = self.nthCoeff.arithmeticFaceValue

        else:
            self.nthCoeff = None

        _UnaryTerm.__init__(self, coeff=coeff, var=var)

        if self.order > 0:
            self.lowerOrderDiffusionTerm = self.__class__(coeff=coeff[1:], var=var)

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

    def _getBoundaryConditions(self, boundaryConditions):
        higherOrderBCs = []
        lowerOrderBCs = []

        for bc in boundaryConditions:
            bcDeriv = bc._getDerivative(self.order - 2)
            if bcDeriv:
                higherOrderBCs.append(bcDeriv)
            else:
                lowerOrderBCs.append(bc)

        return higherOrderBCs, lowerOrderBCs

    def _getRotationTensor(self, mesh):
        if not hasattr(self, 'rotationTensor'):

            rotationTensor = FaceVariable(mesh=mesh, rank=2)

            rotationTensor[:, 0] = self._getNormals(mesh)

            if mesh.dim == 2:
                rotationTensor[:, 1] = rotationTensor[:, 0].dot(((( 0, 1),
                                                                  (-1, 0))))
            elif mesh.dim == 3:
                epsilon = 1e-20

                div = numerix.sqrt(1 - rotationTensor[2, 0]**2)
                flag = numerix.resize(div > epsilon, (mesh.dim, mesh.numberOfFaces))

                rotationTensor[0, 1] = 1
                rotationTensor[:, 1] = numerix.where(flag,
                                                     rotationTensor[:, 0].dot(((( 0, 1, 0),
                                                                                (-1, 0, 0),
                                                                                ( 0, 0, 0)))) / div,
                                                     rotationTensor[:, 1])


                rotationTensor[1, 2] = 1
                rotationTensor[:, 2] = numerix.where(flag,
                                                     rotationTensor[:, 0] * rotationTensor[2, 0] / div,
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

    def _isotropicOrthogonalCoeff(self, coeff, mesh):
        """Geometric coefficient for isotropic diffusion on orthogonal mesh

        Parameters
        ----------
        coeff : array_like
            Diffusion coefficient.
        mesh : ~fipy.meshes.mesh.Mesh
            Geometry and topology.

        Returns
        -------
        ~fipy.variables.faceVariable.FaceVariable
            Contribution to the matrix, consisting of a rank-1 vector at
            each face, composed of the normal contribution.

        Notes
        -----
        .. math::

           \left(\Gamma_0 \frac{A_f}{d_{AP}}\right)^T
        """
        if coeff.shape != () and not isinstance(coeff, FaceVariable):
            # increase dimension of non-scalar coefficient such that
            # it projects to each cell or face it applies to
            coeff = coeff[..., numerix.newaxis]

        return (coeff * FaceVariable(mesh=mesh, value=mesh._faceAreas)
                / mesh._cellDistances)[numerix.newaxis, ...]

    def _anisotropicOrNonorthogonalCoeff(self, coeff, mesh, anisotropicRank):
        """Geometric coefficient for anisotropic diffusion or nonorthogonal mesh

        Parameters
        ----------
        coeff : array_like
            Diffusion coefficient.
        mesh : ~fipy.meshes.mesh.Mesh
            Geometry and topology.
        anisotropyRank : int
            ???

        Returns
        -------
        ~fipy.variables.faceVariable.FaceVariable
            Contribution to the matrix, consisting of a rank-1 vector at
            each face, composed of the normal contribution and the
            tangential contribution(s).

        Notes
        -----
        .. math::

           \hat{n}\cdot\Gamma_0\cdot\mathsf{R} A_f

        where :math:`\mathsf{R}` is the rotation tensor.  The normal
        component of the rotation tensor is scaled by the cell distances.
        """
        if anisotropicRank < 2:
            coeff = coeff * numerix.identity(mesh.dim)

        if anisotropicRank > 0:
            shape = numerix.getShape(coeff)
            if mesh.dim != shape[0] or mesh.dim != shape[1]:
                raise IndexError('diffusion coefficient tensor is not an appropriate shape for this mesh')

        faceNormals = FaceVariable(mesh=mesh, rank=1, value=mesh.faceNormals)
        rotationTensor = self._getRotationTensor(mesh)
        rotationTensor[:, 0] = rotationTensor[:, 0] / mesh._cellDistances

        return faceNormals.dot(coeff).dot(rotationTensor) * mesh._faceAreas

    def _calcGeomCoeff(self, var):
        """Geometric cofficient

        Combination of diffusion coefficient and geometric factor.

        Parameters
        ----------
        var : ~fipy.variables.cellVariable.CellVariable
            Solution variable.

        Returns
        -------
        ~fipy.variables.faceVariable.FaceVariable
            Contribution to the matrix, consisting of a rank-1 vector at
            each face, composed of the normal contribution and (for
            nonorthogonal meshes or anisotropic diffusion) the tangential
            contribution(s).
        """

        mesh = var.mesh
        if self.nthCoeff is not None:

            coeff = self.nthCoeff

            if isinstance(coeff, FaceVariable):
                rank = coeff.rank
            else:
                rank = len(numerix.getShape(coeff))

            if var.rank == 0:
                anisotropicRank = rank
            elif var.rank == 1:
                anisotropicRank = rank - 2
            else:
                raise IndexError('the solution variable has the wrong rank')

            if anisotropicRank == 0 and self._treatMeshAsOrthogonal(mesh):

                return self._isotropicOrthogonalCoeff(coeff, mesh)

            else:

                return self._anisotropicOrNonorthogonalCoeff(coeff, mesh, anisotropicRank)

        else:

            return None

    def _getCoefficientMatrix(self, SparseMatrix, var, coeff):
        mesh = var.mesh

        id1, id2 = mesh._adjacentCellIDs
        interiorFaces = numerix.nonzero(mesh.interiorFaces)[0]

        id1 = numerix.take(id1, interiorFaces)
        id2 = numerix.take(id2, interiorFaces)

        id1 = self._reshapeIDs(var, id1)
        id2 = self._reshapeIDs(var, id2)

        facesPerCell = mesh._facesPerCell[..., mesh._localNonOverlappingCellIDs]
        coefficientMatrix = self._getMatrix(SparseMatrix=SparseMatrix, mesh=mesh, nonZerosPerRow=facesPerCell + 1)
        interiorCoeff = numerix.take(coeff, interiorFaces, axis=-1).ravel()
        coefficientMatrix.addAt(interiorCoeff, id1.ravel(), id1.swapaxes(0, 1).ravel())
        coefficientMatrix.addAt(-interiorCoeff, id1.ravel(), id2.swapaxes(0, 1).ravel())
        coefficientMatrix.addAt(-interiorCoeff, id2.ravel(), id1.swapaxes(0, 1).ravel())
        coefficientMatrix.addAt(interiorCoeff, id2.ravel(), id2.swapaxes(0, 1).ravel())

        return coefficientMatrix

    def _doBCs(self, SparseMatrix, higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB):
        for boundaryCondition in higherOrderBCs:
            LL, bb = boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffs)
            if 'FIPY_DISPLAY_MATRIX' in os.environ:
                self._viewer.title = r"%s %s" % (boundaryCondition.__class__.__name__, self.__class__.__name__)
                self._viewer.plot(matrix=LL, RHSvector=bb)
                from fipy import input
                input()
            coefficientMatrix += LL
            boundaryB += bb

        return coefficientMatrix, boundaryB

    def _constrainValue(self, var):
        """Determine value constraint contributions to matrix and RHS

        Parameters
        ----------
        var : ~fipy.variables.cellVariable.CellVariable
            Constrained solution variable

        Returns
        -------
        L : ~fipy.matrices.sparseMatrix.SparseMatrix
            The NxN sparse matrix contribution.
        b : array_like
            The length-N right-hand-side vector contribution.

        Notes
        -----
        For the variable :math:`\phi`, with its value constrained to
        :math:`\phi\rvert_{\partial\Omega_V} = V` on boundary faces
        :math:`\partial\Omega_V`, determines the matrix contribution

        .. math::

             \begin{align}
                \mathsf{L} &= -\nabla\cdot\left(\frac{\Gamma}{d_{fP}}\hat{n}\right)_{f\in\partial\Omega_V} V_P
                \\
                &\approx -\sum_{f\in\partial\Omega_V}(\frac{\Gamma}{d_{fP}}\hat{n}\cdot\hat{n})_f A_f
             \end{align}

        and the right-hand-side vector contribution

        .. math::

             \begin{align}
                \mathbf{b} &= -\nabla\cdot\left(\frac{\Gamma V}{d_{fP}}\hat{n}\right)_{f\in\partial\Omega_V} V_P
                \\
                &\approx -\sum_{f\in\partial\Omega_V}(\frac{\Gamma V}{d_{fP}}\hat{n}\cdot\hat{n})_f A_f
             \end{align}
        """
        mesh = var.mesh
        normals = FaceVariable(mesh=mesh, rank=1, value=mesh._orientedFaceNormals)

        if len(var.shape) == 1 and len(self.nthCoeff.shape) > 1:
            normalsNthCoeff =  normals.dot(self.nthCoeff)
        else:

            if self.nthCoeff.shape != () and not isinstance(self.nthCoeff, FaceVariable):
                coeff = self.nthCoeff[..., numerix.newaxis]
            else:
                coeff = self.nthCoeff

            s = (slice(0, None, None),) + (numerix.newaxis,) * (len(coeff.shape) - 1) + (slice(0, None, None),)
            normalsNthCoeff = coeff[numerix.newaxis] * normals[s]

        constrainedNormalsDotCoeffOverdAP = var.arithmeticFaceValue.constraintMask * \
                                            normalsNthCoeff / mesh._cellDistances

        L = -constrainedNormalsDotCoeffOverdAP.divergence * mesh.cellVolumes
        b = -(constrainedNormalsDotCoeffOverdAP
              * var.arithmeticFaceValue).divergence * mesh.cellVolumes

        return L, b

    def _constrainGradient(self, var):
        """Determine gradient constraint contributions to matrix and RHS

        Parameters
        ----------
        var : ~fipy.variables.cellVariable.CellVariable
            Constrained solution variable of N cells.

        Returns
        -------
        L : ~fipy.matrices.sparseMatrix.SparseMatrix
            The NxN sparse matrix contribution.
        b : array_like
            The length-N right-hand-side vector contribution

        Notes
        -----
        For the variable :math:`\phi`, with its gradient constrained to
        :math:`\nabla\phi\rvert_{\partial\Omega_G} = \vec{G}` on boundary
        faces :math:`\partial\Omega_G`, determines the matrix contribution

        .. math::

             \begin{align}
                \mathsf{L} &= \mathsf{0}
             \end{align}

        and the right-hand-side vector contribution

        .. math::

             \begin{align}
                 \mathbf{b} &= -\nabla\cdot\left(\Gamma\vec{G}\right)_{f\in\partial\Omega_G} V_P
                 \\
                 &\approx -\sum_{f\in\partial\Omega_G}(\Gamma\vec{G}\cdot\hat{n})_f A_f
             \end{align}
        """
        if len(var.shape) == 1 and len(self.nthCoeff.shape) > 1:
            # var is scalar field and self.nthCoeff is vector (or tensor)
            nthCoeffFaceGrad = var.faceGrad.dot(self.nthCoeff)
        else:
            # var is vector or tensor field or self.nthCoeff is scalar
            if not (self.nthCoeff.shape == () or isinstance(self.nthCoeff, FaceVariable)):
                # self.nthCoeff is not a scalar or a FaceVariable
                coeff = self.nthCoeff[..., numerix.newaxis]
            else:
                # self.nthCoeff is a scalar or a FaceVariable
                coeff = self.nthCoeff

            nthCoeffFaceGrad = coeff[numerix.newaxis] * var.faceGrad[:, numerix.newaxis]

        b = -(var.faceGrad.constraintMask
              * nthCoeffFaceGrad).divergence * mesh.cellVolumes

        return 0, b


    def _calcConstraints(self, var):
        """Determine contributions to matrix and RHS due to constraints on `var`

        Parameters
        ----------
        var : ~fipy.variables.cellVariable.CellVariable
            The constrained variable

        Returns
        -------
        None
        """
        if (var not in self.constraintL) or (var not in self.constraintB):
            LL, bb = self._constrainValue(var)

            self.constraintL[var] = LL
            self.constraintB[var] = bb

            LL, bb = self._constrainGradient(var)

            self.constraintL[var] += LL
            self.constraintB[var] += bb

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        """
        Test to ensure that a changing coefficient influences the boundary conditions.

        >>> from fipy import *
        >>> m = Grid2D(nx=2, ny=2)
        >>> v = CellVariable(mesh=m)
        >>> c0 = Variable(1.)
        >>> v.constrain(c0, where=m.facesLeft)

        Diffusion will only be in the y-direction

        >>> coeff = Variable([[0., 0.], [0., 1.]])
        >>> eq = DiffusionTerm(coeff)
        >>> eq.solve(v, solver=DummySolver())
        >>> print(v)
        [ 0.  0.  0.  0.]

        Change the coefficient.

        >>> coeff[0, 0] = 1.
        >>> eq.solve(v)
        >>> print(v)
        [ 1.  1.  1.  1.]

        Change the constraints.

        >>> c0.setValue(2.)
        >>> v.constrain(3., where=m.facesRight)
        >>> print(v.faceValue.constraintMask)
        [False False False False False False  True False  True  True False  True]
        >>> eq.solve(v)
        >>> print(v)
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
                        coeff = self.nthCoeff[..., numerix.newaxis]
                    else:
                        coeff = self.nthCoeff

                    nthCoeffFaceGrad = coeff[numerix.newaxis] * var.faceGrad[:, numerix.newaxis]
                    s = (slice(0, None, None),) + (numerix.newaxis,) * (len(coeff.shape) - 1) + (slice(0, None, None),)
                    normalsNthCoeff = coeff[numerix.newaxis] * normals[s]

                self.constraintB = -(var.faceGrad.constraintMask * nthCoeffFaceGrad).divergence * mesh.cellVolumes

                constrainedNormalsDotCoeffOverdAP = var.arithmeticFaceValue.constraintMask * \
                                                    normalsNthCoeff / mesh._cellDistances

                self.constraintB -= (constrainedNormalsDotCoeffOverdAP * var.arithmeticFaceValue).divergence * mesh.cellVolumes

                ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))

                self.constraintL = -constrainedNormalsDotCoeffOverdAP.divergence * mesh.cellVolumes

            ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells,
                                                       dtype=INDEX_TYPE))
            L.addAt(self.constraintL.ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
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
            volMatrix = SparseMatrix(mesh=var.mesh, nonZerosPerRow=1)

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
                               mm, numerix.zeros(len(var.ravel()), 'd'))

            del higherOrderBCs
            del mm

            b = numerix.asarray(L * lowerOrderb) + b
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
                               self.__getCoefficientMatrix(SparseMatrix, var, self.coeffDict['cell 1 diag']), numerix.zeros(len(var.ravel()), 'd'))

            if hasattr(self, 'anisotropySource'):
                b -= self.anisotropySource

            del higherOrderBCs


        else:

            L = SparseMatrix(mesh=mesh, nonZerosPerRow=1)
            L.addAtDiagonal(mesh.cellVolumes)
            b = numerix.zeros(len(var.ravel()), 'd')

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


