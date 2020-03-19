from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.terms.faceTerm import FaceTerm
from fipy.variables.meshVariable import _MeshVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.cellVariable import CellVariable
from fipy.terms import AbstractBaseClassError
from fipy.terms import VectorCoeffError
from fipy.tools import numerix

class _AbstractConvectionTerm(FaceTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=1.0, var=None):
        """
        Create a `_AbstractConvectionTerm` object.

            >>> from fipy import *
            >>> m = Grid1D(nx = 2)
            >>> cv = CellVariable(mesh = m)
            >>> fv = FaceVariable(mesh = m)
            >>> vcv = CellVariable(mesh=m, rank=1)
            >>> vfv = FaceVariable(mesh=m, rank=1)
            >>> __ConvectionTerm(coeff = cv) # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = fv) # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> __ConvectionTerm(coeff = vcv)
            __ConvectionTerm(coeff=_ArithmeticCellToFaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = vfv)
            __ConvectionTerm(coeff=FaceVariable(value=array([[ 0.,  0.,  0.]]), mesh=UniformGrid1D(dx=1.0, nx=2)))
            >>> __ConvectionTerm(coeff = (1,))
            __ConvectionTerm(coeff=(1,))
            >>> ExplicitUpwindConvectionTerm(coeff = (0,)).solve(var=cv, solver=DummySolver()) # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            TransientTermError: The equation requires a TransientTerm with explicit convection.
            >>> (TransientTerm(0.) - ExplicitUpwindConvectionTerm(coeff = (0,))).solve(var=cv, solver=DummySolver(), dt=1.)

            >>> (TransientTerm() - ExplicitUpwindConvectionTerm(coeff = 1)).solve(var=cv, solver=DummySolver(), dt=1.) # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            VectorCoeffError: The coefficient must be a vector value.
            >>> m2 = Grid2D(nx=2, ny=1)
            >>> cv2 = CellVariable(mesh=m2)
            >>> vcv2 = CellVariable(mesh=m2, rank=1)
            >>> vfv2 = FaceVariable(mesh=m2, rank=1)
            >>> __ConvectionTerm(coeff=vcv2)
            __ConvectionTerm(coeff=_ArithmeticCellToFaceVariable(value=array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
                   [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]]), mesh=UniformGrid2D(dx=1.0, nx=2, dy=1.0, ny=1)))
            >>> __ConvectionTerm(coeff=vfv2)
            __ConvectionTerm(coeff=FaceVariable(value=array([[ 0.,  0.,  0.,  0.,  0.,  0.,  0.],
                   [ 0.,  0.,  0.,  0.,  0.,  0.,  0.]]), mesh=UniformGrid2D(dx=1.0, nx=2, dy=1.0, ny=1)))
            >>> (TransientTerm() - ExplicitUpwindConvectionTerm(coeff = ((0,), (0,)))).solve(var=cv2, solver=DummySolver(), dt=1.)
            >>> (TransientTerm() - ExplicitUpwindConvectionTerm(coeff = (0, 0))).solve(var=cv2, solver=DummySolver(), dt=1.)


        Parameters
        ----------
        coeff : The `Term`'s coefficient value.
        """
        if self.__class__ is _AbstractConvectionTerm:
            raise AbstractBaseClassError

        self.stencil = None

        if isinstance(coeff, _MeshVariable) and coeff.rank < 1:
            raise VectorCoeffError

        if isinstance(coeff, CellVariable):
            coeff = coeff.arithmeticFaceValue

        FaceTerm.__init__(self, coeff=coeff, var=var)

    def _calcGeomCoeff(self, var):
        mesh = var.mesh

        if not isinstance(self.coeff, FaceVariable):
            shape = numerix.array(self.coeff).shape

            if shape != () and shape != (1,) and shape[-1] == 1:
                shape = shape[:-1]

            self.coeff = FaceVariable(mesh=mesh, elementshape=shape, value=self.coeff)

        projectedCoefficients = self.coeff * mesh._orientedAreaProjections

        return projectedCoefficients.sum(0)

    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        r"""
        Testing that the sign of the equation is taken into account
        when evaluation upwind direction.

        >>> from fipy import Grid1D, CellVariable
        >>> from fipy import TransientTerm, ConvectionTerm, ImplicitSourceTerm

        >>> m = Grid1D(nx=3, dx=0.5)
        >>> v = CellVariable(mesh=m)
        >>> v.constrain(1, m.facesLeft)
        >>> v.constrain(0, m.facesRight)

        >>> (-TransientTerm(coeff=1 / m.x) ==
        ...  ConvectionTerm(coeff=[[1]])
        ...  + ImplicitSourceTerm(coeff=m.x)).solve(v, dt=1.)

        >>> v0 = v.copy()
        >>> v[:] = 0
        >>> (TransientTerm(coeff=1 / m.x) ==
        ...  - ConvectionTerm(coeff=[[1]])
        ...  - ImplicitSourceTerm(coeff=m.x)).solve(v, dt=1.)

        >>> print(numerix.allclose(v, v0))
        True

        """

        if self.stencil is None:

            geomCoeff = self._getGeomCoeff(var)
            large = 1e+20
            pecletLarge = large - (geomCoeff < 0) * (2 * large)
            if numerix.all(self._getDiagonalSign(transientGeomCoeff, diffusionGeomCoeff) < 0):
                pecletLarge = -pecletLarge

            if diffusionGeomCoeff is None or diffusionGeomCoeff[0] is None:
                peclet = pecletLarge
            else:
                diffCoeff = diffusionGeomCoeff[0].numericValue
                diffCoeff = diffCoeff - (diffCoeff == 0) * geomCoeff / pecletLarge
                peclet = -geomCoeff / diffCoeff

            alpha = self._alpha(peclet)

            self.stencil = {'implicit' : {'cell 1 diag'    : alpha,
                                          'cell 1 offdiag' : (1-alpha),
                                          'cell 2 diag'    : -(1-alpha),
                                          'cell 2 offdiag' : -alpha}}

        return self.stencil

    def _checkVar(self, var):
        FaceTerm._checkVar(self, var)
        if not (isinstance(self.coeff, FaceVariable) and self.coeff.rank == 1):
            coeffShape = numerix.getShape(self.coeff)
            if (coeffShape is ()) or (coeffShape[0] != var.mesh.dim):
                raise VectorCoeffError

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):

        var, L, b = FaceTerm._buildMatrix(self, var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

##        if var.rank != 1:

        mesh = var.mesh

        if (not hasattr(self, 'constraintL')) or (not hasattr(self, 'constraintB')):

            weight = self._getWeight(var, transientGeomCoeff, diffusionGeomCoeff)

            if 'implicit' in weight:
                alpha = weight['implicit']['cell 1 diag']
            else:
                alpha = 0.0

            alpha_constraint = numerix.where(var.arithmeticFaceValue.constraintMask, 0.0, alpha)
            alpha_constraint = numerix.where(var.faceGrad.constraintMask, 1.0, alpha)

            def divergence(face_value):
                return (
                    face_value * \
                    (var.faceGrad.constraintMask | var.arithmeticFaceValue.constraintMask) * \
                    self.coeff * mesh.exteriorFaces
                ).divergence * mesh.cellVolumes

            self.constraintL = divergence(alpha_constraint)
            dvar = (var.faceGrad * mesh._cellDistances * mesh.faceNormals).sum(axis=0)
            self.constraintB = divergence(
                (alpha_constraint - 1) * var.arithmeticFaceValue + (alpha - 1)  * dvar * var.faceGrad.constraintMask
            )


        ids = self._reshapeIDs(var, numerix.arange(mesh.numberOfCells))
        L.addAt(numerix.array(self.constraintL).ravel(), ids.ravel(), ids.swapaxes(0, 1).ravel())
        b += numerix.reshape(self.constraintB.value, ids.shape).sum(0).ravel()

        return (var, L, b)


    def _test(self):
        """Test cases for convection with constraints.

        The following tests both a Dirichlet and a Neumann type
        boundary condition with convection using constraints. This
        tests that the Neumann boundary condition converges without
        multiple iterations of the solve step. The Dirichlet boundary
        condition is second order accurate. The Neumann is only first
        order accurate.

        >>> import fipy as fp
        >>> import numpy as np

        >>> def solve(nx, constrain_right):
        ...     var = fp.CellVariable(fp.Grid1D(nx=nx, dx=1. / nx))
        ...     var.constrain(1.0, var.mesh.facesLeft)
        ...     constrain_right(var)
        ...     eqn = fp.CentralDifferenceConvectionTerm((1.0,)) - fp.DiffusionTerm()
        ...     eqn.sweep(var, solver=fp.LinearLUSolver())
        ...     return var

        >>> def error(nx, constrain_right, expected):
        ...     var = solve(nx, constrain_right)
        ...     return np.sqrt((abs(var.value - expected(var))**2 * var.mesh.dx).sum())

        >>> def error_dirichlet(nx):
        ...     return error(
        ...         nx,
        ...         lambda v: v.constrain(0.0, v.mesh.facesRight),
        ...         lambda v: (np.exp(v.mesh.x) - np.exp(1.)) / (1 - np.exp(1.))
        ...     )

        >>> def error_neumann(nx, dphi=1.0):
        ...     return error(
        ...         nx,
        ...         lambda v: v.faceGrad.constrain([dphi], v.mesh.facesRight),
        ...         lambda v: 1 + dphi * (np.exp(v.mesh.x) - 1) / np.exp(1.)
        ...     )

        >>> def slope(ferror, nx0, nx1):
        ...     return np.log(ferror(nx1) / ferror(nx0)) / np.log(nx0 / nx1)

        >>> assert np.allclose(slope(error_dirichlet, 64, 128), 2.0, atol=0.02)
        >>> assert np.allclose(slope(error_neumann, 64, 128), 1.0, atol=0.002)

        """

class __ConvectionTerm(_AbstractConvectionTerm):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
