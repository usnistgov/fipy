#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "cellVariable.py"
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.variables.meshVariable import _MeshVariable
from fipy.tools import numerix

__all__ = ["CellVariable"]

class CellVariable(_MeshVariable):
    """
    Represents the field of values of a variable on a `Mesh`.

    A `CellVariable` can be ``pickled`` to persistent storage (disk) for later use:

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(dx = 1., dy = 1., nx = 10, ny = 10)

    >>> var = CellVariable(mesh = mesh, value = 1., hasOld = 1, name = 'test')
    >>> x, y = mesh.cellCenters
    >>> var.value = (x * y)

    >>> from fipy.tools import dump
    >>> (f, filename) = dump.write(var, extension = '.gz')
    >>> unPickledVar = dump.read(filename, f)

    >>> print var.allclose(unPickledVar, atol = 1e-10, rtol = 1e-10)
    1

    """

    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, unit=None, hasOld=0):
        _MeshVariable.__init__(self, mesh=mesh, name=name, value=value,
                               rank=rank, elementshape=elementshape, unit=unit)

        if hasOld:
            self._old = self.copy()
        else:
            self._old = None

    @property
    def _variableClass(self):
        return CellVariable

    def _OperatorVariableClass(self, baseClass=None):
        """
            >>> from fipy.meshes import Grid1D

            >>> a = CellVariable(mesh=Grid1D(nx=1), value=1)

            >>> c = -a
            >>> b = c.old + 3
            >>> print b
            [2]
            >>> print str(b.getsctype()) == str(numerix.NUMERIX.obj2sctype(numerix.array(1)))
            True

        replacing with the same thing is no problem

            >>> a.value = (3)
            >>> b = c.old + 3
            >>> print b
            [0]

        replacing with multiple copies causes the reference counting problem

            >>> a.value = (3)
            >>> b = (c + c).old + 3
            >>> print b
            [-3]

        the order matters

            >>> b = (c + c).old + 3
            >>> a.value = (2)
            >>> print b
            [-1]
        """
        baseClass = _MeshVariable._OperatorVariableClass(self,
                                                         baseClass=baseClass)

        class _CellOperatorVariable(baseClass):
            @property
            def old(self):
                if self._old is None:
                    oldVar = []
                    for v in self.var:
                        if hasattr(v, "old"):
                            oldVar.append(v.old)
                        else:
                            oldVar.append(v)

                    self._old = self.__class__(op=self.op, var=oldVar,
                                               opShape=self.opShape,
                                               canInline=self.canInline)

                return self._old

        return _CellOperatorVariable

    def copy(self):

        return self._getArithmeticBaseClass()(mesh=self.mesh,
                                              name=self.name + "_old",
                                              value=self.value,
                                              hasOld=False)

    @property
    def _globalNumberOfElements(self):
        return self.mesh.globalNumberOfCells

    @property
    def _globalOverlappingIDs(self):
        return self.mesh._globalOverlappingCellIDs

    @property
    def _localNonOverlappingIDs(self):
        return self.mesh._localNonOverlappingCellIDs

    @property
    def globalValue(self):
        """Concatenate and return values from all processors

        When running on a single processor, the result is identical to
        :attr:`~fipy.variables.variable.Variable.value`.
        """
        return self._getGlobalValue(self.mesh._localNonOverlappingCellIDs,
                                    self.mesh._globalNonOverlappingCellIDs)

    def setValue(self, value, unit = None, where = None):
        _MeshVariable.setValue(self, value=self._globalToLocalValue(value), unit=unit, where=where)

    def __call__(self, points=None, order=0, nearestCellIDs=None):
        r"""
        Interpolates the CellVariable to a set of points using a
        method that has a memory requirement on the order of Ncells by
        Npoints in general, but uses only Ncells when the
        CellVariable's mesh is a UniformGrid object.

        :Parameters:

           - `points`: A point or set of points in the format (X, Y, Z)
           - `order`: The order of interpolation, 0 or 1, default is 0
           - `nearestCellIDs` : Optional argument if user can calculate own
             nearest cell IDs array, shape should be same as points

        Tests

            >>> from fipy import *
            >>> m = Grid2D(nx=3, ny=2)
            >>> v = CellVariable(mesh=m, value=m.cellCenters[0])
            >>> print v(((0., 1.1, 1.2), (0., 1., 1.)))
            [ 0.5  1.5  1.5]
            >>> print v(((0., 1.1, 1.2), (0., 1., 1.)), order=1)
            [ 0.25  1.1   1.2 ]
            >>> m0 = Grid2D(nx=2, ny=2, dx=1., dy=1.)
            >>> m1 = Grid2D(nx=4, ny=4, dx=.5, dy=.5)
            >>> x, y = m0.cellCenters
            >>> v0 = CellVariable(mesh=m0, value=x * y)
            >>> print v0(m1.cellCenters.globalValue)
            [ 0.25  0.25  0.75  0.75  0.25  0.25  0.75  0.75  0.75  0.75  2.25  2.25
              0.75  0.75  2.25  2.25]
            >>> print v0(m1.cellCenters.globalValue, order=1)
            [ 0.125  0.25   0.5    0.625  0.25   0.375  0.875  1.     0.5    0.875
              1.875  2.25   0.625  1.     2.25   2.625]

        """
        if points is not None:

            if nearestCellIDs is None:
                nearestCellIDs = self.mesh._getNearestCellID(points)

            if order == 0:
                return self.globalValue[..., nearestCellIDs]

            elif order == 1:
                ##cellID = self.mesh._getNearestCellID(points)
##                return self[...,self.mesh._getNearestCellID(points)] + numerix.dot(points - self.mesh.cellCenters[...,cellID], self.grad[...,cellID])
                return (self.globalValue[..., nearestCellIDs]
                        + numerix.dot(points - self.mesh.cellCenters.globalValue[...,nearestCellIDs],
                                      self.grad.globalValue[...,nearestCellIDs]))

            else:
                raise ValueError, 'order should be either 0 or 1'

        else:
            return _MeshVariable.__call__(self)

    @property
    def cellVolumeAverage(self):
        r"""
        Return the cell-volume-weighted average of the `CellVariable`:

        .. math::

           <\phi>_\text{vol}
           = \frac{\sum_\text{cells} \phi_\text{cell} V_\text{cell}}
               {\sum_\text{cells} V_\text{cell}}

        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(nx = 3, ny = 1, dx = .5, dy = .1)
        >>> var = CellVariable(value = (1, 2, 6), mesh = mesh)
        >>> print var.cellVolumeAverage
        3.0
        """

        if not hasattr(self, 'volumeAverage'):
            from fipy.variables.cellVolumeAverageVariable import _CellVolumeAverageVariable
            self.volumeAverage = _CellVolumeAverageVariable(self)

        return self.volumeAverage

    @property
    def grad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `CellVariable` (first-order
        gradient).
        """
        return self.gaussGrad

    @property
    def gaussGrad(self):
        r"""
        Return :math:`\frac{1}{V_P} \sum_f \vec{n} \phi_f A_f`
        as a rank-1 `CellVariable` (first-order gradient).

        """
        if not hasattr(self, '_gaussGrad'):
            from fipy.variables.gaussCellGradVariable import _GaussCellGradVariable
            self._gaussGrad = _GaussCellGradVariable(var = self, name = "%s_gauss_grad" % self.name)

        return self._gaussGrad

    @property
    def leastSquaresGrad(self):
        r"""
        Return :math:`\nabla \phi`, which is determined by solving for :math:`\nabla \phi`
        in the following matrix equation,

        .. math::

           \nabla \phi \cdot \sum_f d_{AP}^2 \vec{n}_{AP} \otimes \vec{n}_{AP} =
           \sum_f d_{AP}^2 \left( \vec{n} \cdot \nabla \phi \right)_{AP}

        The matrix equation is derived by minimizing
        the following least squares sum,

        .. math::

           F \left( \phi_x, \phi_y \right) = \sqrt{\sum_f \left( d_{AP}
           \vec{n}_{AP} \cdot \nabla \phi - d_{AP} \left( \vec{n}_{AP} \cdot
           \nabla \phi \right)_{AP} \right)^2 }

        Tests

        >>> from fipy import Grid2D
        >>> m = Grid2D(nx=2, ny=2, dx=0.1, dy=2.0)
        >>> print numerix.allclose(CellVariable(mesh=m, value=(0,1,3,6)).leastSquaresGrad.globalValue, \
        ...                                     [[8.0, 8.0, 24.0, 24.0],
        ...                                      [1.2, 2.0, 1.2, 2.0]])
        True

        >>> from fipy import Grid1D
        >>> print numerix.allclose(CellVariable(mesh=Grid1D(dx=(2.0, 1.0, 0.5)),
        ...                                     value=(0, 1, 2)).leastSquaresGrad.globalValue, [[0.461538461538, 0.8, 1.2]])
        True
        """

        if not hasattr(self, '_leastSquaresGrad'):
            from fipy.variables.leastSquaresCellGradVariable import _LeastSquaresCellGradVariable
            self._leastSquaresGrad = _LeastSquaresCellGradVariable(var = self,
                    name = "%s_least_squares_grad" % self.name)

        return self._leastSquaresGrad

    @property
    def arithmeticFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the arithmetic interpolation
        of the adjacent cells:

        .. math::

           \phi_f = (\phi_1 - \phi_2) \frac{d_{f2}}{d_{12}} + \phi_2

        >>> from fipy.meshes import Grid1D
        >>> from fipy import numerix
        >>> mesh = Grid1D(dx = (1., 1.))
        >>> L = 1
        >>> R = 2
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.arithmeticFaceValue[mesh.interiorFaces.value]
        >>> answer = (R - L) * (0.5 / 1.) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (2., 4.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.arithmeticFaceValue[mesh.interiorFaces.value]
        >>> answer = (R - L) * (1.0 / 3.0) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (10., 100.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.arithmeticFaceValue[mesh.interiorFaces.value]
        >>> answer = (R - L) * (5.0 / 55.0) + L
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        """
        if not hasattr(self, '_arithmeticFaceValue'):
            from fipy.variables.arithmeticCellToFaceVariable import _ArithmeticCellToFaceVariable
            self._arithmeticFaceValue = _ArithmeticCellToFaceVariable(self)

        return self._arithmeticFaceValue

    faceValue = arithmeticFaceValue

    @property
    def minmodFaceValue(self):
        r"""
        Returns a `FaceVariable` with a value that is the minimum of
        the absolute values of the adjacent cells. If the values are
        of opposite sign then the result is zero:

        .. math::

           \phi_f = \begin{cases}
                          \phi_1& \text{when $|\phi_1| \le |\phi_2|$},\\
                          \phi_2& \text{when $|\phi_2| < |\phi_1|$},\\
                          0 & \text{when $\phi1 \phi2 < 0$}
                    \end{cases}

        >>> from fipy import *
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(1, 2)).minmodFaceValue
        [1 1 2]
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(-1, -2)).minmodFaceValue
        [-1 -1 -2]
        >>> print CellVariable(mesh=Grid1D(nx=2), value=(-1, 2)).minmodFaceValue
        [-1  0  2]
        """
        if not hasattr(self, '_minmodFaceValue'):
            from fipy.variables.minmodCellToFaceVariable import _MinmodCellToFaceVariable
            self._minmodFaceValue = _MinmodCellToFaceVariable(self)

        return self._minmodFaceValue

    @property
    def harmonicFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the harmonic interpolation
        of the adjacent cells:

        .. math::

           \phi_f = \frac{\phi_1 \phi_2}{(\phi_2 - \phi_1) \frac{d_{f2}}{d_{12}} + \phi_1}

        >>> from fipy.meshes import Grid1D
        >>> from fipy import numerix
        >>> mesh = Grid1D(dx = (1., 1.))
        >>> L = 1
        >>> R = 2
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.harmonicFaceValue[mesh.interiorFaces.value]
        >>> answer = L * R / ((R - L) * (0.5 / 1.) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (2., 4.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.harmonicFaceValue[mesh.interiorFaces.value]
        >>> answer = L * R / ((R - L) * (1.0 / 3.0) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True

        >>> mesh = Grid1D(dx = (10., 100.))
        >>> var = CellVariable(mesh = mesh, value = (L, R))
        >>> faceValue = var.harmonicFaceValue[mesh.interiorFaces.value]
        >>> answer = L * R / ((R - L) * (5.0 / 55.0) + L)
        >>> print numerix.allclose(faceValue, answer, atol = 1e-10, rtol = 1e-10)
        True
        """
        if not hasattr(self, '_harmonicFaceValue'):
            from fipy.variables.harmonicCellToFaceVariable import _HarmonicCellToFaceVariable
            self._harmonicFaceValue = _HarmonicCellToFaceVariable(self)

        return self._harmonicFaceValue

    @property
    def faceGrad(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` using differencing
        for the normal direction(second-order gradient).
        """
        if not hasattr(self, '_faceGrad'):
            from fipy.variables.faceGradVariable import _FaceGradVariable
            self._faceGrad = _FaceGradVariable(self)

        return self._faceGrad

    @property
    def faceGradAverage(self):
        r"""
        Return :math:`\nabla \phi` as a rank-1 `FaceVariable` using averaging
        for the normal direction(second-order gradient)
        """
        return self.grad.arithmeticFaceValue

    @property
    def old(self):
        """
        Return the values of the `CellVariable` from the previous
        solution sweep.

        Combinations of `CellVariable's` should also return old
        values.

        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx = 2)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> var1 = CellVariable(mesh = mesh, value = (2, 3), hasOld = 1)
        >>> var2 = CellVariable(mesh = mesh, value = (3, 4))
        >>> v = var1 * var2
        >>> print v
        [ 6 12]
        >>> var1.value = ((3,2))
        >>> print v
        [9 8]
        >>> print v.old
        [ 6 12]

        The following small test is to correct for a bug when the
        operator does not just use variables.

        >>> v1 = var1 * 3
        >>> print v1
        [9 6]
        >>> print v1.old
        [6 9]
        """
        if self._old is None:
            return self
        else:
            return self._old
##             import weakref
##          return weakref.proxy(self._old)

    def updateOld(self):
        """
        Set the values of the previous solution sweep to the current
        values.

        >>> from fipy import *
        >>> v = CellVariable(mesh=Grid1D(), hasOld=False)
        >>> v.updateOld()
        Traceback (most recent call last):
           ...
        AssertionError: The updateOld method requires the CellVariable to have an old value. Set hasOld to True when instantiating the CellVariable.

        """
        if self._old is None:
            raise AssertionError, 'The updateOld method requires the CellVariable to have an old value. Set hasOld to True when instantiating the CellVariable.'
        else:
            self._old.value = self.value.copy()

    def _resetToOld(self):
        if self._old is not None:
            self.value = (self._old.value)

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this variable type, given a particular mesh.
        """
        return (mesh.numberOfCells,)

    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return CellVariable

        return _MeshVariable._getArithmeticBaseClass(self, other)

##pickling

    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the
        `CellVariable` to persistent storage.
        """
        return {
            'mesh' : self.mesh,
            'name' : self.name,
            'value' : self.globalValue,
            'unit' : self.unit,
            'old' : self._old
        }

    def __setstate__(self, dict):
        """
        Used internally to create a new `CellVariable` from ``pickled``
        persistent storage.
        """

        import sys
        self._refcount = sys.getrefcount(self)

        hasOld = 0
        if dict['old'] is not None:
            hasOld = 1

        self.__init__(mesh=dict['mesh'], name=dict['name'], value=dict['value'], unit=dict['unit'], hasOld=hasOld)
##         self.__init__(hasOld=hasOld, **dict)
        if self._old is not None:
            self._old.value = (dict['old'].value)

    def constrain(self, value, where=None):
        r"""
        Constrains the `CellVariable` to `value` at a location specified by `where`.

            >>> from fipy import *
            >>> m = Grid1D(nx=3)
            >>> v = CellVariable(mesh=m, value=m.cellCenters[0])
            >>> v.constrain(0., where=m.facesLeft)
            >>> v.faceGrad.constrain([1.], where=m.facesRight)
            >>> print v.faceGrad
            [[ 1.  1.  1.  1.]]
            >>> print v.faceValue
            [ 0.   1.   2.   2.5]

        Changing the constraint changes the dependencies

            >>> v.constrain(1., where=m.facesLeft)
            >>> print v.faceGrad
            [[-1.  1.  1.  1.]]
            >>> print v.faceValue
            [ 1.   1.   2.   2.5]

        Constraints can be `Variable`

            >>> c = Variable(0.)
            >>> v.constrain(c, where=m.facesLeft)
            >>> print v.faceGrad
            [[ 1.  1.  1.  1.]]
            >>> print v.faceValue
            [ 0.   1.   2.   2.5]
            >>> c.value = 1.
            >>> print v.faceGrad
            [[-1.  1.  1.  1.]]
            >>> print v.faceValue
            [ 1.   1.   2.   2.5]

        Constraints can have a `Variable` mask.

            >>> v = CellVariable(mesh=m)
            >>> mask = FaceVariable(mesh=m, value=m.facesLeft)
            >>> v.constrain(1., where=mask)
            >>> print v.faceValue
            [ 1.  0.  0.  0.]
            >>> mask[:] = mask | m.facesRight
            >>> print v.faceValue
            [ 1.  0.  0.  1.]

        """
        from fipy.boundaryConditions.constraint import Constraint
        if not isinstance(value, Constraint):
            value = Constraint(value=value, where=where)

        if numerix.shape(value.where)[-1] == self.mesh.numberOfFaces:

            if not hasattr(self, 'faceConstraints'):
                self.faceConstraints = []
            self.faceConstraints.append(value)
            self._requires(value.value)
            # self._requires(value.where) ???
            self._markStale()
        else:
##            _MeshVariable.constrain(value, where)
            super(CellVariable, self).constrain(value, where)

    def release(self, constraint):
        """Remove `constraint` from `self`

        >>> from fipy import *
        >>> m = Grid1D(nx=3)
        >>> v = CellVariable(mesh=m, value=m.cellCenters[0])
        >>> c = Constraint(0., where=m.facesLeft)
        >>> v.constrain(c)
        >>> print v.faceValue
        [ 0.   1.   2.   2.5]
        >>> v.release(constraint=c)
        >>> print v.faceValue
        [ 0.5  1.   2.   2.5]
        """
        try:
            _MeshVariable.release(self, constraint=constraint)
        except ValueError:
            self.faceConstraints.remove(constraint)

    def _test(self):
        """
        Tests

        >>> from fipy import *
        >>> m = Grid1D(nx=6)
        >>> q = CellVariable(mesh=m, elementshape=(2,))
        >>> print q.faceGrad.globalValue.shape
        (1, 2, 7)
        >>> from fipy import *
        >>> m = Grid2D(nx=3, ny=3)
        >>> x, y = m.cellCenters
        >>> v = CellVariable(mesh=m, elementshape=(3,))
        >>> v[0] = x
        >>> v[1] = y
        >>> v[2] = x**2
        >>> print v.faceGrad.allclose([[[ 0.5, 1.,  0.5, 0.5, 1.,  0.5, 0.5, 1.,  0.5, 0.5, 1.,  0.5, 0.,  1.,  1.,
        ...                               0.,  0.,  1.,  1.,  0.,  0.,  1.,  1.,  0. ],
        ...                             [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        ...                               0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. ],
        ...                             [ 1.,  3.,  2.,  1.,  3.,  2.,  1.,  3.,  2.,  1.,  3.,  2.,  0.,  2.,  4.,
        ...                               0.,  0.,  2.,  4.,  0.,  0.,  2.,  4.,  0. ]],
        ...                            [[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        ...                               0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. ],
        ...                             [ 0.,  0.,  0.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.,  0.5, 0.5,
        ...                               0.5, 0.5, 1.,  1.,  1.,  1.,  0.5, 0.5, 0.5, 0.5],
        ...                             [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        ...                               0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0. ]]])
        True
        >>> print v.grad
        [[[ 0.5  1.   0.5  0.5  1.   0.5  0.5  1.   0.5]
          [ 0.   0.   0.   0.   0.   0.   0.   0.   0. ]
          [ 1.   3.   2.   1.   3.   2.   1.   3.   2. ]]
        <BLANKLINE>
         [[ 0.   0.   0.   0.   0.   0.   0.   0.   0. ]
          [ 0.5  0.5  0.5  1.   1.   1.   0.5  0.5  0.5]
          [ 0.   0.   0.   0.   0.   0.   0.   0.   0. ]]]

        """

class _ReMeshedCellVariable(CellVariable):
    def __init__(self, oldVar, newMesh):
        newValues = oldVar.getValue(points = newMesh.cellCenters)
        CellVariable.__init__(self, newMesh, name = oldVar.name, value = newValues, unit = oldVar.unit)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
