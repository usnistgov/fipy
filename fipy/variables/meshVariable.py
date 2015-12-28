#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 #
 # FILE: "meshVariable.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.variable import Variable
from fipy.variables.constant import _Constant
from fipy.tools import numerix

class _MeshVariable(Variable):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.

    Abstract base class for a `Variable` that is defined on a mesh
    """
    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None,
                 unit=None, cached=1):
        """
        :Parameters:
          - `mesh`: the mesh that defines the geometry of this `Variable`
          - `name`: the user-readable name of the `Variable`
          - `value`: the initial value
          - `rank`: the rank (number of dimensions) of each element of this
            `Variable`. Default: 0
          - `elementshape`: the shape of each element of this variable
             Default: `rank * (mesh.dim,)`
          - `unit`: the physical units of the `Variable`
        """
        if isinstance(value, (list, tuple)):
            value = numerix.array(value)

        if isinstance(value, _MeshVariable):
            if mesh is None:
                mesh = value.mesh
            elif mesh != value.mesh:
                raise ValueError, "The new 'Variable' must use the same mesh as the supplied value"

        self.mesh = mesh
        value = self._globalToLocalValue(value)

        if value is None:
            array = None
        elif not isinstance(value, _Constant) and isinstance(value, Variable):
            name = name or value.name
            unit = None
            if isinstance(value, _MeshVariable):
                if not isinstance(value, self._variableClass):
                    raise TypeError, "A '%s' cannot be cast to a '%s'" % (value._variableClass.__name__,
                                                                          self._variableClass.__name__)
                if elementshape is not None and elementshape != value.shape[:-1]:
                    raise ValueError, "'elementshape' != shape of elements of 'value'"

                if rank is not None and rank != value.rank:
                    raise ValueError, "'rank' != rank of 'value'"

                elementshape = value.shape[:-1]
                array = None

#             value = value._copyValue()

        if elementshape is None:
            valueShape = numerix.getShape(value)
            if valueShape != () and valueShape[-1] == self._getShapeFromMesh(mesh)[-1]:
                if elementshape is not None and elementshape != valueShape[:-1]:
                    raise ValueError, "'elementshape' != shape of elements of 'value'"

                if rank is not None and rank != len(valueShape[:-1]):
                    raise ValueError, "'rank' != rank of 'value'"
                elementshape = valueShape[:-1]
            elif rank is None and elementshape is None:
                elementshape = valueShape

        if rank is None:
            if elementshape is None:
                elementshape = ()
        elif elementshape is None:
            elementshape = rank * (mesh.dim,)
        elif len(elementshape) != rank:
            raise ValueError, 'len(elementshape) != rank'

        self.elementshape = elementshape

        if not "array" in locals():
            if numerix._isPhysical(value):
                dtype = numerix.obj2sctype(value.value)
            else:
                dtype = numerix.obj2sctype(value)
            #print "meshvariable elshape: ",self.elementshape
            #print "meshvariable _getShapeFromMesh: ",self._getShapeFromMesh(mesh)
            array = numerix.zeros(self.elementshape
                                  + self._getShapeFromMesh(mesh),
                                  dtype)
            if numerix._broadcastShape(array.shape, numerix.shape(value)) is None:
                if not isinstance(value, Variable):
                    value = _Constant(value)
                value = value[..., numerix.newaxis]

        Variable.__init__(self, name=name, value=value, unit=unit,
                          array=array, cached=cached)

    def _globalToLocalValue(self, value):
        if value is not None:
            if not isinstance(value, Variable):
                value = _Constant(value)
            valueShape = value.shape
            if valueShape is not () and valueShape[-1] == self._globalNumberOfElements:
                if valueShape[-1] != 0:
                    # workaround for NumPy:ticket:1171
                    value = value[..., self._globalOverlappingIDs]

            value = value.value
        return value

    def _getGlobalValue(self, localIDs, globalIDs):
        localValue = self.value
        if self.mesh.communicator.Nproc > 1:
            if localValue.shape[-1] != 0:
                localValue = localValue[..., localIDs]
            globalIDs = numerix.concatenate(self.mesh.communicator.allgather(globalIDs))

            globalValue = numerix.empty(localValue.shape[:-1] + (max(globalIDs) + 1,),
                                        dtype=numerix.obj2sctype(localValue))
            globalValue[..., globalIDs] = numerix.concatenate(self.mesh.communicator.allgather(localValue), axis=-1)

            return globalValue
        else:
            return localValue

    def __str__(self):
        return str(self.globalValue)

    def __repr__(self):
        if hasattr(self, 'name') and len(self.name) > 0:
            return self.name
        else:
            s = self.__class__.__name__ + '('
            s += 'value=' + repr(self.globalValue)
            s += ')'
            if len(self.name) == 0:
                s = s[:-1] + ', mesh=' + repr(self.mesh) + s[-1]
            return s

    @property
    def constraintMask(self):
        r"""
        Test that `constraintMask` returns a Variable that updates itself whenver the constraints change.

        >>> from fipy import *

        >>> m = Grid2D(nx=2, ny=2)
        >>> x, y = m.cellCenters
        >>> v0 = CellVariable(mesh=m)
        >>> v0.constrain(1., where=m.facesLeft)
        >>> print v0.faceValue.constraintMask
        [False False False False False False  True False False  True False False]
        >>> print v0.faceValue
        [ 0.  0.  0.  0.  0.  0.  1.  0.  0.  1.  0.  0.]
        >>> v0.constrain(3., where=m.facesRight)
        >>> print v0.faceValue.constraintMask
        [False False False False False False  True False  True  True False  True]
        >>> print v0.faceValue
        [ 0.  0.  0.  0.  0.  0.  1.  0.  3.  1.  0.  3.]
        >>> v1 = CellVariable(mesh=m)
        >>> v1.constrain(1., where=(x < 1) & (y < 1))
        >>> print v1.constraintMask
        [ True False False False]
        >>> print v1
        [ 1.  0.  0.  0.]
        >>> v1.constrain(3., where=(x > 1) & (y > 1))
        >>> print v1.constraintMask
        [ True False False  True]
        >>> print v1
        [ 1.  0.  0.  3.]

        """
        if not hasattr(self, '_constraintMask'):
            from fipy.variables.constraintMask import _ConstraintMask
            self._constraintMask = _ConstraintMask(self)
        return self._constraintMask

    def constrain(self, value, where=None):
        super(_MeshVariable, self).constrain(value, where=where)
        if hasattr(self, '_constraintMask'):
            self._constraintMask._requires(self._constraints[-1].where)

    def _getShapeFromMesh(mesh):
        """
        Return the shape of this `MeshVariable` type, given a particular mesh.
        Return `None` if unknown or independent of the mesh.
        """
        return None

    _getShapeFromMesh = staticmethod(_getShapeFromMesh)

    @property
    def shape(self):
        """
        >>> from fipy.meshes import Grid2D
        >>> from fipy.variables.cellVariable import CellVariable
        >>> mesh = Grid2D(nx=2, ny=3)
        >>> var = CellVariable(mesh=mesh)
        >>> print numerix.allequal(var.shape, (6,)) # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(var.arithmeticFaceValue.shape, (17,)) # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(var.grad.shape, (2, 6)) # doctest: +PROCESSOR_0
        True
        >>> print numerix.allequal(var.faceGrad.shape, (2, 17)) # doctest: +PROCESSOR_0
        True

        Have to account for zero length arrays

        >>> from fipy import Grid1D
        >>> m = Grid1D(nx=0)
        >>> v = CellVariable(mesh=m, elementshape=(2,))
        >>> (v * 1).shape
        (2, 0)

        """
        return (Variable._getShape(self)
                or (self.elementshape + self._getShapeFromMesh(self.mesh))
                or ())

    def _dot(a, b, index):
        """
        Workhorse method to calculate the scalar product

        .. math::

           \mathsf{a} \cdot \mathsf{b}

        for all but the last index of `a` and `b`. Both `a` and `b` can be of
        arbitrary rank, but at this point, both must be appropriately broadcast
        `array` objects.
        """
        rankA = len(a.shape) - 1
        rankB = len(b.shape) - 1
        if rankA <= 0 or rankB <= 0:
            # if either a or b is scalar, then just do multiplication
            return a[index] * b
        else:
            return numerix.sum(a[index] * b, axis=rankA - 1)
    _dot = staticmethod(_dot)

    def __dot(A, B, operatorClass):
        """
        Workhorse method to return a `_BinaryOperatorVariable` that will
        dynamically perform the mesh-element--by--mesh-element (cell-by-cell,
        face-by-face, etc.) scalar product

        .. math::

           \mathsf{A} \cdot \mathsf{B}

        Both `A` and `B` can be of arbitrary rank, but at this point, both must
        be appropriately broadcast `_MeshVariable` objects.

        Test for inline bug

        >>> from fipy import *
        >>> m = Grid2D(nx=2, ny=2)
        >>> v = FaceVariable(mesh=m, rank=1, value=m._orientedFaceNormals)

        >>> print len(v.dot(1.).shape)
        2
        >>> print v.dot(1.).globalValue.shape
        (2, 12)
        >>> tmp = m._cellDistances * v.dot(1.)
        >>> print tmp.globalValue.shape
        (2, 12)

        The value shouldn't change shape the second time it's
        evaluated. The second time is inline and the inline code does
        not have the correct shape.

        >>> print tmp.globalValue.shape
        (2, 12)

        More inconsistent shape problems.

        >>> m = Grid2D(nx=3, ny=3)
        >>> v0 = FaceVariable(mesh=m, rank=1, value=m._orientedFaceNormals)
        >>> print len(v0.dot(m.faceCenters[0]).shape)
        2
        >>> print v0.dot(m.faceCenters[0]).globalValue.shape
        (2, 24)

        """
        rankA = len(A.shape) - 1
        rankB = len(B.shape) - 1

        index = (numerix.index_exp[...] + (numerix.newaxis,) * (rankB - 1)
                 + numerix.index_exp[:])
        opShape = numerix._broadcastShape(A[index].shape, B.shape)

        if rankA > 0 and rankB > (rankA - 1):
            opShape = opShape[:rankA-1] + opShape[rankA:]

        return A._BinaryOperatorVariable(lambda a,b: _MeshVariable._dot(a, b, index),
                                         B,
                                         opShape=opShape,
                                         operatorClass=operatorClass,
                                         canInline=False)
    __dot = staticmethod(__dot)

    def dot(self, other, opShape=None, operatorClass=None):
        """
        Return the mesh-element--by--mesh-element (cell-by-cell, face-by-face,
        etc.) scalar product

        .. math::

           \text{self} \cdot \text{other}

        Both `self` and `other` can be of arbitrary rank, and `other` does not
        need to be a `_MeshVariable`.
        """
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)
        opShape, baseClass, other = self._shapeClassAndOther(opShape=None, operatorClass=None, other=other)

        return _MeshVariable.__dot(self, other, self._OperatorVariableClass(baseClass))

    def rdot(self, other, opShape=None, operatorClass=None):
        """
        Return the mesh-element--by--mesh-element (cell-by-cell, face-by-face,
        etc.) scalar product

        .. math::

           \text{other} \cdot \text{self}

        Both `self` and `other` can be of arbitrary rank, and `other` does not
        need to be a `_MeshVariable`.
        """
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)
        opShape, baseClass, other = self._shapeClassAndOther(opShape=None, operatorClass=None, other=other)

        return self.__dot(other, self, self._OperatorVariableClass(baseClass))

    def _maxminparallel_(self, a, axis, default, fn, fnParallel):
        a = a[..., self._localNonOverlappingIDs]

        if numerix.multiply.reduce(a.shape) == 0:
            if axis is None:
                opShape = ()
            else:
                opShape=self.shape[:axis] + self.shape[axis+1:]

            if len(opShape) == 0:
                nodeVal = default
            else:
                nodeVal = numerix.empty(opShape)
                nodeVal[:] = default
        else:
            nodeVal = fn(axis=axis)

        return fnParallel(nodeVal)

    def max(self, axis=None):
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def maxParallel(a):
                return self._maxminparallel_(a=a, axis=axis, default=-numerix.inf,
                                             fn=a.max, fnParallel=self.mesh.communicator.MaxAll)

            return self._axisOperator(opname="maxVar",
                                      op=maxParallel,
                                      axis=axis)
        else:
            return Variable.max(self, axis=axis)

    def min(self, axis=None):
        """
        >>> from fipy import Grid2D, CellVariable
        >>> mesh = Grid2D(nx=5, ny=5)
        >>> x, y = mesh.cellCenters
        >>> v = CellVariable(mesh=mesh, value=x*y)
        >>> print v.min()
        0.25
        """
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def minParallel(a):
                return self._maxminparallel_(a=a, axis=axis, default=numerix.inf,
                                             fn=a.min, fnParallel=self.mesh.communicator.MinAll)

            return self._axisOperator(opname="minVar",
                                      op=minParallel,
                                      axis=axis)
        else:
            return Variable.min(self, axis=axis)

    def all(self, axis=None):
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def allParallel(a):
                a = a[..., self._localNonOverlappingIDs]
                return self.mesh.communicator.all(a, axis=axis)

            return self._axisOperator(opname="allVar",
                                      op=allParallel,
                                      axis=axis)
        else:
            return Variable.all(self, axis=axis)

    def any(self, axis=None):
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def anyParallel(a):
                a = a[..., self._localNonOverlappingIDs]
                return self.mesh.communicator.any(a, axis=axis)

            return self._axisOperator(opname="anyVar",
                                      op=anyParallel,
                                      axis=axis)
        else:
            return Variable.any(self, axis=axis)

    def sum(self, axis=None):
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def sumParallel(a):
                a = a[..., self._localNonOverlappingIDs]
                return self.mesh.communicator.sum(a, axis=axis)

            return self._axisOperator(opname="sumVar",
                                      op=sumParallel,
                                      axis=axis)
        else:
            return Variable.sum(self, axis=axis)

    def allclose(self, other, rtol=1.e-5, atol=1.e-8):
         if self.mesh.communicator.Nproc > 1:
             def allcloseParallel(a, b):
                 return self.mesh.communicator.allclose(a, b, rtol=rtol, atol=atol)

             operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
             return self._BinaryOperatorVariable(allcloseParallel,
                                                 other,
                                                 operatorClass=operatorClass,
                                                 opShape=(),
                                                 canInline=False)
         else:
             return Variable.allclose(self, other, rtol=rtol, atol=atol)

    def allequal(self, other):
         if self.mesh.communicator.Nproc > 1:
             def allequalParallel(a, b):
                 return self.mesh.communicator.allequal(a, b)

             operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
             return self._BinaryOperatorVariable(allequalParallel,
                                                 other,
                                                 operatorClass=operatorClass,
                                                 opShape=(),
                                                 canInline=False)
         else:
             return Variable.allequal(self, other)


    def _shapeClassAndOther(self, opShape, operatorClass, other):
        """
        Determine the shape of the result, the base class of the result, and (if
        necessary) a modified form of `other` that is suitable for the
        operation.

        By default, returns the result of the generic
        `Variable._shapeClassAndOther()`, but if that fails, and if each
        dimension of `other` is exactly the `Mesh` dimension, do what the user
        probably "meant" and project `other` onto the `Mesh`.

        >>> from fipy import *
        >>> mesh = Grid1D(nx=5)
        >>> A = numerix.arange(5)
        >>> B = Variable(1.)
        >>> import warnings
        >>> savedFilters = list(warnings.filters)
        >>> warnings.resetwarnings()
        >>> warnings.simplefilter("error", UserWarning, append=True)
        >>> C = (A * B) * CellVariable(mesh=mesh)
        Traceback (most recent call last):
          ...
        UserWarning: The expression `(Variable(value=array(1.0)) * [0 1 2 3 4])` has been cast to a constant `CellVariable`
        >>> warnings.filters = savedFilters
        """
        otherShape = numerix.getShape(other)
        if (not isinstance(other, _MeshVariable)
            and otherShape is not ()
            and otherShape[-1] == self._globalNumberOfElements):
            if (isinstance(other, Variable) and len(other.requiredVariables) > 0):
                import warnings
                warnings.warn("The expression `%s` has been cast to a constant `%s`"
                              % (repr(other), self._variableClass.__name__),
                              UserWarning, stacklevel=4)
            other = self._variableClass(value=other, mesh=self.mesh)

        newOpShape, baseClass, newOther = Variable._shapeClassAndOther(self, opShape, operatorClass, other)

        if ((newOpShape is None or baseClass is None)
            and numerix.alltrue(numerix.array(numerix.getShape(other)) == self.mesh.dim)):
                newOpShape, baseClass, newOther = Variable._shapeClassAndOther(self, opShape, operatorClass, other[..., numerix.newaxis])

        return (newOpShape, baseClass, newOther)

    def _OperatorVariableClass(self, baseClass=None):
        baseClass = Variable._OperatorVariableClass(self, baseClass=baseClass)

        class _MeshOperatorVariable(baseClass):
            def __init__(self, op, var, opShape=None, canInline=True,
                         *args, **kwargs):
                mesh = reduce(lambda a, b: a or b,
                              [getattr(v, "mesh", None) for v in var])
                for shape in [opShape] + [getattr(v, "opShape", None) for v in var]:
                    if shape is not None:
                        opShape = shape
                        break
##                 opShape = reduce(lambda a, b: a or b,
##                                  [opShape] + [getattr(v, "opShape", None) for v in var])
                if opShape is not None:
                    elementshape = opShape[:-1]
                else:
                    elementshape = reduce(lambda a, b: a or b,
                                          [getattr(v, "elementshape", None) for v in var])

                baseClass.__init__(self, mesh=mesh, op=op, var=var,
                                   opShape=opShape, canInline=canInline,
                                   elementshape=elementshape,
                                   *args, **kwargs)

            @property
            def rank(self):
                return len(self.opShape) - 1

        return _MeshOperatorVariable

    @property
    def rank(self):
        return len(self.shape) - 1

    def setValue(self, value, unit = None, where = None):
        if where is not None:
            shape = numerix.getShape(where)
            if shape != self.shape \
              and shape == self._getShapeFromMesh(mesh=self.mesh):
                for dim in self.elementshape:
                    where = numerix.repeat(where[numerix.newaxis, ...], repeats=dim, axis=0)

        return Variable.setValue(self, value=value, unit=unit, where=where)

    def _axisClass(self, axis):
        """
        if we operate along the mesh elements, then this is no longer a
        `_MeshVariable`, otherwise we get back a `_MeshVariable` of the same
        class, but lower rank.
        """
        if axis is None or axis == len(self.shape) - 1 or axis == -1:
            return Variable._OperatorVariableClass(self, baseClass=Variable)
        else:
            return self._OperatorVariableClass()

    def _getitemClass(self, index):
        if not isinstance(index, tuple):
            if isinstance(index, list):
                index = tuple(index)
            else:
                index = (index,)
        indexshape = numerix._indexShape(index=index, arrayShape=self.shape)
        if (len(indexshape) > 0
            and indexshape[-1] == self.shape[-1]
            and numerix.obj2sctype(index[-1]) != numerix.obj2sctype(bool)):
            return self._OperatorVariableClass()
        else:
            return Variable._OperatorVariableClass(self, baseClass=Variable)

    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the
        `_MeshVariable` to persistent storage.
        """
        return {
            'mesh': self.mesh,
            'name': self.name,
            'value': self.value,
            'unit': self.unit,
        }


def _testDot(self):
    """
    >>> from fipy import *
    >>> mesh = Grid2D(nx=2, ny=3)

    >>> s1 = CellVariable(mesh=mesh, value=2)
    >>> s2 = CellVariable(mesh=mesh, value=3)

    >>> v1 = CellVariable(mesh=mesh, rank=1,
    ...                   value=numerix.array([2,3])[..., numerix.newaxis])
    >>> v2 = CellVariable(mesh=mesh, rank=1,
    ...                   value=numerix.array([3,4])[..., numerix.newaxis])

    >>> t21 = CellVariable(mesh=mesh, rank=2,
    ...                    value=numerix.array([[2, 3],
    ...                                         [4, 5]])[..., numerix.newaxis])
    >>> t22 = CellVariable(mesh=mesh, rank=2,
    ...                    value=numerix.array([[3, 4],
    ...                                         [5, 6]])[..., numerix.newaxis])

    >>> t31 = CellVariable(mesh=mesh, rank=3,
    ...                    value=numerix.array([[[3, 4],
    ...                                         [5, 6]],
    ...                                        [[5, 6],
    ...                                         [7, 8]]])[..., numerix.newaxis])
    >>> t32 = CellVariable(mesh=mesh, rank=3,
    ...                    value=numerix.array([[[2, 3],
    ...                                         [4, 5]],
    ...                                        [[4, 5],
    ...                                         [6, 7]]])[..., numerix.newaxis])

    >>> def P(a):
    ...     a = a.globalValue
    ...     print a[...,0], a.shape

    >>> P(v1.dot(v2))
    18 (6,)
    >>> P(v1.dot(t22))
    [21 26] (2, 6)
    >>> P(v1.dot(t31))
    [[21 26]
     [31 36]] (2, 2, 6)

    >>> P(t21.dot(v1))
    [13 23] (2, 6)
    >>> P(t21.dot(t22))
    [[21 26]
     [37 46]] (2, 2, 6)
    >>> P(t21.dot(t31))
    [[[21 26]
      [31 36]]
    <BLANKLINE>
     [[37 46]
      [55 64]]] (2, 2, 2, 6)

    >>> P(t31.dot(v1))
    [[18 28]
     [28 38]] (2, 2, 6)
    >>> P(t31.dot(t21))
    [[[22 29]
      [34 45]]
    <BLANKLINE>
     [[34 45]
      [46 61]]] (2, 2, 2, 6)
    >>> P(t31.dot(t32))
    [[[[22 29]
       [36 43]]
    <BLANKLINE>
      [[34 45]
       [56 67]]]
    <BLANKLINE>
    <BLANKLINE>
     [[[34 45]
       [56 67]]
    <BLANKLINE>
      [[46 61]
       [76 91]]]] (2, 2, 2, 2, 6)

    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
