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

from fipy.variables.variable import Variable
from fipy.variables.constant import _Constant
from fipy.tools import numerix
from fipy.tools.decorators import getsetDeprecated

class _MeshVariable(Variable):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    
    Abstract base class for a `Variable` that is defined on a mesh
    """
    def __init__(self, mesh, name='', value=0., rank=None, elementshape=None, 
                 unit=None, cached=1, _bootstrap=False):
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
          - `cached`: whether to cache or always recalculate the value
          - `_bootstrap`: if `True`, accept supplied value as given, without 
            attempting validation. (only useful during unpickling and `Mesh` creation). 
            Default: `False`
        """
        self.mesh = mesh
        if _bootstrap:
            array = value
        else:
            value = self._globalToLocalValue(value)
            if value is None:
                array = None
    ##         else:
    ##             array = numerix.zeros(self.elementshape 
    ##                                   + self._getShapeFromMesh(mesh),
    ##                                   numerix.obj2sctype(value))
            elif isinstance(value, Variable):
                name = name or value.name
                unit = None
                if isinstance(value, _MeshVariable):
                    if not isinstance(value, self._variableClass):
                        raise TypeError, "A '%s' cannot be cast to a '%s'" % (value._variableClass.__name__, 
                                                                              self._variableClass.__name__)
                    if mesh is None:
                        mesh = value.mesh
                    elif mesh != value.mesh:
                        raise ValueError, "The new 'Variable' must use the same mesh as the supplied value"
                        
                    if elementshape is not None and elementshape != value.shape[:-1]:
                        raise ValueError, "'elementshape' != shape of elements of 'value'"

                    if rank is not None and rank != value.rank:
                        raise ValueError, "'rank' != rank of 'value'"

                    elementshape = value.shape[:-1]
                    array = None

                value = value._copyValue()
            elif isinstance(value, numerix.ndarray) or numerix.MA.isMaskedArray(value):
                if value.shape != () and value.shape[-1] == self._getShapeFromMesh(mesh)[-1]:
                    if elementshape is not None and elementshape != value.shape[:-1]:
                        raise ValueError, "'elementshape' != shape of elements of 'value'"

                    if rank is not None and rank != len(value.shape[:-1]):
                        raise ValueError, "'rank' != rank of 'value'"
                    elementshape = value.shape[:-1]
                elif rank is None and elementshape is None:
                    elementshape = value.shape

            if rank is None:
                if elementshape is None:
                    elementshape = ()
            elif elementshape is None:
                elementshape = rank * (mesh.dim,)
            elif len(elementshape) != rank:
                raise ValueError, 'len(elementshape) != rank'
                    
            self.elementshape = elementshape
            
            if not locals().has_key("array"):
                if numerix.MA.isMaskedArray(value):
                    arrayMaker = numerix.MA.zeros
                else:
                    arrayMaker = numerix.zeros
                    
                array = arrayMaker(self.elementshape 
                                   + self._getShapeFromMesh(mesh),
                                   numerix.obj2sctype(value))
                                   
                if numerix._broadcastShape(array.shape, numerix.shape(value)) is None:
                    if not isinstance(value, Variable):
                        value = _Constant(value)
                    value = value[..., numerix.newaxis]
                                      
    ##         if isinstance(value, _MeshVariable):
    ##             mesh = mesh or value.mesh
            
        self.mesh._subscribe(self)

        Variable.__init__(self, name=name, value=value, unit=unit, 
                          array=array, cached=cached)

    def _isSolvable(self):
        """
        Is this `Variable` suitable for passing to `solve()` or `sweep()`?
        """
        return False
        
    def copy(self):
        return self._variableClass(mesh=self.mesh, 
                                   name=self.name, 
                                   value=self)

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
        
    @getsetDeprecated
    def _getGlobalNumberOfElements(self):
        return self._globalNumberOfElements

    @getsetDeprecated
    def _getGlobalOverlappingIDs(self):
        return self._globalOverlappingIDs

    @getsetDeprecated
    def _getLocalNonOverlappingIDs(self):
        return self._localNonOverlappingIDs
        
    @getsetDeprecated
    def getGlobalValue(self):
        """Concatenate and return values from all processors
        """
        return self.globalValue

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

                            
    @getsetDeprecated
    def getMesh(self):
        return self.mesh
        
    def __str__(self):
        return str(self.globalValue)
        
    def __repr__(self):
        if hasattr(self, 'name') and len(self.name) > 0:
            return self.name
        else:
            s = self.__class__.__name__ + '('
            s += 'value=' + `self.globalValue`
            s += ')'
            if len(self.name) == 0:
                s = s[:-1] + ', mesh=' + `self.mesh` + s[-1]
            return s
        
    @property
    def constraintMask(self):
        if hasattr(self, 'constraints'):
            returnMask = False
            for value, mask in self.constraints:
                returnMask = returnMask | numerix.array(mask)
            return self._variableClass(mesh=self.mesh, rank=0, value=returnMask)
        else:
            return None

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
        >>> from fipy.tools import parallel
        >>> print parallel.procID > 0 or numerix.allequal(var.shape, (6,))
        True
        >>> print parallel.procID > 0 or numerix.allequal(var.arithmeticFaceValue.shape, (17,))
        True
        >>> print parallel.procID > 0 or numerix.allequal(var.grad.shape, (2, 6))
        True
        >>> print parallel.procID > 0 or numerix.allequal(var.faceGrad.shape, (2, 17))
        True
        """
        return (Variable._getShape(self)
                or (self.elementshape + self._getShapeFromMesh(self.mesh)) 
                or ())

    def _dot(a, b, index, omit=()):
        """
        Workhorse method to calculate the scalar product
        
        .. math::
        
           \mathsf{a} \cdot \mathsf{b}

        for all but the last index of `a` and `b`. Both `a` and `b` can be of
        arbitrary rank, but at this point, both must be appropriately broadcast
        `array` objects.
        """
        Ashape = a.shape
        Bshape = b.shape
        for axis in omit:
            Ashape = Ashape[:axis] + Ashape[axis+1:]
            Bshape = Bshape[:axis] + Bshape[axis+1:]
        rankA = len(Ashape) - 1
        rankB = len(Bshape) - 1
        if rankA <= 0 or rankB <= 0:
            # if either a or b is scalar, then just do multiplication
            return a[index] * b
        else:
            return numerix.sum(a[index] * b, axis=rankA - 1)
    _dot = staticmethod(_dot)

    def __dot(A, B, operatorClass, omit=()):
        """
        Workhorse method to return a `_BinaryOperatorVariable` that will
        dynamically perform the mesh-element--by--mesh-element (cell-by-cell,
        face-by-face, etc.) scalar product
        
        .. math::
        
           \mathsf{A} \cdot \mathsf{B}
           
        Both `A` and `B` can be of arbitrary rank, but at this point, both must
        be appropriately broadcast `_MeshVariable` objects.
        """
        Ashape = A.shape
        Bshape = B.shape
        for axis in omit:
            Ashape = Ashape[:axis] + Ashape[axis+1:]
            Bshape = Bshape[:axis] + Bshape[axis+1:]
        rankA = len(Ashape) - 1
        rankB = len(Bshape) - 1
        
        index = (numerix.index_exp[...] + (numerix.newaxis,) * (rankB - 1) 
                 + numerix.index_exp[:])
        opShape = numerix._broadcastShape(A[index].shape, Bshape)
        if rankA > 0:
            opShape = opShape[:rankA-1] + opShape[rankA:]
        
        return A._BinaryOperatorVariable(lambda a,b: _MeshVariable._dot(a, b, index, omit=omit), 
                                         B, 
                                         opShape=opShape,
                                         operatorClass=operatorClass,
                                         canInline=False)
    __dot = staticmethod(__dot)

    def dot(self, other, omit=()):
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
        
        return self.__dot(self, other, operatorClass=self._OperatorVariableClass(baseClass), omit=omit)

    def rdot(self, other, omit=()):
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
        
        return self.__dot(other, self, operatorClass=self._OperatorVariableClass(baseClass), omit=omit)

    def cross(self, other, axisa=-1, axisb=-1, axisc=-1, axis=None):
        """
        Return the mesh-element--by--mesh-element (cell-by-cell, face-by-face,
        etc.) vector product
        
        .. raw:: latex
        
           \[ \text{self} \times \text{other} \]
           
        Both `self` and `other` can be of arbitrary rank, and `other` does not
        need to be a `_MeshVariable`.
        
        :Parameters:
          - `axisa`: *unused*
          - `axisa`: *unused*
          - `axisa`: *unused*
          - `axis`: *unused*

        """
        return Variable.cross(self, other, axis=0)

    def rcross(self, other, axisa=-1, axisb=-1, axisc=-1, axis=None):
        """
        Return the mesh-element--by--mesh-element (cell-by-cell, face-by-face,
        etc.) vector product
        
        .. raw:: latex
        
           \[ \text{other} \times \text{self} \]
           
        Both `self` and `other` can be of arbitrary rank, and `other` does not
        need to be a `_MeshVariable`.
        
        :Parameters:
          - `axisa`: *unused*
          - `axisa`: *unused*
          - `axisa`: *unused*
          - `axis`: *unused*
        """
        return Variable.rcross(other, self, axis=0)
        
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
                                             fn=a.max, fnParallel=self.mesh.communicator.epetra_comm.MaxAll)
                
            return self._axisOperator(opname="maxVar", 
                                      op=maxParallel, 
                                      axis=axis)
        else:
            return Variable.max(self, axis=axis)
                                  
    def min(self, axis=None):
        if self.mesh.communicator.Nproc > 1 and (axis is None or axis == len(self.shape) - 1):
            def minParallel(a):
                return self._maxminparallel_(a=a, axis=axis, default=numerix.inf, 
                                             fn=a.min, fnParallel=self.mesh.communicator.epetra_comm.MinAll)
                
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
        >>> warnings.simplefilter("error", UserWarning, append=True)
        >>> C = CellVariable(mesh=mesh) * (A * B)
        Traceback (most recent call last):
          ...
        UserWarning: The expression `(multiply([0 1 2 3 4], Variable(value=array(1.0))))` has been cast to a constant `CellVariable`
        >>> junk = warnings.filters.pop()
        """
        if other is not None:
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
        
        from fipy.variables.indexVariable import _IndexVariable_
        
        if ((newOpShape is None or baseClass is None)
            and not isinstance(other, _IndexVariable_)
            and numerix.alltrue(numerix.array(numerix.getShape(other)) == self.mesh.dim)):
                newOpShape, baseClass, newOther = Variable._shapeClassAndOther(self, opShape, operatorClass, other[..., numerix.newaxis])

        if (baseClass is not None
            and issubclass(baseClass, _MeshVariable) 
            and (len(newOpShape) == 0
                 or newOpShape[-1] != baseClass._getShapeFromMesh(self.mesh)[-1])):
            baseClass = None
            
        return (newOpShape, baseClass, newOther)

    def _OperatorVariableClass(self, baseClass=None):
        baseClass = Variable._OperatorVariableClass(self, baseClass=baseClass)
                                     
        class _MeshOperatorVariable(baseClass):
            def __init__(self, op, var, opShape=None, canInline=True, mesh=None,
                         *args, **kwargs):
                mesh = mesh or reduce(lambda a, b: a or b, 
                                      [getattr(v, "mesh", None) for v in var])
                for shape in [opShape] + [getattr(v, "opShape", None) for v in var]:
                    if shape is not None:
                        opShape = shape
                        break
                if opShape is not None:
                    elementshape = opShape[:-1]
                else:
                    elementshape = reduce(lambda a, b: a or b, 
                                          [getattr(v, "elementshape", None) for v in var])

                baseClass.__init__(self, mesh=mesh, op=op, var=var, 
                                   opShape=opShape, canInline=canInline,
                                   elementshape=elementshape,
                                   *args, **kwargs)
                                 
            def _isSolvable(self):
                """
                Is this `Variable` suitable for passing to `solve()` or `sweep()`?
                """
                return False
                
            @getsetDeprecated
            def getRank(self):
                return self.rank

            @property
            def rank(self):
                return len(self.opShape) - 1
                
        return _MeshOperatorVariable
                          
    def _reshapeClass(self, opShape):
        if opShape[-1] == self.shape[-1]:
            return self._OperatorVariableClass()

    @getsetDeprecated
    def getRank(self):
        return self.rank

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
        ...                   value=array([2,3])[..., newaxis])
        >>> v2 = CellVariable(mesh=mesh, rank=1, 
        ...                   value=array([3,4])[..., newaxis])
        
        >>> t21 = CellVariable(mesh=mesh, rank=2, 
        ...                    value=array([[2, 3],
        ...                                 [4, 5]])[..., newaxis])
        >>> t22 = CellVariable(mesh=mesh, rank=2, 
        ...                    value=array([[3, 4],
        ...                                 [5, 6]])[..., newaxis])

        >>> t31 = CellVariable(mesh=mesh, rank=3, 
        ...                    value=array([[[3, 4],
        ...                                  [5, 6]],
        ...                                 [[5, 6],
        ...                                  [7, 8]]])[..., newaxis])
        >>> t32 = CellVariable(mesh=mesh, rank=3, 
        ...                    value=array([[[2, 3],
        ...                                  [4, 5]],
        ...                                 [[4, 5],
        ...                                  [6, 7]]])[..., newaxis])

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
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
