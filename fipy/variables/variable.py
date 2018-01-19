#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "variable.py"
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

import os

from fipy.tools.dimensions import physicalField
from fipy.tools import numerix
from fipy.tools import parser
from fipy.tools import inline

__all__ = ["Variable"]

class Variable(object):
    """
    Lazily evaluated quantity with units.

    Using a :class:`~fipy.variables.variable.Variable` in a mathematical expression will create an
    automatic dependency :class:`~fipy.variables.variable.Variable`, e.g.,

    >>> a = Variable(value=3)
    >>> b = 4 * a
    >>> b
    (Variable(value=array(3)) * 4)
    >>> b()
    12

    Changes to the value of a :class:`~fipy.variables.variable.Variable` will automatically trigger
    changes in any dependent :class:`~fipy.variables.variable.Variable` objects

    >>> a.setValue(5)
    >>> b
    (Variable(value=array(5)) * 4)
    >>> print b()
    20
    """

    _cacheAlways = (os.getenv("FIPY_CACHE") is not None) or False
    if parser.parse("--no-cache", action="store_true"):
        _cacheAlways = False
    if parser.parse("--cache", action="store_true"):
        _cacheAlways = True

    _cacheNever = False

    def __new__(cls, *args, **kwds):
        return object.__new__(cls)

    def __init__(self, value=0., unit=None, array=None, name='', cached=1):
        """
        Create a `Variable`.

            >>> Variable(value=3)
            Variable(value=array(3))
            >>> Variable(value=3, unit="m")
            Variable(value=PhysicalField(3,'m'))
            >>> Variable(value=3, unit="m", array=numerix.zeros((3,2), 'l'))
            Variable(value=PhysicalField(array([[3, 3],
                   [3, 3],
                   [3, 3]]),'m'))

        :Parameters:
          - `value`: the initial value
          - `unit`: the physical units of the `Variable`
          - `array`: the storage array for the `Variable`
          - `name`: the user-readable name of the `Variable`
          - `cached`: whether to cache or always recalculate the value
        """

        self.requiredVariables = []
        self.subscribedVariables = []

        if isinstance(value, Variable):
            value = value.value
            if hasattr(value, 'copy'):
                value = value.copy()
            unit = None

        self._setValueInternal(value=value, unit=unit, array=array)

        self._name = name

        self._cached = cached

        self.stale = 1
        self._markFresh()

##    __array_priority__ and __array_wrap__ are required to override
##    the default behavior of numpy. If a numpy array and a Variable
##    are in a binary operation and numpy is first, then numpy will,
##    by default, try and do everything it can to get a a raw numpy
##    array out of Variable. __array_wrap__ seems to have been
##    introduced into masked array to fix this issue. __array_wrap__ is
##    called after the operation is done so it could hurt efficiency badly.
##    Something else needs to be done to stop the initial evaluation.

    __array_priority__ = 100.0

    def __array_wrap__(self, arr, context=None):
        """
        Required to prevent numpy not calling the reverse binary operations.
        Both the following tests are examples ufuncs.

        >>> print type(numerix.array([1.0, 2.0]) * Variable([1.0, 2.0])) # doctest: +ELLIPSIS
        <class 'fipy.variables.binaryOperatorVariable...binOp'>

        >>> from scipy.special import gamma as Gamma # doctest: +SCIPY
        >>> print type(Gamma(Variable([1.0, 2.0]))) # doctest: +SCIPY +ELLIPSIS
        <class 'fipy.variables.unaryOperatorVariable...unOp'>
        """
        result = arr

        if context is not None:
            from fipy.variables.constant import _Constant
            (func, args, _) = context
            def __makeVariable(v):
                if not isinstance(v, Variable):
                    v = _Constant(v)
                return v
            args = [__makeVariable(arg) for arg in args]

            cannotInline = ["expi", "logical_and", "logical_or", "logical_not", "logical_xor", "sign",
                            "conjugate", "dot", "allclose", "allequal"]
            if len(args) == 1:
                result = args[0]._UnaryOperatorVariable(op=func, opShape=arr.shape, canInline=func.__name__ not in cannotInline)
            elif len(args) == 2:
                result = args[0]._BinaryOperatorVariable(op=func, other=args[1], opShape=arr.shape, canInline=func.__name__ not in cannotInline)
            else:
                result = NotImplemented

        return result

    def __array__(self, t=None):
        """
        Attempt to convert the `Variable` to a numerix `array` object

            >>> v = Variable(value=[2,3])
            >>> print numerix.array(v)
            [2 3]

        A dimensional `Variable` will convert to the numeric value in its base units

            >>> v = Variable(value=[2,3], unit="mm")
            >>> numerix.array(v)
            array([ 0.002,  0.003])
        """

        return numerix.array(self.value, t)

##    def _get_array_interface(self):
##        return self._array.__array_interface__

##    def _set_array_interface(self, value):
##        self._array.__array_interface__ = value

##    def _del_array_interface(self):
##        del self._array.__array_interface__

##    __array_interface__ = property(_get_array_interface,
##                                   _set_array_interface,
##                                   _del_array_interface,
##                                   "the '__array_inteface__'")

    def copy(self):
        """
        Make an duplicate of the `Variable`

            >>> a = Variable(value=3)
            >>> b = a.copy()
            >>> b
            Variable(value=array(3))

        The duplicate will not reflect changes made to the original

            >>> a.setValue(5)
            >>> b
            Variable(value=array(3))

        Check that this works for arrays.

            >>> a = Variable(value=numerix.array((0,1,2)))
            >>> b = a.copy()
            >>> b
            Variable(value=array([0, 1, 2]))
            >>> a[1] = 3
            >>> b
            Variable(value=array([0, 1, 2]))

        """
        return self._getArithmeticBaseClass()(value=self)

    @property
    def _unitAsOne(self):
        unit = self.unit
        if unit is physicalField._unity:
            return 1.
        else:
            return physicalField.PhysicalField(value=1, unit=unit)

    def _extractUnit(self, value):
        if isinstance(value, physicalField.PhysicalField):
            return value.unit
        else:
            return physicalField._unity

    def _getUnit(self):
        """
        Return the unit object of `self`.

            >>> Variable(value="1 m").unit
            <PhysicalUnit m>
        """
        return self._extractUnit(self.value)

    def _setUnit(self, unit):
        """
        Change the unit object of `self` to `unit`

            >>> a = Variable(value="1 m")
            >>> a.unit = "m**2/s"
            >>> print a
            1.0 m**2/s
        """
        if self._value is None:
            self.value

        if isinstance(self._value, physicalField.PhysicalField):
            self._value.unit = unit
        else:
            self._value = physicalField.PhysicalField(value=self._value, unit=unit)

    unit = property(_getUnit, _setUnit)

    def inBaseUnits(self):
        """
        Return the value of the `Variable` with all units reduced to
        their base SI elements.

            >>> e = Variable(value="2.7 Hartree*Nav")
            >>> print e.inBaseUnits().allclose("7088849.01085 kg*m**2/s**2/mol")
            1
        """
        value = self.value
        if isinstance(value, physicalField.PhysicalField):
            return value.inBaseUnits()
        else:
            return value

    def inUnitsOf(self, *units):
        """
        Returns one or more `Variable` objects that express the same
        physical quantity in different units.  The units are specified by
        strings containing their names.  The units must be compatible with
        the unit of the object.  If one unit is specified, the return value
        is a single `Variable`.

        >>> freeze = Variable('0 degC')
        >>> print freeze.inUnitsOf('degF').allclose("32.0 degF")
        1

        If several units are specified, the return value is a tuple of
        `Variable` instances with with one element per unit such that
        the sum of all quantities in the tuple equals the the original
        quantity and all the values except for the last one are integers.
        This is used to convert to irregular unit systems like
        hour/minute/second.  The original object will not be changed.

        >>> t = Variable(value=314159., unit='s')
        >>> print numerix.allclose([e.allclose(v) for (e, v) in zip(t.inUnitsOf('d','h','min','s'),
        ...                                                         ['3.0 d', '15.0 h', '15.0 min', '59.0 s'])],
        ...                        True)
        1
        """
        value = self.value
        if isinstance(value, physicalField.PhysicalField):
            return value.inUnitsOf(*units)
        else:
            return value

##     def __getitem__(self, index):
##         """
##         "Evaluate" the `Variable` and return the specified element
##
##            >>> ## a = Variable(value=((3.,4.),(5.,6.)), unit="m") + "4 m"
##            >>> ## print a[1,1]
##             10.0 m
##
##         It is an error to slice a `Variable` whose `value` is not sliceable
##
##            >>> ## Variable(value=3)[2]
##             Traceback (most recent call last):
##                   ...
##             IndexError: 0-d arrays can't be indexed
##
##         """
##         return (self.value)[index]

    def _getName(self):
        return self._name

    def _setName(self, name):
        self._name = name

    name = property(_getName, _setName)

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        if hasattr(self, 'name') and len(self.name) > 0:
            return self.name
        else:
            s = self.__class__.__name__ + '('
            s += 'value=' + repr(self.value)
            s += ')'
            return s

    def _getCIndexString(self, shape):
        r"""
        Test for inline issue. (1, ni) shapes were not handled correctly.

        >>> from fipy import *
        >>> mesh = Tri2D(dx=1., dy=1., nx=1, ny=1)
        >>> diffCoeff = FaceVariable(mesh = mesh, value = 1.0)
        >>> normals = FaceVariable(mesh=mesh, rank=1, value=mesh._orientedFaceNormals)
        >>> normalsNthCoeff = diffCoeff[numerix.newaxis] * normals

        First variable value access does not provoke inlining.

        >>> value = normalsNthCoeff.value
        >>> print (normalsNthCoeff == value).all()
        True

        """

        dimensions = len(shape)
        if dimensions == 1:
            return '[i]'
        elif dimensions == 2:
            if shape[-1] == 1:
                return '[j]'
            elif shape[0] == 1:
                return '[i]'
            else:
                return '[i + j * ni]'
        elif dimensions == 3:
            if shape[-1] == 1:
                if shape[-2] == 1:
                    return '[k]'
                else:
                    return '[j + k * nj]'
            elif shape[-2] == 1:
                return '[i + k * ni * nj]'
            else:
                return '[i + j * ni + k * ni * nj]'


    def _getCstring(self, argDict={}, id="", freshen=None):
         """
         Generate the string and dictionary to be used in inline
             >>> (Variable((1)))._getCstring(argDict={})
             'var'

             >>> (Variable((1,2,3,4)))._getCstring(argDict={})
             'var[i]'

             >>> (Variable(((1,2),(3,4))))._getCstring(argDict={})
             'var[i + j * ni]'

             >>> Variable((((1,2),(3,4)),((5,6),(7,8))))._getCstring(argDict={})
             'var[i + j * ni + k * ni * nj]'

             >>> (Variable(1) * Variable((1,2,3)))._getCstring(argDict={})
             '(var0 * var1[i])'

         freshen is ignored

         Testing when a cell variable multiplies an array that has a
         shape, but has olny one element. This works regullarly,
         but fails when inlining.

         >>> from fipy import *
         >>> m = Grid1D(nx=3)
         >>> x = m.cellCenters[0]
         >>> tmp = m.cellCenters[0] * numerix.array(((0.,), (1.,)))[1]
         >>> print numerix.allclose(tmp, x)
         True
         >>> print numerix.allclose(tmp, x)
         True

         """

         identifier = 'var%s' % (id)

         v = self.value

         if type(v) not in (type(numerix.array(1)),):
             varray = numerix.array(v)
         else:
             varray = v

         if len(varray.shape) == 0:
             if varray.dtype in (numerix.array(1).dtype,):
                 argDict[identifier] = int(varray)
             elif varray.dtype in (numerix.array(1.).dtype,):
                 argDict[identifier] = float(varray)
             else:
                 argDict[identifier] = varray
         else:
             argDict[identifier] = varray

         try:
             shape = self.opShape
         except AttributeError:
             shape = self.shape

         if len(shape) == 0:
             return identifier
         elif len(shape) == 1 and shape[0] == 1:
             return identifier + '[0]'
         else:
             return identifier + self._getCIndexString(shape)

    def tostring(self, max_line_width=75, precision=8, suppress_small=False, separator=' '):
        return numerix.tostring(self.value,
                                max_line_width=max_line_width,
                                precision=precision,
                                suppress_small=suppress_small,
                                separator=separator)

    def __setitem__(self, index, value):
        if self._value is None:
            self._getValue()
        self._value[index] = value
        self._markFresh()

    def itemset(self, value):
        if self._value is None:
            self._getValue()
        self._value.itemset(value)
        self._markFresh()

    def put(self, indices, value):
        if self._value is None:
            self._getValue()
        numerix.put(self._value, indices, value)
        self._markFresh()

    def __call__(self):
        """
        "Evaluate" the `Variable` and return its value

            >>> a = Variable(value=3)
            >>> print a()
            3
            >>> b = a + 4
            >>> b
            (Variable(value=array(3)) + 4)
            >>> b()
            7
        """
        return self.value

    def _getValue(self):
        """
        "Evaluate" the `Variable` and return its value (longhand)

            >>> a = Variable(value=3)
            >>> print a.value
            3
            >>> b = a + 4
            >>> b
            (Variable(value=array(3)) + 4)
            >>> b.value
            7

        """

        if self.stale or not self._isCached() or self._value is None:
            value = self._calcValue()
            if self._isCached():
                self._setValueInternal(value=value)
            else:
                self._setValueInternal(value=None)
            self._markFresh()
        else:
            value = self._value

        if len(self.constraints) > 0:
            value = value.copy()
            for constraint in self.constraints:
                if constraint.where is None:
                    value[:] = constraint.value
                else:
                    mask = constraint.where
                    if not hasattr(mask, 'dtype') or mask.dtype != bool:
                        mask = numerix.array(mask, dtype=numerix.NUMERIX.bool)

                    if 0 not in value.shape:
                        try:
                            value[..., mask] = constraint.value
                        except:
                            value[..., mask] = numerix.array(constraint.value)[..., mask]

        return value

    def _setValueProperty(self, newVal):
        """Since `self.setValue` contains optional, named parameters, we will
        punt the property's set method off to that."""
        self.setValue(newVal)

    value = property(_getValue, _setValueProperty)

    @property
    def constraints(self):
        if not hasattr(self, "_constraints"):
            self._constraints = []
        return self._constraints

    def constrain(self, value, where=None):
        """
        Constrain the `Variable` to have a `value` at an index or mask location specified by `where`.

        >>> v = Variable((0,1,2,3))
        >>> v.constrain(2, numerix.array((True, False, False, False)))
        >>> print v
        [2 1 2 3]
        >>> v[:] = 10
        >>> print v
        [ 2 10 10 10]
        >>> v.constrain(5, numerix.array((False, False, True, False)))
        >>> print v
        [ 2 10  5 10]
        >>> v[:] = 6
        >>> print v
        [2 6 5 6]
        >>> v.constrain(8)
        >>> print v
        [8 8 8 8]
        >>> v[:] = 10
        >>> print v
        [8 8 8 8]
        >>> del v.constraints[2]
        >>> print v
        [ 2 10  5 10]

        >>> from fipy.variables.cellVariable import CellVariable
        >>> from fipy.meshes import Grid2D
        >>> m = Grid2D(nx=2, ny=2)
        >>> x, y = m.cellCenters
        >>> v = CellVariable(mesh=m, rank=1, value=(x, y))
        >>> v.constrain(((0.,), (-1.,)), where=m.facesLeft)
        >>> print v.faceValue
        [[ 0.5  1.5  0.5  1.5  0.5  1.5  0.   1.   1.5  0.   1.   1.5]
         [ 0.5  0.5  1.   1.   1.5  1.5 -1.   0.5  0.5 -1.   1.5  1.5]]

        :Parameters:
          - `value`: the value of the constraint
          - `where`: the constraint mask or index specifying the location of the constraint

        """

        from fipy.boundaryConditions.constraint import Constraint
        if not isinstance(value, Constraint):
            value = Constraint(value=value, where=where)

        if not hasattr(self, "_constraints"):
            self._constraints = []
        self._constraints.append(value)
        self._requires(value.value)
        self._markStale()

    def release(self, constraint):
        """Remove `constraint` from `self`

        >>> v = Variable((0,1,2,3))
        >>> v.constrain(2, numerix.array((True, False, False, False)))
        >>> v[:] = 10
        >>> from fipy.boundaryConditions.constraint import Constraint
        >>> c1 = Constraint(5, numerix.array((False, False, True, False)))
        >>> v.constrain(c1)
        >>> v[:] = 6
        >>> v.constrain(8)
        >>> v[:] = 10
        >>> del v.constraints[2]
        >>> v.release(constraint=c1)
        >>> print v
        [ 2 10 10 10]
        """
        self.constraints.remove(constraint)

    def _isCached(self):
        return self._cacheAlways or (self._cached and not self._cacheNever)

    def cacheMe(self, recursive=False):
        self._cached = True
        if recursive:
            for var in self.requiredVariables:
                var.cacheMe(recursive=True)

    def dontCacheMe(self, recursive=False):
        self._cached = False
        if recursive:
            for var in self.requiredVariables:
                var.dontCacheMe(recursive=False)

    def _setValueInternal(self, value, unit=None, array=None):
        self._value = self._makeValue(value=value, unit=unit, array=array)

    def _makeValue(self, value, unit=None, array=None):

        ## --inline code often returns spurious results with noncontiguous
        ## arrays. A test case was put in _execInline(). The best fix turned out to
        ## be here.

        if (inline.doInline
            and hasattr(value, 'iscontiguous') and not value.iscontiguous()):
            value = value.copy()

        if isinstance(value, Variable):
            value = value.value

        PF = physicalField.PhysicalField

        if not isinstance(value, PF):

            if getattr(self, '_value', None) is not None:
                v = self._value
                if isinstance(v, PF):
                    v = self._value.value
                if type(value) in (type(1), type(1.)):
                    if type(v) is type(numerix.array(1)):
                        if v.shape is not ():
##                        if len(v) > 1:
                            value = numerix.resize(value, v.shape).astype(v.dtype)

            if unit is not None or type(value) in [type(''), type(()), type([])]:
                value = PF(value=value, unit=unit, array=array)
            elif array is not None:
                array[:] = value
                value = array
            elif type(value) not in (type(None), type(numerix.array(1)), type(numerix.MA.array(1))):
                value = numerix.array(value)
##                 # numerix does strange things with really large integers.
##                 # Even though Python knows how to do arithmetic with them,
##                 # Numeric converts them to 'O' objects that it then doesn't understand.
##                 if value.typecode() == 'O':
##                     value = numerix.array(float(value))

        if isinstance(value, PF) and value.unit.isDimensionless():
            value = value.numericValue

        return value

    def setValue(self, value, unit=None, where=None):
        """
        Set the value of the Variable. Can take a masked array.

            >>> a = Variable((1,2,3))
            >>> a.setValue(5, where=(1, 0, 1))
            >>> print a
            [5 2 5]

            >>> b = Variable((4,5,6))
            >>> a.setValue(b, where=(1, 0, 1))
            >>> print a
            [4 2 6]
            >>> print b
            [4 5 6]
            >>> a.value = 3
            >>> print a
            [3 3 3]

            >>> b = numerix.array((3,4,5))
            >>> a.value = b
            >>> a[:] = 1
            >>> print b
            [3 4 5]

            >>> a.setValue((4,5,6), where=(1, 0)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ....
            ValueError: shape mismatch: objects cannot be broadcast to a single shape

        """
        if where is not None:
            tmp = numerix.empty(numerix.getShape(where), self.getsctype())
            tmp[:] = value
            tmp = numerix.where(where, tmp, self.value)
        else:
            if hasattr(value, 'copy'):
                tmp = value.copy()
            else:
                tmp = value

        value = self._makeValue(value=tmp, unit=unit, array=None)

        if numerix.getShape(self._value) == ():
            self._value.itemset(value)
        else:
            self._value[:] = value

        self._markFresh()

    def _setNumericValue(self, value):
        if isinstance(self._value, physicalField.PhysicalField):
            self._value.value = value
        else:
            self._value = value

    @property
    def _array(self):
        if isinstance(self._value, physicalField.PhysicalField):
            return self._value._array
        else:
            return self._value

    @property
    def numericValue(self):
        value = self.value
        if isinstance(value, physicalField.PhysicalField):
            return value.numericValue
        else:
            return value

    def _getShape(self):
        """
        Tuple of array dimensions.

        >>> Variable(value=3).shape
        ()
        >>> Variable(value=(3,)).shape
        (1,)
        >>> Variable(value=(3,4)).shape
        (2,)

        >>> Variable(value="3 m").shape
        ()
        >>> Variable(value=(3,), unit="m").shape
        (1,)
        >>> Variable(value=(3,4), unit="m").shape
        (2,)
        """
        if self._value is not None:
            return numerix.getShape(self._value)
        else:
            return ()

    shape = property(_getShape)

    def getsctype(self, default=None):
        """

        Returns the Numpy sctype of the underlying array.

            >>> Variable(1).getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1))
            True
            >>> Variable(1.).getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1.))
            True
            >>> Variable((1,1.)).getsctype() == numerix.NUMERIX.obj2sctype(numerix.array((1., 1.)))
            True

        """

        if not hasattr(self, 'typecode'):
            self.typecode = numerix.obj2sctype(rep=self.numericValue, default=default)

        return self.typecode

    @property
    def itemsize(self):
        return self.value.itemsize

    def _calcValue(self):
        return self._value

    def _calcValueNoInline(self):
        raise NotImplementedError

    def _calcValueInline(self):
        raise NotImplementedError

    def _getSubscribedVariables(self):
        self._subscribedVariables = [sub for sub in self._subscribedVariables if sub() is not None]

        return self._subscribedVariables

    def _setSubscribedVariables(self, sVars):
        self._subscribedVariables = sVars

    subscribedVariables = property(_getSubscribedVariables,
                                   _setSubscribedVariables)

    def __markStale(self):
        for subscriber in self.subscribedVariables:
            if subscriber() is not None:
                ## Even though getSubscribedVariables() strips out dead
                ## references, subscriber() might still be dead due to the
                ## vagaries of garbage collection and the possibility that
                ## later subscribedVariables were removed, changing the
                ## dependencies of this subscriber.
                ## See <https://github.com/usnistgov/fipy/issues/103> for more explanation.
                subscriber()._markStale()

    def _markFresh(self):
        self.stale = 0
        self.__markStale()

    def _markStale(self):
        if not self.stale:
            self.stale = 1
            self.__markStale()

    def _requires(self, var):
        if isinstance(var, Variable):
            self.requiredVariables.append(var)
            var._requiredBy(self)
        else:
            from fipy.variables.constant import _Constant
            var = _Constant(value=var)
        self._markStale()
        return var

    def _requiredBy(self, var):
        assert isinstance(var, Variable)

        # we retain a weak reference to avoid a memory leak
        # due to circular references between the subscriber
        # and the subscribee
        import weakref
        self.subscribedVariables.append(weakref.ref(var))

    @property
    def _variableClass(self):
        return Variable

    def _execInline(self, comment=None):
        """
        Gets the stack from _getCstring() which calls _getRepresentation()

            >>> (Variable((1,2,3,4)) * Variable((5,6,7,8)))._getCstring()
            '(var0[i] * var1[i])'
            >>> (Variable(((1,2),(3,4))) * Variable(((5,6),(7,8))))._getCstring()
            '(var0[i + j * ni] * var1[i + j * ni])'
            >>> (Variable((1,2)) * Variable((5,6)) * Variable((7,8)))._getCstring()
            '((var00[i] * var01[i]) * var1[i])'

        The following test was implemented due to a problem with
        contiguous arrays.  The `mesh.cellCenters[1]` command
        introduces a non-contiguous array into the `Variable` and this
        causes the inline routine to return senseless results.

            >>> from fipy import Grid2D, CellVariable
            >>> mesh = Grid2D(dx=1., dy=1., nx=2, ny=2)
            >>> var = CellVariable(mesh=mesh, value=0.)
            >>> Y =  mesh.cellCenters[1]
            >>> var.value = (Y + 1.0)
            >>> print var - Y
            [ 1.  1.  1.  1.]
        """

        from fipy.tools import inline
        argDict = {}
        string = self._getCstring(argDict=argDict, freshen=True) + ';'

        try:
            shape = self.opShape
        except AttributeError:
            shape = self.shape

        dimensions = len(shape)

        if dimensions == 0:
            string = 'result[0] = ' + string
            dim = ()
        else:
            string = 'result' + self._getCIndexString(shape) + ' = ' + string
            ni = self.opShape[-1]
            argDict['ni'] = ni
            if dimensions == 1:
                dim = (ni)
            else:
                nj = self.opShape[-2]
                argDict['nj'] = nj
                if dimensions == 2:
                    dim =(nj,ni)
                elif dimensions == 3:
                    nk = self.opShape[-3]
                    dim = (nk,nj,ni)
                    argDict['nk'] = nk
                else:
                    raise DimensionError, 'Impossible Dimensions'

        ## Following section makes sure that the result array has a
        ## valid typecode. If self._value is None then a typecode is
        ## assigned to the Variable by running the calculation without
        ## inlining. The non-inlined result is thus used the first
        ## time through.


        if self._value is None and not hasattr(self, 'typecode'):
            self.canInline = False
            argDict['result'] = self.value
            self.canInline = True
            self.typecode = numerix.obj2sctype(argDict['result'])
        else:
            if self._value is None:
                if self.getsctype() == numerix.bool_:
                    argDict['result'] = numerix.empty(dim, numerix.int8)
                else:
                    argDict['result'] = numerix.empty(dim, self.getsctype())
            else:
                argDict['result'] = self._value

            resultShape = argDict['result'].shape

            if resultShape == ():
                argDict['result'] = numerix.reshape(argDict['result'], (1,))

            inline._runInline(string, converters=None, comment=comment, **argDict)

            if resultShape == ():
                argDict['result'] = numerix.reshape(argDict['result'], resultShape)

            if self.getsctype() == numerix.bool_:
                argDict['result'] = numerix.asarray(argDict['result'], dtype=self.getsctype())

        return argDict['result']

    def _broadcastShape(self, other):
        ignore, ignore, broadcastshape = numerix._broadcastShapes(self.shape, numerix.getShape(other))

        return broadcastshape

##         selfshape = self.shape
##         othershape = other.shape
##
##         if len(selfshape) > len(othershape):
##             othershape = (1,) * (len(selfshape) - len(othershape)) + othershape
##         elif len(selfshape) < len(othershape):
##             selfshape = (1,) * (len(othershape) - len(selfshape)) + selfshape
##
##         if numerix.logical_and.reduce([(s == o or s == 1 or o == 1) for s,o in zip(selfshape, othershape)]):
##             return tuple([max(s,o) for s,o in zip(selfshape, othershape)])
##         else:
##             return None

    def _getArithmeticBaseClass(self, other=None):
        """
        Given `self` and `other`, return the desired base class for an operation
        result.
        """
        if other is None:
            return Variable

        if self._broadcastShape(other) is not None:
            # If self and other have the same base class, result has that base class.
            # If self derives from other, result has self's base class.
            # If other derives from self, result has other's base class.
            # If self and other don't have a common base, we don't know how to combine them.
            from fipy.variables.constant import _Constant
            if isinstance(self, other._getArithmeticBaseClass()) or isinstance(other, _Constant):
                return self._getArithmeticBaseClass()
            else:
                return None
        else:
            # If self and other have un-broadcastable shapes, we don't know how to combine them.
            return None

    def _OperatorVariableClass(self, baseClass=None):
        from fipy.variables import operatorVariable

        baseClass = baseClass or self._variableClass
        return operatorVariable._OperatorVariableClass(baseClass=baseClass)

    def _UnaryOperatorVariable(self, op, operatorClass=None, opShape=None, canInline=True, unit=None,
                               valueMattersForUnit=False):
        """
        Check that unit works for unOp

            >>> (-Variable(value="1 m")).unit
            <PhysicalUnit m>

        :Parameters:
          - `op`: the operator function to apply (takes one argument for `self`)
          - `operatorClass`: the `Variable` class that the binary operator should inherit from
          - `opShape`: the shape that should result from the operation
          - `valueMattersForUnit`: whether value of `self` should be used when determining unit,
                                    e.g., ???
        """
        operatorClass = operatorClass or self._OperatorVariableClass()
        from fipy.variables import unaryOperatorVariable
        unOp = unaryOperatorVariable._UnaryOperatorVariable(operatorClass)

        # If the caller has not specified a shape for the result, determine the
        # shape from the base class or from the inputs
        if opShape is None:
            opShape = self.shape

        if opShape is None:
            return NotImplemented

        if not self.unit.isDimensionless():
            canInline = False

        return unOp(op=op, var=[self], opShape=opShape, canInline=canInline, unit=unit,
                    inlineComment=inline._operatorVariableComment(canInline=canInline),
                    valueMattersForUnit=[valueMattersForUnit])

    def _shapeClassAndOther(self, opShape, operatorClass, other):
        """
        Determine the shape of the result, the base class of the result, and (if
        necessary) a modified form of `other` that is suitable for the
        operation.
        """
        # If the caller has not specified a base class for the binop,
        # check if the member Variables know what type of Variable should
        # result from the operation.
        baseClass = operatorClass or self._getArithmeticBaseClass(other)

        # If the caller has not specified a shape for the result, determine the
        # shape from the base class or from the inputs
        if opShape is None:
            opShape = self._broadcastShape(other)

        return (opShape, baseClass, other)

    def _BinaryOperatorVariable(self, op, other, operatorClass=None, opShape=None, canInline=True, unit=None,
                                value0mattersForUnit=False, value1mattersForUnit=False):
        """
        :Parameters:
          - `op`: the operator function to apply (takes two arguments for `self` and `other`)
          - `other`: the quantity to be operated with
          - `operatorClass`: the `Variable` class that the binary operator should inherit from
          - `opShape`: the shape that should result from the operation
          - `value0mattersForUnit`: whether value of `self` should be used when determining unit,
                                    e.g., `__rpow__`
          - `value1mattersForUnit`: whether value of `other` should be used when determining unit,
                                    e.g., `__pow__`
        """
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)

        opShape, baseClass, other = self._shapeClassAndOther(opShape, operatorClass, other)

        if opShape is None or baseClass is None:
            return NotImplemented

        for v in [self, other]:
            if not v.unit.isDimensionless() or len(v.shape) > 3:
                canInline = False

        # obtain a general operator class with the desired base class
        operatorClass = operatorClass or self._OperatorVariableClass(baseClass)
        from fipy.variables import binaryOperatorVariable
        binOp = binaryOperatorVariable._BinaryOperatorVariable(operatorClass)

        return binOp(op=op, var=[self, other], opShape=opShape, canInline=canInline, unit=unit,
                     inlineComment=inline._operatorVariableComment(canInline=canInline),
                     valueMattersForUnit=[value0mattersForUnit, value1mattersForUnit])

    def __add__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return other + self
        else:
            return self._BinaryOperatorVariable(lambda a,b: a+b, other)

    __radd__ = __add__

    def __sub__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return -other + self
        else:
            return self._BinaryOperatorVariable(lambda a,b: a-b, other)

    def __rsub__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: b-a, other)

    def __mul__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return other * self
        else:
            return self._BinaryOperatorVariable(lambda a,b: a*b, other)

    __rmul__ = __mul__

    def __mod__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: numerix.fmod(a, b), other)

    def __pow__(self, other):
        """return self**other, or self raised to power other

        >>> print Variable(1, "mol/l")**3
        1.0 mol**3/l**3
        >>> print (Variable(1, "mol/l")**3).unit
        <PhysicalUnit mol**3/l**3>
        """
        return self._BinaryOperatorVariable(lambda a,b: pow(a,b), other, value1mattersForUnit=True)

    def __rpow__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: pow(b,a), other, value0mattersForUnit=True)

    def __truediv__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: a/b, other)

    __div__ = __truediv__

    def __rtruediv__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: b/a, other)

    __rdiv__ = __rtruediv__

    def __neg__(self):
        return self._UnaryOperatorVariable(lambda a: -a)

    def __pos__(self):
        return self

    def __abs__(self):
        """

        Following test it to fix a bug with C inline string using
        abs() instead of fabs()

            >>> print abs(Variable(2.3) - Variable(1.2))
            1.1

        """
        return self._UnaryOperatorVariable(lambda a: numerix.fabs(a))

    def __invert__(self):
        """
        Returns logical "not" of the `Variable`

            >>> a = Variable(value=True)
            >>> print ~a
            False
        """
        return self._UnaryOperatorVariable(lambda a: ~a)

    def __lt__(self,other):
        """
        Test if a `Variable` is less than another quantity

            >>> a = Variable(value=3)
            >>> b = (a < 4)
            >>> b
            (Variable(value=array(3)) < 4)
            >>> b()
            1
            >>> a.value = 4
            >>> print b()
            0
            >>> print 1000000000000000000 * Variable(1) < 1.
            0
            >>> print 1000 * Variable(1) < 1.
            0


        Python automatically reverses the arguments when necessary

            >>> 4 > Variable(value=3)
            (Variable(value=array(3)) < 4)
        """
        return self._BinaryOperatorVariable(lambda a,b: a<b, other)

    def __le__(self,other):
        """
        Test if a `Variable` is less than or equal to another quantity

            >>> a = Variable(value=3)
            >>> b = (a <= 4)
            >>> b
            (Variable(value=array(3)) <= 4)
            >>> b()
            1
            >>> a.value = 4
            >>> print b()
            1
            >>> a.value = 5
            >>> print b()
            0
        """
        return self._BinaryOperatorVariable(lambda a,b: a<=b, other)

    def __eq__(self,other):
        """
        Test if a `Variable` is equal to another quantity

            >>> a = Variable(value=3)
            >>> b = (a == 4)
            >>> b
            (Variable(value=array(3)) == 4)
            >>> b()
            0
        """
        return self._BinaryOperatorVariable(lambda a,b: a==b, other)

    __hash__ = object.__hash__

    def __ne__(self,other):
        """
        Test if a `Variable` is not equal to another quantity

            >>> a = Variable(value=3)
            >>> b = (a != 4)
            >>> b
            (Variable(value=array(3)) != 4)
            >>> b()
            1
        """
        return self._BinaryOperatorVariable(lambda a,b: a!=b, other)

    def __gt__(self,other):
        """
        Test if a `Variable` is greater than another quantity

            >>> a = Variable(value=3)
            >>> b = (a > 4)
            >>> b
            (Variable(value=array(3)) > 4)
            >>> print b()
            0
            >>> a.value = 5
            >>> print b()
            1
        """
        return self._BinaryOperatorVariable(lambda a,b: a>b, other)

    def __ge__(self,other):
        """
        Test if a `Variable` is greater than or equal to another quantity

            >>> a = Variable(value=3)
            >>> b = (a >= 4)
            >>> b
            (Variable(value=array(3)) >= 4)
            >>> b()
            0
            >>> a.value = 4
            >>> print b()
            1
            >>> a.value = 5
            >>> print b()
            1
        """
        return self._BinaryOperatorVariable(lambda a,b: a>=b, other)

    def __and__(self, other):
        """
        This test case has been added due to a weird bug that was appearing.

        >>> a = Variable(value=(0, 0, 1, 1))
        >>> b = Variable(value=(0, 1, 0, 1))
        >>> print numerix.equal((a == 0) & (b == 1), [False,  True, False, False]).all()
        True
        >>> print a & b
        [0 0 0 1]
        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx=4)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> a = CellVariable(value=(0, 0, 1, 1), mesh=mesh)
        >>> b = CellVariable(value=(0, 1, 0, 1), mesh=mesh)
        >>> print numerix.allequal((a == 0) & (b == 1), [False,  True, False, False])
        True
        >>> print a & b
        [0 0 0 1]
        """
        return self._BinaryOperatorVariable(lambda a,b: a & b, other, canInline=False)

    def __or__(self, other):
        """
        This test case has been added due to a weird bug that was appearing.

        >>> a = Variable(value=(0, 0, 1, 1))
        >>> b = Variable(value=(0, 1, 0, 1))
        >>> print numerix.equal((a == 0) | (b == 1), [True,  True, False, True]).all()
        True
        >>> print a | b
        [0 1 1 1]
        >>> from fipy.meshes import Grid1D
        >>> mesh = Grid1D(nx=4)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> a = CellVariable(value=(0, 0, 1, 1), mesh=mesh)
        >>> b = CellVariable(value=(0, 1, 0, 1), mesh=mesh)
        >>> print numerix.allequal((a == 0) | (b == 1), [True,  True, False, True])
        True
        >>> print a | b
        [0 1 1 1]
        """
        return self._BinaryOperatorVariable(lambda a,b: a | b, other, canInline=False)

    def __iter__(self):
        return iter(self.value)

    def __len__(self):
        return len(self.value)

    def __float__(self):
        return float(self.value)

    def __int__(self):
        return int(self.value)

    def __nonzero__(self):
        """
            >>> print bool(Variable(value=0))
            0
            >>> print bool(Variable(value=(0, 0, 1, 1)))
            Traceback (most recent call last):
                ...
            ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
        """
        return bool(self.value)

    def any(self, axis=None):
        """
            >>> print Variable(value=0).any()
            0
            >>> print Variable(value=(0, 0, 1, 1)).any()
            1
        """
        operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
        return self._UnaryOperatorVariable(lambda a: a.any(axis=axis),
                                           operatorClass=operatorClass,
                                           opShape=(),
                                           canInline=False)

    def all(self, axis=None):
        """
            >>> print Variable(value=(0, 0, 1, 1)).all()
            0
            >>> print Variable(value=(1, 1, 1, 1)).all()
            1
        """
        operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
        return self._UnaryOperatorVariable(lambda a: a.all(axis=axis),
                                           operatorClass=operatorClass,
                                           opShape=(),
                                           canInline=False)

    def dot(self, other, opShape=None, operatorClass=None, axis=0):
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)
        if opShape is None:
            opShape = self._broadcastShape(other)
        return self._BinaryOperatorVariable(lambda a,b: numerix.dot(a,b, axis=axis),
                                            other,
                                            opShape=opShape[:axis]+opShape[axis+1:],
                                            operatorClass=operatorClass,
                                            canInline=False)

    def ravel(self):
        return self.value.ravel()

    def _axisClass(self, axis):
        return self._OperatorVariableClass()

    def _axisOperator(self, opname, op, axis=None):
        if not hasattr(self, opname):
            setattr(self, opname, {})

        opdict = getattr(self, opname)
        if axis not in opdict:
            if axis is None:
                opShape = ()
            else:
                opShape=self.shape[:axis] + self.shape[axis+1:]

            opdict[axis] = self._UnaryOperatorVariable(op,
                                                       operatorClass=self._axisClass(axis=axis),
                                                       opShape=opShape,
                                                       canInline=False)

        return opdict[axis]

    def sum(self, axis=None):
        return self._axisOperator(opname="sumVar",
                                  op=lambda a: numerix.sum(a, axis=axis),
                                  axis=axis)

    def max(self, axis=None):
        return self._axisOperator(opname="maxVar",
                                  op=lambda a: a.max(axis=axis),
                                  axis=axis)

    def min(self, axis=None):
        return self._axisOperator(opname="minVar",
                                  op=lambda a: a.min(axis=axis),
                                  axis=axis)

    def _getitemClass(self, index):
        return self._OperatorVariableClass()

    def __getitem__(self, index):
        """
        "Evaluate" the `Variable` and return the specified element

            >>> a = Variable(value=((3.,4.),(5.,6.)), unit="m") + "4 m"
            >>> print a[1,1]
            10.0 m

        It is an error to slice a `Variable` whose `value` is not sliceable

            >>> Variable(value=3)[2]
            Traceback (most recent call last):
                  ...
            IndexError: 0-d arrays can't be indexed

        """
        return self._UnaryOperatorVariable(lambda a: a[index],
                                           operatorClass=self._getitemClass(index=index),
                                           opShape=numerix._indexShape(index=index, arrayShape=self.shape),
                                           unit=self.unit,
                                           canInline=False)

    def take(self, ids, axis=0):
        return numerix.take(self.value, ids, axis)

    def allclose(self, other, rtol=1.e-5, atol=1.e-8):
        """
           >>> var = Variable((1, 1))
           >>> print var.allclose((1, 1))
           1
           >>> print var.allclose((1,))
           1

        The following test is to check that the system does not run
        out of memory.

           >>> from fipy.tools import numerix
           >>> var = Variable(numerix.ones(10000))
           >>> print var.allclose(numerix.zeros(10000, 'l'))
           False

        """
        operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
        return self._BinaryOperatorVariable(lambda a,b: numerix.allclose(a, b, atol=atol, rtol=rtol),
                                            other,
                                            operatorClass=operatorClass,
                                            opShape=(),
                                            canInline=False)

    def allequal(self, other):
        operatorClass = Variable._OperatorVariableClass(self, baseClass=Variable)
        return self._BinaryOperatorVariable(lambda a,b: numerix.allequal(a,b),
                                            other,
                                            operatorClass=operatorClass,
                                            opShape=(),
                                            canInline=False)

    @property
    def mag(self):
        if not hasattr(self, "_mag"):
            self._mag = numerix.sqrt(self.dot(self))

        return self._mag

    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the
        `Variable` to persistent storage.
        """
        return {
            'value': self.value,
            'unit': self.unit,
            'array': None,
            'name': self.name,
        }

    def __setstate__(self, dict):
        """
        Used internally to create a new `Variable` from ``pickled``
        persistent storage.
        """

        import sys
        self._refcount = sys.getrefcount(self)

        self.__init__(**dict)

    def _test(self):
        """
        Inverse cosine of :math:`x`, :math:`\cos^{-1} x`

        >>> from fipy.tools.numerix import *

        >>> arccos(Variable(value=(0,0.5,1.0)))
        arccos(Variable(value=array([ 0. ,  0.5,  1. ])))

        .. attention::

           the next should really return radians, but doesn't

        >>> print tostring(arccos(Variable(value=(0,0.5,1.0))), precision=3)
        [ 1.571  1.047  0.   ]

        Inverse hyperbolic cosine of :math:`x`, :math:`\cosh^{-1} x`

        >>> arccosh(Variable(value=(1,2,3)))
        arccosh(Variable(value=array([1, 2, 3])))
        >>> print tostring(arccosh(Variable(value=(1,2,3))), precision=3)
        [ 0.     1.317  1.763]

        Inverse sine of :math:`x`, :math:`\sin^{-1} x`

        >>> arcsin(Variable(value=(0,0.5,1.0)))
        arcsin(Variable(value=array([ 0. ,  0.5,  1. ])))

        .. attention::

           the next should really return radians, but doesn't

        >>> print tostring(arcsin(Variable(value=(0,0.5,1.0))), precision=3)
        [ 0.     0.524  1.571]

        Inverse hyperbolic sine of :math:`x`, :math:`\sinh^{-1} x`

        >>> arcsinh(Variable(value=(1,2,3)))
        arcsinh(Variable(value=array([1, 2, 3])))
        >>> print tostring(arcsinh(Variable(value=(1,2,3))), precision=3)
        [ 0.881  1.444  1.818]

        Inverse tangent of :math:`x`, :math:`\tan^{-1} x`

        >>> arctan(Variable(value=(0,0.5,1.0)))
        arctan(Variable(value=array([ 0. ,  0.5,  1. ])))

        .. attention::

           the next should really return radians, but doesn't

        >>> print tostring(arctan(Variable(value=(0,0.5,1.0))), precision=3)
        [ 0.     0.464  0.785]

        Inverse tangent of a ratio :math:`x/y`, :math:`\tan^{-1} \frac{x}{y}`

        >>> arctan2(Variable(value=(0, 1, 2)), 2)
        (arctan2(Variable(value=array([0, 1, 2])), 2))

        .. attention::

           the next should really return radians, but doesn't

        >>> print tostring(arctan2(Variable(value=(0, 1, 2)), 2), precision=3)
        [ 0.     0.464  0.785]

        Inverse hyperbolic tangent of :math:`x`, :math:`\tanh^{-1} x`

        >>> arctanh(Variable(value=(0,0.25,0.5)))
        arctanh(Variable(value=array([ 0.  ,  0.25,  0.5 ])))
        >>> print tostring(arctanh(Variable(value=(0,0.25,0.5))), precision=3)
        [ 0.     0.255  0.549]

        Cosine of :math:`x`, :math:`\cos x`

        >>> cos(Variable(value=(0,2*pi/6,pi/2), unit="rad"))
        cos(Variable(value=PhysicalField(array([ 0.        ,  1.04719755,  1.57079633]),'rad')))
        >>> print tostring(cos(Variable(value=(0,2*pi/6,pi/2), unit="rad")), suppress_small=1)
        [ 1.   0.5  0. ]

        Hyperbolic cosine of :math:`x`, :math:`\cosh x`

        >>> cosh(Variable(value=(0,1,2)))
        cosh(Variable(value=array([0, 1, 2])))
        >>> print tostring(cosh(Variable(value=(0,1,2))), precision=3)
        [ 1.     1.543  3.762]

        Tangent of :math:`x`, :math:`\tan x`

        >>> tan(Variable(value=(0,pi/3,2*pi/3), unit="rad"))
        tan(Variable(value=PhysicalField(array([ 0.        ,  1.04719755,  2.0943951 ]),'rad')))
        >>> print tostring(tan(Variable(value=(0,pi/3,2*pi/3), unit="rad")), precision=3)
        [ 0.     1.732 -1.732]

        Hyperbolic tangent of :math:`x`, :math:`\tanh x`

        >>> tanh(Variable(value=(0,1,2)))
        tanh(Variable(value=array([0, 1, 2])))
        >>> print tostring(tanh(Variable(value=(0,1,2))), precision=3)
        [ 0.     0.762  0.964]

        Base-10 logarithm of :math:`x`, :math:`\log_{10} x`

        >>> log10(Variable(value=(0.1,1,10)))
        log10(Variable(value=array([  0.1,   1. ,  10. ])))
        >>> print log10(Variable(value=(0.1,1,10)))
        [-1.  0.  1.]

        Sine of :math:`x`, :math:`\sin x`

        >>> sin(Variable(value=(0,pi/6,pi/2), unit="rad"))
        sin(Variable(value=PhysicalField(array([ 0.        ,  0.52359878,  1.57079633]),'rad')))
        >>> print sin(Variable(value=(0,pi/6,pi/2), unit="rad"))
        [ 0.   0.5  1. ]

        Hyperbolic sine of :math:`x`, :math:`\sinh x`

        >>> sinh(Variable(value=(0,1,2)))
        sinh(Variable(value=array([0, 1, 2])))
        >>> print tostring(sinh(Variable(value=(0,1,2))), precision=3)
        [ 0.     1.175  3.627]

        Square root of :math:`x`, :math:`\sqrt{x}`

        >>> sqrt(Variable(value=(1, 2, 3), unit="m**2"))
        sqrt(Variable(value=PhysicalField(array([1, 2, 3]),'m**2')))
        >>> print tostring(sqrt(Variable(value=(1, 2, 3), unit="m**2")), precision=3)
        [ 1.     1.414  1.732] m

        The largest integer :math:`\le x`, :math:`\lfloor x \rfloor`

        >>> floor(Variable(value=(-1.5,2,2.5), unit="m**2"))
        floor(Variable(value=PhysicalField(array([-1.5,  2. ,  2.5]),'m**2')))
        >>> print floor(Variable(value=(-1.5,2,2.5), unit="m**2"))
        [-2.  2.  2.] m**2

        The largest integer :math:`\ge x`, :math:`\lceil x \rceil`

        >>> ceil(Variable(value=(-1.5,2,2.5), unit="m**2"))
        ceil(Variable(value=PhysicalField(array([-1.5,  2. ,  2.5]),'m**2')))
        >>> print ceil(Variable(value=(-1.5,2,2.5), unit="m**2"))
        [-1.  2.  3.] m**2

        Natural logarithm of :math:`x`, :math:`\ln x \equiv \log_e x`

        >>> log(Variable(value=(0.1,1,10)))
        log(Variable(value=array([  0.1,   1. ,  10. ])))
        >>> print tostring(log(Variable(value=(0.1,1,10))), precision=3)
        [-2.303  0.     2.303]

        Complex conjugate of :math:`z = x + i y`, :math:`z^\star = x - i y`

        >>> var = conjugate(Variable(value=(3 + 4j, -2j, 10), unit="ohm"))
        >>> print var.unit
        <PhysicalUnit ohm>
        >>> print allclose(var.numericValue, (3 - 4j, 2j, 10))
        1
        """
        pass


def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
