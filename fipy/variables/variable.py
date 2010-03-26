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

import sys
import os
import inspect

from fipy.tools.dimensions import physicalField
from fipy.tools import numerix
from fipy.tools import parser
from fipy.tools import inline

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
    >>> b()
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
            >>> Variable(value=3, unit="m", array=numerix.zeros((3,2)))
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
            value = value.getValue()
            if hasattr(value, 'copy'):
                value = value.copy()
            unit = None
            
        self._setValue(value=value, unit=unit, array=array)
        
        self.name = name
                
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
        
           >>> print type(numerix.array([1.0, 2.0]) * Variable([1.0, 2.0]))
           <class 'fipy.variables.binaryOperatorVariable.binOp'>

           >>> from scipy.special import gamma as Gamma
           >>> print type(Gamma(Variable([1.0, 2.0])))
           <class 'fipy.variables.unaryOperatorVariable.unOp'>

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

            if len(args) == 1:
                result = args[0]._UnaryOperatorVariable(op=func, opShape=arr.shape)
            elif len(args) == 2:
                result = args[0]._BinaryOperatorVariable(op=func, other=args[1], opShape=arr.shape)
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

        return numerix.array(self.getValue(), t)

##    def _get_array_interface(self):
##        return self._getArray().__array_interface__
     
##    def _set_array_interface(self, value):
##        self._getArray().__array_interface__ = value
         
##    def _del_array_interface(self):
##        del self._getArray().__array_interface__
  
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


    def _getUnitAsOne(self):
        unit = self.getUnit()
        if unit is physicalField._unity:
            return 1.
        else:
            return physicalField.PhysicalField(value=1, unit=unit)

    def _extractUnit(self, value):
        if isinstance(value, physicalField.PhysicalField):
            return value.getUnit()
        else:
            return physicalField._unity 

    def getUnit(self):
        """
        Return the unit object of `self`.
        
            >>> Variable(value="1 m").getUnit()
            <PhysicalUnit m>
        """
        return self._extractUnit(self.getValue())
        
    def setUnit(self, unit):
        """
        Change the unit object of `self` to `unit`
        
            >>> a = Variable(value="1 m")
            >>> a.setUnit("m**2/s")
            >>> print a
            1.0 m**2/s
        """
        if self.value is None:
            self.getValue()

        if isinstance(self.value, physicalField.PhysicalField):
            self.value.setUnit(unit)
        else:
            self.value = physicalField.PhysicalField(value=self.value, unit=unit)

    def inBaseUnits(self):
        """
        Return the value of the `Variable` with all units reduced to 
        their base SI elements.
        
            >>> e = Variable(value="2.7 Hartree*Nav")
            >>> print e.inBaseUnits()
            7088849.01085 kg*m**2/s**2/mol
        """
        value = self.getValue()
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
            >>> print freeze.inUnitsOf('degF')
            32.0 degF
        
        If several units are specified, the return value is a tuple of
        `Variable` instances with with one element per unit such that
        the sum of all quantities in the tuple equals the the original
        quantity and all the values except for the last one are integers.
        This is used to convert to irregular unit systems like
        hour/minute/second.  The original object will not be changed.
        
            >>> t = Variable(value=314159., unit='s')
            >>> [str(element) for element in t.inUnitsOf('d','h','min','s')]
            ['3.0 d', '15.0 h', '15.0 min', '59.0 s']
        """
        value = self.getValue()
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
##         return (self.getValue())[index]
                            
    def getName(self):
        return self.name
        
    def setName(self, name):
        self.name = name
    
    def __str__(self):
        return str(self.getValue())
            
    def __repr__(self):
        if hasattr(self, 'name') and len(self.name) > 0:
            return self.name
        else:
            s = self.__class__.__name__ + '('
            s += 'value=' + `self.getValue()`
            s += ')'
            return s

    def _getCIndexString(self, shape):
        dimensions = len(shape)
        if dimensions == 1:
            return '[i]'
        elif dimensions == 2:
            if shape[-1] == 1:
                return '[j]'
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
         """
         
         identifier = 'var%s' % (id)

         v = self.getValue()

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
##             return identifier + '(0)'         
             return identifier
         else:
             return identifier + self._getCIndexString(shape)

    def tostring(self, max_line_width=75, precision=8, suppress_small=False, separator=' '):
        return numerix.tostring(self.getValue(), 
                                max_line_width=max_line_width,
                                precision=precision, 
                                suppress_small=suppress_small, 
                                separator=separator)
        
    def __setitem__(self, index, value):
        if self.value is None:
            self.getValue()
        self.value[index] = value
        self._markFresh()
        
    def itemset(self, value):
        if self.value is None:
            self.getValue()
        self.value.itemset(value)
        self._markFresh()
        
    def put(self, indices, value):
        if self.value is None:
            self.getValue()
        numerix.put(self.value, indices, value)
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
        return self.getValue()

    def getValue(self):
        """
        "Evaluate" the `Variable` and return its value (longhand)
        
            >>> a = Variable(value=3)
            >>> print a.getValue()
            3
            >>> b = a + 4
            >>> b
            (Variable(value=array(3)) + 4)
            >>> b.getValue()
            7

        """
        
        if self.stale or not self._isCached() or self.value is None:
            value = self._calcValue()
            if self._isCached():
                self._setValue(value=value)
            else:
                self._setValue(value=None)
            self._markFresh()
        else:
            value = self.value
            
        return value

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

    def _setValue(self, value, unit=None, array=None):
        self.value = self._makeValue(value=value, unit=unit, array=array)
     
    def _makeValue(self, value, unit=None, array=None):

        ## --inline code often returns spurious results with noncontiguous
        ## arrays. A test case was put in _execInline(). The best fix turned out to
        ## be here.
        
        if (inline.doInline 
            and hasattr(value, 'iscontiguous') and not value.iscontiguous()):
            value = value.copy()
            
        if isinstance(value, Variable):
            value = value.getValue()
            
        PF = physicalField.PhysicalField

        if not isinstance(value, PF):
            
            if getattr(self, 'value', None) is not None:
                v = self.value
                if isinstance(v, PF):
                    v = self.value.value
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

        if isinstance(value, PF) and value.getUnit().isDimensionless():
            value = value.getNumericValue()
            
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
            >>> a.setValue(3)
            >>> print a
            [3 3 3]

            >>> b = numerix.array((3,4,5))
            >>> a.setValue(b)
            >>> a[:] = 1
            >>> print b
            [3 4 5]

            >>> a.setValue((4,5,6), where=(1, 0))
            Traceback (most recent call last):
                ....
            ValueError: shape mismatch: objects cannot be broadcast to a single shape
            
        """
        if where is not None:
            tmp = numerix.empty(numerix.getShape(where), self.getsctype())
            tmp[:] = value
            tmp = numerix.where(where, tmp, self.getValue())
        else:
            if hasattr(value, 'copy'):
                tmp = value.copy()
            else:
                tmp = value

        value = self._makeValue(value=tmp, unit=unit, array=None)

        if numerix.getShape(self.value) == ():
            self.value.itemset(value)
        else:
            self.value[:] = value
            
        self._markFresh()
        
    def _setNumericValue(self, value):
        if isinstance(self.value, physicalField.PhysicalField):
            self.value.value = value
        else:
            self.value = value
        
    def _getArray(self):
        if isinstance(self.value, physicalField.PhysicalField):
            return self.value._getArray()
        else:
            return self.value
            
    def getNumericValue(self):
        value = self.getValue()
        if isinstance(value, physicalField.PhysicalField):
            return value.getNumericValue()
        else:
            return value
            
    def getShape(self):
        """
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
        if self.value is not None:
            return numerix.getShape(self.value)
        else:
            return ()
            
    shape = property(fget=lambda self: self.getShape(), doc="Tuple of array dimensions.")

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
            self.typecode = numerix.obj2sctype(rep=self.getNumericValue(), default=default)
        
        return self.typecode
    
    def _calcValue(self):
        return self.value
        
    def getSubscribedVariables(self):
        self.subscribedVariables = [sub for sub in self.subscribedVariables if sub() is not None]
        
        return self.subscribedVariables
        
    def __markStale(self):
        for subscriber in self.getSubscribedVariables():
            if subscriber() is not None:
                ## Even though getSubscribedVariables() strips out dead 
                ## references, subscriber() might still be dead due to the 
                ## vagaries of garbage collection and the possibility that 
                ## later subscribedVariables were removed, changing the 
                ## dependencies of this subscriber. 
                ## See <https://www.matforge.org/fipy/ticket/118> for more explanation.
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
            self._markStale()
        else:
            from fipy.variables.constant import _Constant
            var = _Constant(value=var)
            
        return var
            
    def _requiredBy(self, var):
        assert isinstance(var, Variable)
        
        # we retain a weak reference to avoid a memory leak 
        # due to circular references between the subscriber
        # and the subscribee
        import weakref
        self.subscribedVariables.append(weakref.ref(var))
        
    def _getVariableClass(self):
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
        contiguous arrays.  The `mesh.getCellCenters()[1]` command
        introduces a non-contiguous array into the `Variable` and this
        causes the inline routine to return senseless results.
        
            >>> from fipy import Grid2D, CellVariable
            >>> mesh = Grid2D(dx=1., dy=1., nx=2, ny=2)
            >>> var = CellVariable(mesh=mesh, value=0.)
            >>> Y =  mesh.getCellCenters()[1]
            >>> var.setValue(Y + 1.0)
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
        ## valid typecode. If self.value is None then a typecode is
        ## assigned to the Variable by running the calculation without
        ## inlining. The non-inlined result is thus used the first
        ## time through.

        
        if self.value is None and not hasattr(self, 'typecode'):
            self.canInline = False
            argDict['result'] = self.getValue()
            self.canInline = True
            self.typecode = numerix.obj2sctype(argDict['result'])
        else:
            if self.value is None:
                if self.getsctype() == numerix.bool_:
                    argDict['result'] = numerix.empty(dim, numerix.int8)
                else:
                    argDict['result'] = numerix.empty(dim, self.getsctype())
            else:
                argDict['result'] = self.value

            resultShape = argDict['result'].shape

            if resultShape == ():
                argDict['result'] = numerix.reshape(argDict['result'], (1,))

            inline._runInline(string, converters=None, comment=comment, **argDict)

            if resultShape == ():
                argDict['result'] = numerix.reshape(argDict['result'], resultShape)

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
        
        baseClass = baseClass or self._getVariableClass()
        return operatorVariable._OperatorVariableClass(baseClass=baseClass)
            
    def _UnaryOperatorVariable(self, op, operatorClass=None, opShape=None, canInline=True, unit=None):
        """
        Check that getUnit() works for unOp

            >>> (-Variable(value="1 m")).getUnit()
            <PhysicalUnit m>
            
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

        if not self.getUnit().isDimensionless():
            canInline = False

        return unOp(op=op, var=[self], opShape=opShape, canInline=canInline, unit=unit, 
                    inlineComment=inline._operatorVariableComment(canInline=canInline))

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
        
    def _BinaryOperatorVariable(self, op, other, operatorClass=None, opShape=None, canInline=True, unit=None):
        """
        :Parameters:
          - `op`: the operator function to apply (takes two arguments for `self` and `other`)
          - `other`: the quantity to be operated with
          - `operatorClass`: the `Variable` class that the binary operator should inherit from 
          - `opShape`: the shape that should result from the operation
        """
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value=other)

        opShape, baseClass, other = self._shapeClassAndOther(opShape, operatorClass, other)
        
        if opShape is None or baseClass is None:
            return NotImplemented
    
        for v in [self, other]:
            if not v.getUnit().isDimensionless():
                canInline = False
                
        # obtain a general operator class with the desired base class
        operatorClass = operatorClass or self._OperatorVariableClass(baseClass)
        from fipy.variables import binaryOperatorVariable
        binOp = binaryOperatorVariable._BinaryOperatorVariable(operatorClass)
        
        return binOp(op=op, var=[self, other], opShape=opShape, canInline=canInline, unit=unit, 
                     inlineComment=inline._operatorVariableComment(canInline=canInline))
    
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
        return self._BinaryOperatorVariable(lambda a,b: a%b, other)
            
    def __pow__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: pow(a,b), other)
            
    def __rpow__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: pow(b,a), other)
            
    def __div__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: a/b, other)
        
    def __rdiv__(self, other):
        return self._BinaryOperatorVariable(lambda a,b: b/a, other)
            
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
        
        fabs = abs
        return self._UnaryOperatorVariable(lambda a: fabs(a))

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
            >>> a.setValue(4)
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
            >>> a.setValue(4)
            >>> print b()
            1
            >>> a.setValue(5)
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
            >>> a.setValue(5)
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
            >>> a.setValue(4)
            >>> print b()
            1
            >>> a.setValue(5)
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
        >>> from fipy.meshes.grid1D import Grid1D
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
        >>> from fipy.meshes.grid1D import Grid1D
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
        return iter(self.getValue())
        
    def __len__(self):
        return len(self.getValue())
        
    def __float__(self):
        return float(self.getValue())
        
    def __int__(self):
        return int(self.getValue())
        
    def __nonzero__(self):
        """
            >>> print bool(Variable(value=0))
            0
            >>> print bool(Variable(value=(0, 0, 1, 1)))
            Traceback (most recent call last):
                ...
            ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
        """
        return bool(self.getValue())
    
    def any(self, axis=None):
        """
            >>> print Variable(value=0).any()
            0
            >>> print Variable(value=(0, 0, 1, 1)).any()
            1
        """
        return self._UnaryOperatorVariable(lambda a: a.any(axis=axis))

    def all(self, axis=None):
        """
            >>> print Variable(value=(0, 0, 1, 1)).all()
            0
            >>> print Variable(value=(1, 1, 1, 1)).all()
            1
        """
        return self._UnaryOperatorVariable(lambda a: a.all(axis=axis))

    def arccos(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arccos(a))

    def arccosh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arccosh(a))

    def arcsin(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arcsin(a))

    def arcsinh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arcsinh(a))

    def sqrt(self):
        """
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh= Grid1D(nx=3)

            >>> from fipy.variables.cellVariable import CellVariable
            >>> var = CellVariable(mesh=mesh, value=((0., 2., 3.),), rank=1)
            >>> print (var.dot(var)).sqrt()
            [ 0.  2.  3.]
            
        """
        return self._UnaryOperatorVariable(lambda a: numerix.sqrt(a))
        
    def tan(self):
        return self._UnaryOperatorVariable(lambda a: numerix.tan(a))

    def tanh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.tanh(a))

    def arctan(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arctan(a))

    def arctanh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.arctanh(a))
            
    def exp(self):
        return self._UnaryOperatorVariable(lambda a: numerix.exp(a))

    def log(self):
        return self._UnaryOperatorVariable(lambda a: numerix.log(a))

    def log10(self):
        return self._UnaryOperatorVariable(lambda a: numerix.log10(a))

    def sin(self):
        return self._UnaryOperatorVariable(lambda a: numerix.sin(a))
                
    def sinh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.sinh(a))

    def cos(self):
        return self._UnaryOperatorVariable(lambda a: numerix.cos(a))
        
    def cosh(self):
        return self._UnaryOperatorVariable(lambda a: numerix.cosh(a))

    def floor(self):
        return self._UnaryOperatorVariable(lambda a: numerix.floor(a))

    def ceil(self):
        return self._UnaryOperatorVariable(lambda a: numerix.ceil(a))

    def sign(self):
        return self._UnaryOperatorVariable(lambda a: numerix.sign(a), canInline=False)
        
    def conjugate(self):
        return self._UnaryOperatorVariable(lambda a: numerix.conjugate(a), canInline=False)

    def arctan2(self, other):
        return self._BinaryOperatorVariable(lambda a,b: numerix.arctan2(a,b), other)
                
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
        
    def reshape(self, shape):
        return self._BinaryOperatorVariable(lambda a,b: numerix.reshape(a,b), shape, opShape=shape, canInline=False)
        
    def transpose(self):
        """
        .. attention: This routine is deprecated. 
           It is not longer needed.
        """
        import warnings
        warnings.warn("transpose() is no longer needed", DeprecationWarning, stacklevel=2)
        return self

    def _axisClass(self, axis):
        return self._OperatorVariableClass()

    def _axisOperator(self, opname, op, axis=None):
        if not hasattr(self, opname):
            setattr(self, opname, {})
            
        opdict = getattr(self, opname)
        if not opdict.has_key(axis):
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
                                           unit=self.getUnit(),
                                           canInline=False)

    def take(self, ids, axis=0):
        return numerix.take(self.getValue(), ids, axis)

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
           >>> print var.allclose(numerix.zeros(10000))
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

    def getMag(self):
        if not hasattr(self, "mag"):
            self.mag = self.dot(self).sqrt()
            
        return self.mag
    
    def __getstate__(self):
        """
        Used internally to collect the necessary information to ``pickle`` the 
        `Variable` to persistent storage.
        """
        return {
            'value': self.getValue(),
            'unit': self.getUnit(),
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
        

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
