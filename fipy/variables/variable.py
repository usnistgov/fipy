#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "variable.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 10/26/06 {9:19:38 AM} 
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import sys
import os
from fipy.tools import numerix

from fipy.meshes.meshIterator import MeshIterator

from fipy.tools.dimensions import physicalField

 
from fipy.tools import numerix
from fipy.tools import parser

class Variable(object):
    _cacheAlways = (os.getenv("FIPY_CACHE") is not None) or False
    if parser.parse("--no-cache", action = "store_true"):
        _cacheAlways = False
    if parser.parse("--cache", action = "store_true"):
        _cacheAlways = True

    _cacheNever = False
    
    """
    Lazily evaluated quantity with units. 
    
    Using a `Variable` in a mathematical expression will create an automatic
    dependency `Variable`, e.g.,
    
	>>> a = Variable(value = 3)
	>>> b = 4 * a
	>>> b
	(Variable(value = 3) * 4)
	>>> b()
	12
	
    Changes to the value of a `Variable` will automatically trigger changes in
    any dependent `Variable` objects
    
	>>> a.setValue(5)
	>>> b
	(Variable(value = 5) * 4)
	>>> b()
	20
	
    """
    
    def __new__(cls, *args, **kwds):
        return object.__new__(cls)
    
    def __init__(self, value=0., unit = None, array = None, name = '', mesh = None, cached = 1):
	"""
	Create a `Variable`.
	
	    >>> Variable(value = 3)
	    Variable(value = array(3))
	    >>> Variable(value = 3, unit = "m")
	    Variable(value = PhysicalField(3,'m'))
	    >>> Variable(value = 3, unit = "m", array = numerix.zeros((3,2)))
	    Variable(value = PhysicalField(array([[3, 3],
                   [3, 3],
	           [3, 3]]),'m'))

	:Parameters:
	  - `value`: the initial value
	  - `unit`: the physical units of the `Variable`
	  - `array`: the storage array for the `Variable`
	  - `name`: the user-readable name of the `Variable`
	  - `mesh`: the mesh that defines the geometry of this `Variable`

	"""
            
	self.requiredVariables = []
	self.subscribedVariables = []

	if isinstance(value, Variable):
	    name = value.name
	    mesh = value.mesh
	    value = value.getValue()
            if hasattr(value, 'copy'):
                value = value.copy()
	    unit = None
	    array = None
	    
	self._setValue(value = value, unit = unit, array = array)
	
        self.name = name
	self.mesh = mesh
		
        self._cached = cached

	self.stale = 1
	self._markFresh()
        
	self.sumVar = {}
	self.faceDifferences = {}
	self.laplacian = {}
        self.mag = None
        self.sliceVars = {}

    def getMesh(self):
	return self.mesh
	
    def __array__(self, t = None):
	"""
        Attempt to convert the `Variable` to a numerix `array` object
        
            >>> v = Variable(value = [2,3])
            >>> print numerix.array(v)
            [ 2.  3.]
            
        It is an error to convert a dimensional `Variable` to a 
        Numeric `array`
        
            >>> v = Variable(value = [2,3], unit = "m")
            >>> numerix.array(v)
            Traceback (most recent call last):
                ...
            TypeError: Numeric array value must be dimensionless

        Convert a list of 1 element Variables to an array

            >>> numerix.array([Variable(0), Variable(0)])
            [[0,]
             [0,]]
            >>> print Variable(0) + Variable(0)
            0
            >>> numerix.array([Variable(0) + Variable(0), Variable(0)])

            >>> numerix.array([Variable(0), Variable(0) + Variable(0)])
            [[0,]
             [0,]]
        
	"""
	return numerix.array(self.getValue(), t)

    def _get_array_interface(self):
        return self._getArray().__array_interface__
 	
    def _set_array_interface(self, value):
        self._getArray().__array_interface__ = value
 	       
    def _del_array_interface(self):
        del self._getArray().__array_interface__
 	
    __array_interface__ = property(_get_array_interface,
                                   _set_array_interface,
                                   _del_array_interface,
                                   "the '__array_inteface__'")
	
    def copy(self):
	"""
	Make an duplicate of the `Variable`
	
	    >>> a = Variable(value = 3)
	    >>> b = a.copy()
	    >>> b
	    Variable(value = array(3))

	The duplicate will not reflect changes made to the original
	                  
	    >>> a.setValue(5)
	    >>> b
	    Variable(value = array(3))

        Check that this works for arrays.

            >>> a = Variable(value = numerix.array((0,1,2)))
            >>> b = a.copy()
            >>> b
            Variable(value = array([0, 1, 2]))
            >>> a[1] = 3
            >>> b
            Variable(value = array([0, 1, 2]))
            
	"""
	return Variable(value = self)


    def _getUnitAsOne(self):
        if self.getUnit() is physicalField._unity:
            return 1.
        else:
            return physicalField.PhysicalField(value=1, unit=self.getUnit())

    def _extractUnit(self, value):
        if isinstance(value, physicalField.PhysicalField):
	    return value.getUnit()
	else:
            return physicalField._unity 

    def getUnit(self):
	"""
	Return the unit object of `self`.	
	    >>> Variable(value = "1 m").getUnit()
	    <PhysicalUnit m>
	"""
        return self._extractUnit(self.getValue())
	
    def inBaseUnits(self):
	"""
	Return the value of the `Variable` with all units reduced to 
	their base SI elements.
	
	    >>> e = Variable(value = "2.7 Hartree*Nav")
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
        
            >>> t = Variable(value = 314159., unit = 's')
            >>> [str(element) for element in t.inUnitsOf('d','h','min','s')]
            ['3.0 d', '15.0 h', '15.0 min', '59.0 s']
        """
        value = self.getValue()
        if isinstance(value, physicalField.PhysicalField):
            return value.inUnitsOf(*units)
        else:
            return value

    def __getitem__(self, index):
        """    
        "Evaluate" the `Variable` and return the specified element
        
            >>> a = Variable(value = ((3.,4.),(5.,6.)), unit = "m") + "4 m"
            >>> print a[1,1]
            10.0 m

        It is an error to slice a `Variable` whose `value` is not sliceable

            >>> Variable(value = 3)[2]
            Traceback (most recent call last):
                  ...
            IndexError: 0-d arrays can't be indexed

        """
        if isinstance(index, MeshIterator):
            assert index.getMesh() == self.getMesh()
            return self.take(index)
        else:
            return (self.getValue())[index]
                           
    def getName(self):
        return self.name
        
    def setName(self, name):
        self.name = name
    
    def __str__(self):
	return str(self.getValue())
	    
    def __repr__(self):
        s = self.__class__.__name__ + '('
	if len(self.name) > 0:
	    s += 'name = "' + self.name + '", '
	s += 'value = ' + `self.getValue()`
	if self.mesh:
	    s += ', mesh = ' + `self.mesh`
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
                return '[j * ni + i]'
        elif dimensions == 3:
            if shape[-1] == 1:
                if shape[-2] == 1:
                    return '[k]'
                else:
                    return '[k + j * nk]'
            elif shape[-2] == 1:
                return '[k + i * nj * nk]'
            else:
                return '[k + j * nk + i * nj * nk]'

    def _getCstring(self, argDict={}, id = "", freshen=None):
         """
         Generate the string and dictionary to be used in inline
             >>> (Variable((1)))._getCstring(argDict = {})
             'var[0]'
           
             >>> (Variable((1,2,3,4)))._getCstring(argDict = {})
             'var[i]'
       
             >>> (Variable(((1,2),(3,4))))._getCstring(argDict = {})
             'var[j * ni + i]'

             >>> Variable((((1,2),(3,4)),((5,6),(7,8))))._getCstring(argDict = {})
             'var[k + j * nk + i * nj * nk]'

             >>> (Variable(1) * Variable((1,2,3)))._getCstring(argDict = {})
             '(var0[0] * var1[i])'

         freshen is ignored
         """
         
         identifier = 'var%s' % (id)

         v = self.getValue()
         if type(v) not in (type(numerix.array(1)),):
             argDict[identifier] = numerix.array(v)
         else:
             argDict[identifier] = v
             
         try:
             shape = self.opShape
         except AttributeError:
             shape = self.getShape()

         if len(shape) == 0:         
             return identifier + '[0]'
         else:
             return identifier + self._getCIndexString(shape)

    def tostring(self, max_line_width = None, precision = None, suppress_small = None, separator = ' '):
        return numerix.tostring(self.getValue(), 
                                max_line_width = max_line_width,
                                precision = precision, 
                                suppress_small = suppress_small, 
                                separator = separator)
        
    def __setitem__(self, index, value):
        if isinstance(index, MeshIterator):
            assert index.getMesh() == self.getMesh()
            self.put(indices=index, value=value)
        else:
            if self.value is None:
                self.getValue()
            self.value[index] = value
            self._markFresh()
        
    def put(self, indices, value):
        if self.value is None:
            self.getValue()
        numerix.put(self.value, indices, value)
        self._markFresh()
	
    def __call__(self):
	"""
	"Evaluate" the `Variable` and return its value
	
	    >>> a = Variable(value = 3)
	    >>> print a()
	    3
	    >>> b = a + 4
	    >>> b
	    (Variable(value = array(3)) + 4)
	    >>> b()
	    7
	"""
	return self.getValue()

    def getValue(self):
        """
        "Evaluate" the `Variable` and return its value (longhand)
        
            >>> a = Variable(value = 3)
            >>> print a.getValue()
            3
            >>> b = a + 4
            >>> b
            (Variable(value = array(3)) + 4)
            >>> b.getValue()
            7

        """
        
        if self.stale or not self._isCached() or self.value is None:
            value = self._calcValue()
            if self._isCached():
                self._setValue(value = value)
            else:
                self._setValue(value = None)
            self._markFresh()
        else:
            value = self.value
            
        return value

    def _isCached(self):
        return self._cacheAlways or (self._cached and not self._cacheNever)
        
    def cacheMe(self, recursive = False):
        self._cached = True
        if recursive:
            for var in self.requiredVariables:
                var.cacheMe(recursive = True)
                
    def dontCacheMe(self, recursive = False):
        self._cached = False
        if recursive:
            for var in self.requiredVariables:
                var.dontCacheMe(recursive = False)

    def _setValue(self, value, unit = None, array = None):
        self.value = self._makeValue(value = value, unit = unit, array = array)
     
    def _makeValue(self, value, unit = None, array = None):

        ## --inline code often returns spurious results with noncontiguous
        ## arrays. A test case was put in _execInline(). The best fix turned out to
        ## be here.
        
        if hasattr(value, 'iscontiguous') and not value.iscontiguous():
            value = value.copy()
            
        PF = physicalField.PhysicalField

        from fipy.tools.numerix import MA

        if not isinstance(value, PF):
            
            if getattr(self, 'value', None) is not None:
                v = self.value
                if isinstance(v, PF):
                    v = self.value.value
                if type(value) in (type(1), type(1.)):
                    if type(v) is type(numerix.array(1)):
                        if v.shape is not ():
##                        if len(v) > 1:
                            value = numerix.resize(float(value), (len(v),))
                    
            if unit is not None or type(value) in [type(''), type(()), type([])]:
                value = PF(value = value, unit = unit, array = array)
            elif array is not None:
                array[:] = value
                value = array
            elif type(value) not in (type(None), type(numerix.array(1)), type(MA.array(1))):
                value = numerix.array(value)
##                 # numerix does strange things with really large integers.
##                 # Even though Python knows how to do arithmetic with them,
##                 # Numeric converts them to 'O' objects that it then doesn't understand.
##                 if value.typecode() == 'O':
##                     value = numerix.array(float(value))

        if isinstance(value, PF) and value.getUnit().isDimensionless():
            value = value.getNumericValue()
            
        return value

    def setValue(self, value, unit = None, array = None, where = None):
        """
        Set the value of the Variable. Can take a masked array.

            >>> a = Variable((1,2,3))
            >>> a.setValue(5, where = (1, 0, 1))
            >>> print a
            [ 5.  2.  5.]

            >>> b = Variable((4,5,6))
            >>> a.setValue(b, where = (1, 0, 1))
            >>> print a
            [ 4.  2.  6.]
            >>> print b
            [ 4.  5.  6.]
            >>> a.setValue(3)
            >>> print a
            [ 3.  3.  3.]

            >>> b = numerix.array((3,4,5))
            >>> a.setValue(b)
            >>> a[:] = 1
            >>> print b
            [3 4 5]

            >>> a.setValue((4,5,6), where = (1, 0))
            Traceback (most recent call last):
                ....
            ValueError: array dimensions must agree
            
        """

        if hasattr(value, 'copy'):
            tmp = value.copy()
        else:
            tmp = value

        if where is not None:
            tmp = numerix.where(where, tmp, self.getValue())
	self._setValue(value = tmp, unit = unit, array = array)
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
            return self.getValue()
            
    def getNumericValue(self):
	value = self.getValue()
	if isinstance(value, physicalField.PhysicalField):
	    return value.getNumericValue()
	else:
	    return value
            
    def _getShapeFromMesh(mesh):
        """
        Return the shape of this `Variable` type, given a particular mesh.
        Return None if unknown or independent of the mesh.
        """
        return None
    _getShapeFromMesh = staticmethod(_getShapeFromMesh)
        
    def getShape(self):
        """
            >>> Variable(value = 3).getShape()
            ()
            >>> Variable(value = (3,)).getShape()
            (1,)
            >>> Variable(value = (3,4)).getShape()
            (2,)
            
            >>> Variable(value = "3 m").getShape()
            ()
            >>> Variable(value = (3,), unit = "m").getShape()
            (1,)
            >>> Variable(value = (3,4), unit = "m").getShape()
            (2,)

            >>> from fipy.meshes.grid2D import Grid2D
            >>> from fipy.variables.cellVariable import CellVariable
            >>> mesh = Grid2D(nx = 2, ny = 3)
            >>> var = CellVariable(mesh = mesh)
            >>> var.getShape()
            (6,)
            >>> var.getArithmeticFaceValue().getShape()
            (17,)
            >>> var.getGrad().getShape()
            (6, 2)
            >>> var.getFaceGrad().getShape()
            (17, 2)
        """
        if self.value is not None:
            return numerix.getShape(self.value)
        else:
            return self._getShapeFromMesh(self.getMesh()) or ()

    def getTypecode(self):
        """

        Returns the Numpy typecode of the underlying array.

            >>> Variable(1).getTypecode()
            'l'
            >>> Variable(1.).getTypecode()
            'd'
            >>> Variable((1,1.)).getTypecode()
            'd'
            
        """
        
        if not hasattr(self, 'typecode'):
            self.typecode = numerix.getTypecode(self.getValue())
        
        return self.typecode

    def _calcValue(self):
        return self.value
        
    def getSubscribedVariables(self):
        self.subscribedVariables = [sub for sub in self.subscribedVariables if sub() is not None]
        
        return self.subscribedVariables
        
    def __markStale(self):
        for subscriber in self.getSubscribedVariables():
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
        
    def _getOperatorVariableClass(self, baseClass = None, canInline = True):
        """
            >>> a = Variable(value = 1)

            >>> c = -a
            >>> b = c.getOld() + 3
            >>> print b
            2
            >>> b.getTypecode()
            'l'
            
        replacing with the same thing is no problem
        
            >>> a.setValue(3)
            >>> b = c.getOld() + 3
            >>> print b
            0
            
        replacing with multiple copies causes the reference counting problem
        
            >>> a.setValue(3)
            >>> b = (c + c).getOld() + 3
            >>> print b
            -3
            
        the order matters
        
            >>> b = (c + c).getOld() + 3
            >>> a.setValue(2)
            >>> print b
            -1

        
        Test of _getRepresentation

            >>> v1 = Variable(numerix.array((1,2,3,4)))
            >>> v2 = Variable(numerix.array((5,6,7,8)))
            >>> v3 = Variable(numerix.array((9,10,11,12)))
            >>> v4 = Variable(numerix.array((13,14,15,16)))

            >>> (v1 * v2)._getRepresentation()
            '(Variable(value = array([1, 2, 3, 4])) * Variable(value = array([5, 6, 7, 8])))'
            
            >>> (v1 * v2)._getRepresentation(style='C', id = "")
            '(var0[i] * var1[i])'
            
            >>> (v1 * v2 + v3 * v4)._getRepresentation(style='C', id = "")
            '((var00[i] * var01[i]) + (var10[i] * var11[i]))'
            
            >>> (v1 - v2)._getRepresentation(style='C', id = "")
            '(var0[i] - var1[i])'

            >>> (v1 / v2)._getRepresentation(style='C', id = "")
            '(var0[i] / var1[i])'

            >>> (v1 - 1)._getRepresentation(style='C', id = "")
            '(var0[i] - var1[0])'
                
            >>> (5 * v2)._getRepresentation(style='C', id = "")
            '(var0[i] * var1[0])'

            >>> (v1 / v2 - v3 * v4 + v1 * v4)._getRepresentation(style='C', id = "")
            '(((var000[i] / var001[i]) - (var010[i] * var011[i])) + (var10[i] * var11[i]))'

        Check that getUnit() works for a binOp

            >>> (Variable(value = "1 m") * Variable(value = "1 s")).getUnit()
            <PhysicalUnit s*m>

            >>> (Variable(value = "1 m") / Variable(value = "0 s")).getUnit()
            <PhysicalUnit m/s>

            >>> a = -((Variable() * Variable()).sin())

        Check that getTypeCode() works as expected.

            >>> a = Variable(1.) * Variable(1)
            >>> a.getTypecode()
            'd'

        The following test is to correct an `--inline` bug that was
        being thrown by the Cahn-Hilliard example. The fix for this
        bug was to add the last line to the following code in
        `_getRepresentation()`.
        
            >>> ##elif style == "C":
            >>> ##    counter = _popIndex()
            >>> ##    if not self.var[counter]._isCached():
            >>> ##        stack.append(self.var[counter]._getCstring(argDict, id = id + str(counter), freshen=freshen))
            >>> ##        self.var[counter].value = None

        This is the test that fails if the last line above is removed
        from `_getRepresentation()`, the `binOp.getValue()` statement
        below will return `1.0` and not `0.5`.
            
            >>> from fipy import numerix
            >>> def doBCs(binOp):
            ...     unOp1 = -binOp
            ...     print binOp.getValue()
            >>> var = Variable(1.)
            >>> binOp = 1. * var
            >>> unOp = -binOp
            >>> print unOp.getValue()
            -1.0
            >>> doBCs(binOp)
            1.0
            >>> var.setValue(0.5)
            >>> print unOp.getValue()
            -0.5
            >>> unOp2 = -binOp
            >>> print binOp.getValue()
            0.5

        """
	if baseClass is None:
            baseClass = self._getVariableClass()
	class OperatorVariable(baseClass):
	    def __init__(self, op, var, mesh = None, opShape = (), canInline = canInline):
                mesh = mesh or var[0].getMesh() or (len(var) > 1 and var[1].getMesh())
		self.op = op
		self.var = var
                self.opShape = opShape
                self.canInline = canInline  #allows for certain functions to opt out of --inline
                baseClass.__init__(self, value = None, mesh = mesh)
                self.name = ''
                for var in self.var:    #C does not accept units
                    if not var.getUnit().isDimensionless():
                        self.canInline = False
                        break

		for aVar in self.var:
		    self._requires(aVar)
                
                self.old = None
                self.dontCacheMe()

            def _calcValue(self):
                from fipy.tools.inline import inline
                #if not self._isCached():
                if not self.canInline:
                    return self._calcValuePy()
                else:
                    return inline._optionalInline(self._calcValueIn, self._calcValuePy)

            def _calcValueIn(self):
                return self._execInline()

            def _calcValuePy(self):
                pass

            def _isCached(self):
                return (Variable._isCached(self) 
                        or (len(self.subscribedVariables) > 1 and not self._cacheNever))

            def getOld(self):
                if self.old is None:
                    oldVar = []
                    for v in self.var:
                        from fipy.variables.cellVariable import CellVariable
                        if isinstance(v, CellVariable):
                            oldVar.append(v.getOld())
                        else:
                            oldVar.append(v)
                    
                    self.old = self.__class__(op = self.op, var = oldVar, mesh = self.getMesh(), opShape=self.opShape, canInline=self.canInline)
                                  
                return self.old
         
            def _getCstring(self, argDict = {}, id = "", freshen=False):
                if self.canInline: # and not self._isCached():
                    s = self._getRepresentation(style = "C", argDict = argDict, id = id, freshen=freshen)
                else:
                    s = baseClass._getCstring(self, argDict=argDict, id=id)
                if freshen:
                    self._markFresh()
                  
                return s
                    
	    def _getRepresentation(self, style = "__repr__", argDict = {}, id = id, freshen=False):
                """
                :Parameters:
                    
                  - `style`: one of `'__repr__'`, `'name'`, `'TeX'`, `'C'`

                """
                import opcode
                
		bytecodes = [ord(byte) for byte in self.op.func_code.co_code]
	
		def _popIndex():
		    return bytecodes.pop(0) + bytecodes.pop(0) * 256
		
		stack = []
		    
		unop = {
		    10: "+", 11: "-", 12: "not ", 15: "~"
		}
		
		binop = {
		    19: "**", 20: "*", 21: "/", 22: "%", 23: "+", 24: "-", 26: "//", 27: "/",
			    62: "<<", 63: ">>", 64: "&", 65: "^", 66: "|", 106: "=="
		}
		
   		while len(bytecodes) > 0:
		    bytecode = bytecodes.pop(0)
		    if opcode.opname[bytecode] == 'UNARY_CONVERT':
			stack.append("`" + stack.pop() + "`")
		    elif opcode.opname[bytecode] == 'BINARY_SUBSCR':
			stack.append(stack.pop(-2) + "[" + stack.pop() + "]")
		    elif opcode.opname[bytecode] == 'RETURN_VALUE':
                        s = stack.pop()
                        if style == 'C':
                            return s.replace('numerix.', '').replace('arc', 'a')
                        else:
                            return s
                    elif opcode.opname[bytecode] == 'LOAD_CONST':
                        stack.append(self.op.func_code.co_consts[_popIndex()])
		    elif opcode.opname[bytecode] == 'LOAD_ATTR':
			stack.append(stack.pop() + "." + self.op.func_code.co_names[_popIndex()])
		    elif opcode.opname[bytecode] == 'COMPARE_OP':
			stack.append(stack.pop(-2) + " " + opcode.cmp_op[_popIndex()] + " " + stack.pop())
		    elif opcode.opname[bytecode] == 'LOAD_GLOBAL':
                        counter = _popIndex()
			stack.append(self.op.func_code.co_names[counter])
		    elif opcode.opname[bytecode] == 'LOAD_FAST':
                        if style == "__repr__":
                            stack.append(repr(self.var[_popIndex()]))
                        elif style == "name":
                            v = self.var[_popIndex()]
                            if isinstance(v, Variable):
                                name = v.getName()
                                if len(name) > 0:
                                    stack.append(name)
                                else:
                                    # The string form of a variable
                                    # would probably be too long and messy.
                                    # Just give shorthand.
                                    stack.append("%s(...)" % v.__class__.__name__)
                            elif type(v) in (type(1), type(1.)):
                                stack.append(repr(v))
                            else:
                                # The string form of anything but a
                                # number would be too long and messy.
                                # Just give shorthand.
                                stack.append("<...>")
                        elif style == "TeX":
                            raise Exception, "TeX style not yet implemented"
                        elif style == "C":
                            counter = _popIndex()
                            if not self.var[counter]._isCached():
                                stack.append(self.var[counter]._getCstring(argDict, id = id + str(counter), freshen=freshen))
                                self.var[counter].value = None
                            else:
                                stack.append(self.var[counter]._getVariableClass()._getCstring(self.var[counter], argDict, \
                                                                                               id = id + str(counter),\
                                                                                               freshen=False))
                        else:
                            raise SyntaxError, "Unknown style: %s" % style
                    elif opcode.opname[bytecode] == 'CALL_FUNCTION':    
			args = []
			for j in range(bytecodes.pop(1)):
			    # keyword parameters
			    args.insert(0, stack.pop(-2) + " = " + stack.pop())
			for j in range(bytecodes.pop(0)):
			    # positional parameters
			    args.insert(0, stack.pop())
			stack.append(stack.pop() + "(" + ", ".join(args) + ")")
                    elif opcode.opname[bytecode] == 'LOAD_DEREF':
                        free = self.op.func_code.co_cellvars + self.op.func_code.co_freevars
                        stack.append(free[_popIndex()])
		    elif unop.has_key(bytecode):
			stack.append(unop[bytecode] + stack.pop())
		    elif binop.has_key(bytecode):
			stack.append(stack.pop(-2) + " " + binop[bytecode] + " " + stack.pop())
		    else:
			raise SyntaxError, "Unknown bytecode: %s in %s: %s" % (`bytecode`, `[ord(byte) for byte in self.op.func_code.co_code]`,`"FIXME"`)
                    
            def __repr__(self):
                return self._getRepresentation()
            
            def getName(self):
                name = baseClass.getName(self)
                if len(name) == 0:
                    name = self._getRepresentation(style = "name")
                return name
                
	    def copy(self):
                return self.__class__(
                    op = self.op,
                    var = self.var,
                    mesh = self.getMesh(),
                    opShape = self.opShape,
                    canInline = self.canInline)
                    #myType = self.myType

            def getShape(self):
                return baseClass.getShape(self) or self.opShape

	return OperatorVariable
	
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return Variable
            
        if other._getArithmeticBaseClass().__name__ is self._getArithmeticBaseClass().__name__ \
        or other.getShape() in ((), (1,)):
            # operating with a scalar results in the same base
            # class as self.
            return self._getArithmeticBaseClass()
        elif self.getShape() == other.getShape():
            # If self and other have the same base class, result has that base class.
            # If self derives from other, result has self's base class.
            # If other derives from self, result has other's base class.
            # If self and other don't have a common base, we don't know how to combine them.
            from fipy.variables.constant import _Constant
            if isinstance(self, other._getArithmeticBaseClass()) or isinstance(other, _Constant):
                return self._getArithmeticBaseClass()
            elif isinstance(other, self._getArithmeticBaseClass()) or isinstance(self, _Constant):
                return other._getArithmeticBaseClass()
            else:
                return None
        else:
            # If self and other have different shapes, we don't know how to combine them.
            return None

    def _execInline(self):
        """
        Gets the stack from _getCstring() which calls _getRepresentation()
        
            >>> (Variable((1,2,3,4)) * Variable((5,6,7,8)))._getCstring()
            '(var0[i] * var1[i])'
            >>> (Variable(((1,2),(3,4))) * Variable(((5,6),(7,8))))._getCstring()
            '(var0[j * ni + i] * var1[j * ni + i])'
            >>> (Variable((1,2)) * Variable((5,6)) * Variable((7,8)))._getCstring()
            '((var00[i] * var01[i]) * var1[i])'

        The following test was implemented due to a problem with
        contiguous arrays.  The `mesh.getCellCenters()[:,1]` command
        introduces a non-contiguous array into the `Variable` and this
        causes the inline routine to return senseless results.
        
            >>> from fipy import Grid2D, CellVariable
            >>> mesh = Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
            >>> var = CellVariable(mesh = mesh, value = 0.)
            >>> Y =  mesh.getCellCenters()[:,1]
            >>> var.setValue(Y + 1.0)
            >>> print var - Y
            [ 1.  1.  1.  1.]




                                                           
        """
    
        from fipy.tools.inline import inline
        argDict = {}
        string = self._getCstring(argDict = argDict, freshen=True) + ';'
        
        try:
            shape = self.opShape
        except AttributeError:
            shape = self.getShape()

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
            self.typecode = numerix.getTypecode(argDict['result'])
        else:
            if self.value is None:
                argDict['result'] = numerix.empty(dim, self.getTypecode())
            else:
                argDict['result'] = self.value

            inline._runInline(string, converters=None, **argDict)
                
        return argDict['result']

    
    def _getUnaryOperatorVariable(self, op, baseClass = None, canInline = True):
	"""
        Check that getUnit() works fot unOp

            >>> (-Variable(value = "1 m")).getUnit()
            <PhysicalUnit m>
            
	"""
        
	class unOp(self._getOperatorVariableClass(baseClass)):
            def _calcValuePy(self):
                return self.op(self.var[0].getValue())

            def getUnit(self):
                try:
                    return self._extractUnit(self.op(self.var[0]._getUnitAsOne()))
                except:
                    return self._extractUnit(self.op(self._calcValue()))
                
        if not self.getUnit().isDimensionless():
            canInline = False

	return unOp(op = op, var = [self], opShape = self.getShape(), canInline = canInline)

    def _getArrayAsOnes(object, valueMattersForShape = ()): 
        """ 
        For the purposes of assembling the binop, we are only
        interested in the shape of the operation result, not the result
        itself.  Some operations (e.g. division) will fail if an input
        happens to have been initialized with zeros, even though it
        will not actually contain zeros by the time a value is
        requested.  Setting the arrays of `self` and `other` to 1
        should always pass?
        
        reshape() is one case where the value cannot be substituted, so
        the shape must be included in valueMattersForShape.
        """
        
        if object not in valueMattersForShape:
            a = numerix.ones(object.getShape())
        else:
            a = object._getArray()
            
        return a

    _getArrayAsOnes = staticmethod(_getArrayAsOnes)

    def _rotateShape(op, var0, var1, var0array, var1array, opShape):
        """
        A scalar `Variable` multiplying/dividing a vector `Variable` will
        fail because the scalar field has shape (N,) and the vector field has shape (N, D)
        This manipulation will give the scalar field shape (N, 1), which will
        allow the desired operator shape of (N, D).
        
        We *only* do this rotation if var1array is rank 1.
        """
        if len(var1array.shape) == 1:
            try:
                if numerix.getShape(op(var0array, var1array[..., numerix.NewAxis])) != opShape:
                    raise ValueError
                from fipy.variables.newAxisVariable import _NewAxisVariable
                var1 = _NewAxisVariable(var1)
            except (ValueError, IndexError):
                raise SyntaxError
        else:
            raise SyntaxError
                
        return (var0, var1)
    _rotateShape = staticmethod(_rotateShape)
    
    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass, rotateShape = True):
        try:
            # check if applying the operation to the inputs will produce the desired shape
            if numerix.getShape(op(var0Array, var1Array)) != opShape:
                raise ValueError

        except ValueError:
            if rotateShape:
                try:
                    # check if changing var1 from a row variable to a column variable
                    # will produce the desired shape
                    (var0, var1) = self._rotateShape(op, var0, var1, var0Array, var1Array, opShape)
                except SyntaxError:
                    if not (otherClass and issubclass(otherClass, Variable)):
                        # check if changing var0 from a row variable to a column variable
                        # will produce the desired shape
                        (var1, var0) = self._rotateShape(op, var1, var0, var1Array, var0Array, opShape)
                    else:
                        raise SyntaxError
            else:
                raise ValueError
            
        return (var0, var1)

    def _getBinaryOperatorVariable(self, op, other, baseClass = None, opShape = None, valueMattersForShape = (), rotateShape = True, canInline = True):
        """
        :Parameters:
          - `op`: the operator function to apply (takes two arguments for `self` and `other`)
          - `other`: the quantity to be operated with
          - `baseClass`: the `Variable` class that the binary operator should inherit from 
          - `opShape`: the shape that should result from the operation
          - `valueMattersForShape`: tuple of elements that must have a particular value for the operation to succeed.
          - `rotateShape`: whether the operator should permit rotation of the variable's shape to allow the operation to complete. This is required because some Numeric operators such as allclose() run out of memory.
        """
        
        # for convenience, we want to be able to treat `other` as a Variable
        # so we record its original class for later reference
        if type(other) is type(numerix.array(1)):
            otherClass = None
        else:
            otherClass = other.__class__
        
        if not isinstance(other, Variable):
            from fipy.variables.constant import _Constant
            other = _Constant(value = other)

        # If the caller has not specified a base class for the binop, 
        # check if the member Variables know what type of Variable should
        # result from the operation.
        baseClass = baseClass or self._getArithmeticBaseClass(other) or other._getArithmeticBaseClass(self)


        # This operation is unknown. Fall back on Python's reciprocal operation or error.
        if baseClass is None:
            return NotImplemented

        mesh = self.getMesh() or other.getMesh()
        
        # If the caller has not specified a shape for the result, determine the 
        # shape from the base class or from the inputs
        opShape = opShape or baseClass._getShapeFromMesh(mesh) or self.getShape() or other.getShape()

        # the magic value of "number" specifies that the operation should result in a single value,
        # regardless of the shapes of the inputs. This hack is necessary because "() or ..." is treated
        # identically to "None or ...".
        if opShape == "number":
            opShape = ()

        var0 = self
        var1 = other
        
        selfArray = self._getArrayAsOnes(object = self, valueMattersForShape = valueMattersForShape)
        otherArray = self._getArrayAsOnes(object = other, valueMattersForShape = valueMattersForShape)

        try:
            var0, var1 = self._verifyShape(op, var0, var1, selfArray, otherArray, opShape, otherClass, rotateShape)
        except SyntaxError:
            return NotImplemented
            
        
        # obtain a general operator class with the desired base class
	operatorClass = self._getOperatorVariableClass(baseClass)
        
        
        # declare a binary operator class with the desired base class
	class binOp(operatorClass):

            def _calcValuePy(self):
                if isinstance(self.var[1], Variable):
                    val1 = self.var[1].getValue()
                else:
                    if type(self.var[1]) is type(''):
                        self.var[1] = physicalField.PhysicalField(value = self.var[1])
                    val1 = self.var[1]

                return self.op(self.var[0].getValue(), val1)

            def getUnit(self):
                try:
                    return self._extractUnit(self.op(self.var[0]._getUnitAsOne(), self.var[1]._getUnitAsOne()))
                except:
                    return self._extractUnit(self._calcValuePy())

            def _getRepresentation(self, style = "__repr__", argDict = {}, id = id, freshen=False):
                self.id = id
                return "(" + operatorClass._getRepresentation(self, style = style, argDict = argDict, id = id, freshen=freshen) + ")"
        
        var = [var0, var1]
        for v in var:
            if not v.getUnit().isDimensionless():
                canInline = False
        tmpBop = binOp(op = op, var = [var0, var1], opShape = opShape, canInline = canInline)
        return tmpBop
    
    def __add__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return other + self
        else:
            return self._getBinaryOperatorVariable(lambda a,b: a+b, other)
	
    __radd__ = __add__

    def __sub__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return -other + self
        else:
            return self._getBinaryOperatorVariable(lambda a,b: a-b, other)
	
    def __rsub__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: b-a, other)
	    
    def __mul__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: a*b, other)
	
    __rmul__ = __mul__
	    
    def __mod__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: a%b, other)
	    
    def __pow__(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: pow(a,b), other)
	#return self._getBinaryOperatorVariable(lambda a,b: a**b, other, canInline = False)
	    
    def __rpow__(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: pow(b,a), other)
        #return self._getBinaryOperatorVariable(lambda a,b: b**a, other, canInline = False)
	    
    def __div__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: a/b, other)
	
    def __rdiv__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: b/a, other)
	    
    def __neg__(self):
	return self._getUnaryOperatorVariable(lambda a: -a)
	
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
	return self._getUnaryOperatorVariable(lambda a: fabs(a))

    def __lt__(self,other):
	"""
	Test if a `Variable` is less than another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a < 4)
	    >>> b
	    (Variable(value = array(3)) < 4)
	    >>> b()
	    1
	    >>> a.setValue(4)
	    >>> b()
	    0
            >>> print 1000000000000000000 * Variable(1) < 1.
            0
            >>> print 1000 * Variable(1) < 1.
            0


	Python automatically reverses the arguments when necessary
	
	    >>> 4 > Variable(value = 3)
	    (Variable(value = array(3)) < 4)
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a<b, other)

    def __le__(self,other):
	"""
	Test if a `Variable` is less than or equal to another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a <= 4)
	    >>> b
	    (Variable(value = array(3)) <= 4)
	    >>> b()
	    1
	    >>> a.setValue(4)
	    >>> b()
	    1
	    >>> a.setValue(5)
	    >>> b()
	    0
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a<=b, other)
	
    def __eq__(self,other):
	"""
	Test if a `Variable` is equal to another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a == 4)
	    >>> b
	    (Variable(value = array(3)) == 4)
	    >>> b()
	    0
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a==b, other)
	
    def __ne__(self,other):
	"""
	Test if a `Variable` is not equal to another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a != 4)
	    >>> b
	    (Variable(value = array(3)) != 4)
	    >>> b()
	    1
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a!=b, other)
	
    def __gt__(self,other):
	"""
	Test if a `Variable` is greater than another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a > 4)
	    >>> b
	    (Variable(value = array(3)) > 4)
	    >>> b()
	    0
	    >>> a.setValue(5)
	    >>> b()
	    1
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a>b, other)
	
    def __ge__(self,other):
	"""
	Test if a `Variable` is greater than or equal to another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a >= 4)
	    >>> b
	    (Variable(value = array(3)) >= 4)
	    >>> b()
	    0
	    >>> a.setValue(4)
	    >>> b()
	    1
	    >>> a.setValue(5)
	    >>> b()
	    1
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a>=b, other)

    def __and__(self, other):
        """
        This test case has been added due to a weird bug that was appearing.

            >>> a = Variable(value = (0, 0, 1, 1))
            >>> b = Variable(value = (0, 1, 0, 1))
            >>> print (a == 0) & (b == 1)
            [0 1 0 0]
            >>> print a & b
            [0 0 0 1]
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 4)
            >>> from fipy.variables.cellVariable import CellVariable
            >>> a = CellVariable(value = (0, 0, 1, 1), mesh = mesh)
            >>> b = CellVariable(value = (0, 1, 0, 1), mesh = mesh)
            >>> print (a == 0) & (b == 1)
            [0 1 0 0]
            >>> print a & b
            [0 0 0 1]

        """
        return self._getBinaryOperatorVariable(lambda a,b: a.astype('h') & b.astype('h'), other, canInline = False)

    def __or__(self, other):
        """
        This test case has been added due to a weird bug that was appearing.

            >>> a = Variable(value = (0, 0, 1, 1))
            >>> b = Variable(value = (0, 1, 0, 1))
            >>> print (a == 0) | (b == 1)
            [1 1 0 1]
            >>> print a | b
            [0 1 1 1]
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 4)
            >>> from fipy.variables.cellVariable import CellVariable
            >>> a = CellVariable(value = (0, 0, 1, 1), mesh = mesh)
            >>> b = CellVariable(value = (0, 1, 0, 1), mesh = mesh)
            >>> print (a == 0) | (b == 1)
            [1 1 0 1]
            >>> print a | b
            [0 1 1 1]
            
        """
        
        return self._getBinaryOperatorVariable(lambda a,b: a.astype('h') | b.astype('h'), other, canInline = False)
        
    def __len__(self):
        return len(self.getValue())
	
    def __float__(self):
        return float(self.getValue())

    def arccos(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arccos(a))

    def arccosh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arccosh(a))

    def arcsin(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arcsin(a))

    def arcsinh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arcsinh(a))

    def sqrt(self):
        """
        
            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh= Grid1D(nx=3)

            >>> from fipy.variables.vectorCellVariable import VectorCellVariable
            >>> var = VectorCellVariable(mesh=mesh, value=((0.,),(2.,),(3.,)))
            >>> print (var.dot(var)).sqrt()
            [ 0.  2.  3.]
            
        """
	return self._getUnaryOperatorVariable(lambda a: numerix.sqrt(a))
	
    def tan(self):
	return self._getUnaryOperatorVariable(lambda a: numerix.tan(a))

    def tanh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.tanh(a))

    def arctan(self):
	return self._getUnaryOperatorVariable(lambda a: numerix.arctan(a))

    def arctanh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arctanh(a))
            
    def exp(self):
	return self._getUnaryOperatorVariable(lambda a: numerix.exp(a))

    def log(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.log(a))

    def log10(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.log10(a))

    def sin(self):
	return self._getUnaryOperatorVariable(lambda a: numerix.sin(a))
		
    def sinh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.sinh(a))

    def cos(self):
	return self._getUnaryOperatorVariable(lambda a: numerix.cos(a))
        
    def cosh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.cosh(a))

    def floor(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.floor(a))

    def ceil(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.ceil(a))
        
    def conjugate(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.conjugate(a), canInline = False)

    def arctan2(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.arctan2(a,b), other)
		
    def dot(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: numerix.dot(a,b), other, canInline = False)
        
    def reshape(self, shape):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.reshape(a,b), shape, valueMattersForShape = (shape,), canInline = False)
        
    def transpose(self):
        """
        .. attention: This routine is deprecated. 
           It is not longer needed.
        """
        import warnings
        warnings.warn("transpose() is no longer needed", DeprecationWarning, stacklevel=2)
        return self

    def sum(self, index = 0):
	if not self.sumVar.has_key(index):
	    from sumVariable import _SumVariable
	    self.sumVar[index] = _SumVariable(self, index)
	
	return self.sumVar[index]

    def take(self, ids, axis = 0):
	return numerix.take(self.getValue(), ids, axis)

    def _take(self, ids, axis = 0):
        """
        
        Same as take() but returns a `Variable` subclass.  This function
        has not yet been implemented as a binary operator but is a
        unary operator.  As a unary operator it has to return the same
        shape as the `Variable` it is acting on.  This is not a
        particular useful implementation of take as it stands. It is
        good for axis permutations.
        

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 1, ny = 1)
           >>> from fipy.variables.vectorFaceVariable import VectorFaceVariable
           >>> var = VectorFaceVariable(value = ( (1, 2), (2, 3), (3, 4), (4, 5) ), mesh = mesh)
           >>> v10 = var._take((1, 0), axis = 1)
           >>> print v10
           [[ 2.  1.]
            [ 3.  2.]
            [ 4.  3.]
            [ 5.  4.]]
           >>> var[3, 0] = 1
           >>> print v10
           [[ 2.  1.]
            [ 3.  2.]
            [ 4.  3.]
            [ 5.  1.]]
           >>> isinstance(var, VectorFaceVariable)
           True
           >>> v0 = var._take((0,))
           Traceback (most recent call last):
              ...
           IndexError: _take() must take ids that return a Variable of the same shape
           
        """

        ## Binary operator doesn't work because ids is turned into a _Constant Variable
        ## which contains floats and not integers. Numeric.take needs integers for ids.
        ## return self._getBinaryOperatorVariable(lambda a, b: numerix.take(a, b, axis = axis), ids) 

        if numerix.take(self.getValue(), ids, axis = axis).shape == self.getShape():
            return self._getUnaryOperatorVariable(lambda a: numerix.take(a, ids, axis = axis), canInline = False)
        else:
            raise IndexError, '_take() must take ids that return a Variable of the same shape'
            
    def allclose(self, other, rtol = 1.e-10, atol = 1.e-10):
        """
           >>> var = Variable((1, 1))
           >>> print var.allclose((1, 1))
           1
           >>> print var.allclose((1,))
           1
           >>> print var.allclose((1,1,1))
           Traceback (most recent call last):
               ...
           ValueError

        The following test is to check that the system does not run
        out of memory.

           >>> from fipy.tools import numerix
           >>> var = Variable(numerix.ones(10000))
           >>> var.allclose(numerix.ones(10001))
           Traceback (most recent call last):
               ...
           ValueError
           
        """

        ## This operation passes `rotateShape = False` to stop the variable being rotated. This
        ## is due to the following strange behaviour in Numeric.allclose. The following code snippet runs
        ## out of memory.
        ##
        ##    >>> from fipy.tools import numerix
        ##    >>> a = numerix.ones(10000)
        ##    >>> b = numeri`x.ones(10001)
        ##    >>> b = b[...,numerix.NewAxis]
        ##    >>> numerix.allclose(a, b)
        ##    Traceback (most recent call last):
        ##    ...
        ##    MemoryError: can't allocate memory for array
        ##
    


        return self._getBinaryOperatorVariable(lambda a,b: numerix.allclose(a, b, atol = atol, rtol = rtol), 
                                               other, 
                                               baseClass = Variable,
                                               opShape = "number",
                                               rotateShape = False,
                                               canInline = False)
        
    def allequal(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.allequal(a,b), 
                                               other,
                                               baseClass = Variable,
                                               opShape = "number",
                                               canInline = False)

    def getMag(self):
        if self.mag is None:
	    self.mag = self.dot(self).sqrt()
	    
	return self.mag
        
    def _testBinOp(self):
        """
            >>> from fipy.variables.cellVariable import CellVariable
            >>> from fipy.variables.faceVariable import FaceVariable
            >>> from fipy.variables.vectorCellVariable import VectorCellVariable
            >>> from fipy.variables.vectorFaceVariable import VectorFaceVariable
            
            >>> from fipy.meshes.grid2D import Grid2D
            >>> mesh = Grid2D(nx = 3)
            
            
        `CellVariable` * CellVariable
        
            >>> cv = CellVariable(mesh = mesh, value = (0, 1, 2))
            >>> cvXcv = cv * cv
            >>> print cvXcv
            [ 0.  1.  4.]
            >>> print isinstance(cvXcv, CellVariable)
            1
        
        `CellVariable` * FaceVariable
        
            >>> fv = FaceVariable(mesh = mesh, value = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
            >>> fvXcv = fv * cv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> cvXfv = cv * fv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            
        `CellVariable` * VectorCellVariable
        
            >>> vcv = VectorCellVariable(mesh = mesh, value = ((0,1),(1,2),(2,3)))
            >>> vcvXcv = vcv * cv
            >>> print vcvXcv
            [[ 0.  0.]
             [ 1.  2.]
             [ 4.  6.]]
            >>> print isinstance(vcvXcv, VectorCellVariable)
            1
            >>> cvXvcv = cv * vcv
            >>> print cvXvcv
            [[ 0.  0.]
             [ 1.  2.]
             [ 4.  6.]]
            >>> print isinstance(cvXvcv, VectorCellVariable)
            1

        `CellVariable` * VectorFaceVariable

            >>> vfv = VectorFaceVariable(mesh = mesh, value = ((0,1),(1,2),(2,3),(3,4),(1,3),(2,4),(3,5),(6,9),(2,6),(1,3)))
            >>> vfvXcv = vfv * cv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> cvXvfv = cv * vfv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        `CellVariable` * Scalar
        
            >>> cvXs = cv * 3
            >>> print cvXs
            [ 0.  3.  6.]
            >>> print isinstance(cvXs, CellVariable)
            1
            >>> sXcv = 3 * cv
            >>> print sXcv
            [ 0.  3.  6.]
            >>> print isinstance(sXcv, CellVariable)
            1

        `CellVariable` * Vector
        
            >>> cvXv2 = cv * (3,2)
            >>> print cvXv2
            [[ 0.  0.]
             [ 3.  2.]
             [ 6.  4.]]
            >>> print isinstance(cvXv2, VectorCellVariable)
            1
            >>> v2Xcv = (3,2) * cv
            >>> print v2Xcv
            [[ 0.  0.]
             [ 3.  2.]
             [ 6.  4.]]
            >>> print isinstance(v2Xcv, VectorCellVariable)
            1
            
            >>> cvXv3 = cv * (3,2,1)
            >>> print cvXv3
            [ 0.  2.  2.]
            >>> print isinstance(cvXv3, CellVariable)
            1
            >>> v3Xcv = (3,2,1) * cv
            >>> print v3Xcv
            [ 0.  2.  2.]
            >>> print isinstance(v3Xcv, CellVariable)
            1
            
            >>> cvXv4 = cv * (3,2,1,0) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int
            >>> v4Xcv = (3,2,1,0) * cv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int


        `CellVariable` * `Variable` Scalar
        
            >>> cvXsv = cv * Variable(value = 3)
            >>> print cvXsv
            [ 0.  3.  6.]
            >>> print isinstance(cvXsv, CellVariable)
            1
            >>> svXcv = Variable(value = 3) * cv
            >>> print svXcv
            [ 0.  3.  6.]
            >>> print isinstance(svXcv, CellVariable)
            1
        
        `binOp` `CellVariable` * `binOp` `Variable` Scalar

            >>> cvcvXsvsv = (cv * cv) * (Variable(value = 3) * Variable(value = 3))
            >>> print cvcvXsvsv
            [  0.   9.  36.]
            >>> print isinstance(cvcvXsvsv, CellVariable)
            1
            >>> svsvXcvcv = (Variable(value = 3) * Variable(value = 3)) * (cv * cv)
            >>> print svsvXcvcv
            [  0.   9.  36.]
            >>> print isinstance(svsvXcvcv, CellVariable)
            1
            
        `CellVariable` * `Variable` Vector
            
            >>> cvXv2v = cv * Variable(value = (3,2))
            >>> print cvXv2v
            [[ 0.  0.]
             [ 3.  2.]
             [ 6.  4.]]
            >>> print isinstance(cvXv2v, VectorCellVariable)
            1
            >>> v2vXcv = Variable(value = (3,2)) * cv
            >>> print v2vXcv
            [[ 0.  0.]
             [ 3.  2.]
             [ 6.  4.]]
            >>> print isinstance(v2vXcv, VectorCellVariable)
            1
            
            >>> cvXv3v = cv * Variable(value = (3,2,1))
            >>> print cvXv3v
            [ 0.  2.  2.]
            >>> print isinstance(cvXv3v, CellVariable)
            1
            >>> v3vXcv = Variable(value = (3,2,1)) * cv
            >>> print v3vXcv
            [ 0.  2.  2.]
            >>> print isinstance(v3vXcv, CellVariable)
            1

            >>> cvXv4v = cv * Variable(value = (3,2,1,0)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v4vXcv = Variable(value = (3,2,1,0)) * cv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            

        `CellVariable` * CellGradVariable
        
            >>> cvXcgv = cv * cv.getGrad()
            >>> print cvXcgv
            [[ 0.  0.]
             [ 1.  0.]
             [ 1.  0.]]
            >>> print isinstance(cvXcgv, VectorCellVariable)
            1
            
        `FaceVariable` * FaceVariable

            >>> fvXfv = fv * fv
            >>> print fvXfv
            [  0.   1.   4.   9.  16.  25.  36.  49.  64.  81.]
            >>> print isinstance(fvXfv, FaceVariable)
            1

        `FaceVariable` * VectorCellVariable

            >>> vcvXfv = vcv * fv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> fvXvcv = fv * vcv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        `FaceVariable` * VectorFaceVariable

            >>> vfvXfv = vfv * fv
            >>> print vfvXfv
            [[  0.   0.]
             [  1.   2.]
             [  4.   6.]
             [  9.  12.]
             [  4.  12.]
             [ 10.  20.]
             [ 18.  30.]
             [ 42.  63.]
             [ 16.  48.]
             [  9.  27.]]
            >>> print isinstance(vfvXfv, VectorFaceVariable)
            1
            >>> fvXvfv = fv * vfv
            >>> print fvXvfv
            [[  0.   0.]
             [  1.   2.]
             [  4.   6.]
             [  9.  12.]
             [  4.  12.]
             [ 10.  20.]
             [ 18.  30.]
             [ 42.  63.]
             [ 16.  48.]
             [  9.  27.]]
            >>> print isinstance(fvXvfv, VectorFaceVariable)
            1

        `FaceVariable` * Scalar

            >>> fvXs = fv * 3
            >>> print fvXs
            [  0.   3.   6.   9.  12.  15.  18.  21.  24.  27.]
            >>> print isinstance(fvXs, FaceVariable)
            1
            >>> sXfv = 3 * fv
            >>> print sXfv
            [  0.   3.   6.   9.  12.  15.  18.  21.  24.  27.]
            >>> print isinstance(sXfv, FaceVariable)
            1

        `FaceVariable` * Vector

            >>> fvXv2 = fv * (3,2)
            >>> print fvXv2
            [[  0.   0.]
             [  3.   2.]
             [  6.   4.]
             [  9.   6.]
             [ 12.   8.]
             [ 15.  10.]
             [ 18.  12.]
             [ 21.  14.]
             [ 24.  16.]
             [ 27.  18.]]
            >>> print isinstance(fvXv2, VectorFaceVariable)
            1
            >>> v2Xfv = (3,2) * fv
            >>> print v2Xfv
            [[  0.   0.]
             [  3.   2.]
             [  6.   4.]
             [  9.   6.]
             [ 12.   8.]
             [ 15.  10.]
             [ 18.  12.]
             [ 21.  14.]
             [ 24.  16.]
             [ 27.  18.]]
            >>> print isinstance(v2Xfv, VectorFaceVariable)
            1
            
            >>> fvXv3 = fv * (3,2,1) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3Xfv = (3,2,1) * fv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

            >>> fvXv10 = fv * (9,8,7,6,5,4,3,2,1,0)
            >>> print fvXv10
            [  0.   8.  14.  18.  20.  20.  18.  14.   8.   0.]
            >>> print isinstance(fvXv10, FaceVariable)
            1
            >>> v10Xfv = (9,8,7,6,5,4,3,2,1,0) * fv
            >>> print v10Xfv
            [  0.   8.  14.  18.  20.  20.  18.  14.   8.   0.]
            >>> print isinstance(v10Xfv, FaceVariable)
            1

        `FaceVariable` * `Variable` Scalar

            >>> fvXsv = fv * Variable(value = 3)
            >>> print fvXsv
            [  0.   3.   6.   9.  12.  15.  18.  21.  24.  27.]
            >>> print isinstance(fvXsv, FaceVariable)
            1
            >>> svXfv = Variable(value = 3) * fv
            >>> print svXfv
            [  0.   3.   6.   9.  12.  15.  18.  21.  24.  27.]
            >>> print isinstance(svXfv, FaceVariable)
            1

        `FaceVariable` * `Variable` Vector
            
            >>> fvXv2v = fv * Variable(value = (3,2))
            >>> print fvXv2v
            [[  0.   0.]
             [  3.   2.]
             [  6.   4.]
             [  9.   6.]
             [ 12.   8.]
             [ 15.  10.]
             [ 18.  12.]
             [ 21.  14.]
             [ 24.  16.]
             [ 27.  18.]]
            >>> print isinstance(fvXv2v, VectorFaceVariable)
            1
            >>> v2vXfv = Variable(value = (3,2)) * fv
            >>> print v2vXfv
            [[  0.   0.]
             [  3.   2.]
             [  6.   4.]
             [  9.   6.]
             [ 12.   8.]
             [ 15.  10.]
             [ 18.  12.]
             [ 21.  14.]
             [ 24.  16.]
             [ 27.  18.]]
            >>> print isinstance(v2vXfv, VectorFaceVariable)
            1
            
            >>> fvXv3v = fv * Variable(value = (3,2,1)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3vXfv = Variable(value = (3,2,1)) * fv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

            >>> fvXv10v = fv * Variable(value = (9,8,7,6,5,4,3,2,1,0))
            >>> print fvXv10v
            [  0.   8.  14.  18.  20.  20.  18.  14.   8.   0.]
            >>> print isinstance(fvXv10v, FaceVariable)
            1
            >>> v10vXfv = Variable(value = (9,8,7,6,5,4,3,2,1,0)) * fv
            >>> print v10vXfv
            [  0.   8.  14.  18.  20.  20.  18.  14.   8.   0.]
            >>> print isinstance(v10vXfv, FaceVariable)
            1

            
            
        `VectorCellVariable` * VectorCellVariable

            >>> vcvXvcv = vcv * vcv
            >>> print vcvXvcv
            [[ 0.  1.]
             [ 1.  4.]
             [ 4.  9.]]
            >>> print isinstance(vcvXvcv, VectorCellVariable)
            1

        `VectorCellVariable` * VectorFaceVariable

            >>> vfvXvcv = vfv * vcv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> vcvXvfv = vcv * vfv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        `VectorCellVariable` * Scalar

            >>> vcvXs = vcv * 3
            >>> print vcvXs
            [[ 0.  3.]
             [ 3.  6.]
             [ 6.  9.]]
            >>> print isinstance(vcvXs, VectorCellVariable)
            1
            >>> sXvcv = 3 * vcv
            >>> print sXvcv
            [[ 0.  3.]
             [ 3.  6.]
             [ 6.  9.]]
            >>> print isinstance(vcvXs, VectorCellVariable)
            1

        `VectorCellVariable` * Vector

            >>> vcvXv2 = vcv * (3,2)
            >>> print vcvXv2
            [[ 0.  2.]
             [ 3.  4.]
             [ 6.  6.]]
            >>> print isinstance(vcvXv2, VectorCellVariable)
            1
            >>> v2Xvcv = (3,2) * vcv
            >>> print v2Xvcv
            [[ 0.  2.]
             [ 3.  4.]
             [ 6.  6.]]
            >>> print isinstance(v2Xvcv, VectorCellVariable)
            1
            
            >>> vcvXv3 = vcv * (3,2,1)
            >>> print vcvXv3
            [[ 0.  3.]
             [ 2.  4.]
             [ 2.  3.]]
            >>> isinstance(vcvXv3, VectorCellVariable)
            1
            >>> v3Xvcv = (3,2,1) * vcv 
            >>> print v3Xvcv
            [[ 0.  3.]
             [ 2.  4.]
             [ 2.  3.]]
            >>> isinstance(v3Xvcv, VectorCellVariable)
            1

            >>> vcvXv4 = vcv * (3,2,1,0) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v4Xvcv = (3,2,1,0) * vcv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        `VectorCellVariable` * `Variable` Scalar

            >>> vcvXsv = vcv * Variable(value = 3)
            >>> print vcvXsv
            [[ 0.  3.]
             [ 3.  6.]
             [ 6.  9.]]
            >>> print isinstance(vcvXsv, VectorCellVariable)
            1
            >>> svXvcv = Variable(value = 3) * vcv
            >>> print svXvcv
            [[ 0.  3.]
             [ 3.  6.]
             [ 6.  9.]]
            >>> print isinstance(svXvcv, VectorCellVariable)
            1

        `VectorCellVariable` * `Variable` Vector
            
            >>> vcvXv2v = vcv * Variable(value = (3,2))
            >>> print vcvXv2v
            [[ 0.  2.]
             [ 3.  4.]
             [ 6.  6.]]
            >>> print isinstance(vcvXv2v, VectorCellVariable)
            1
            >>> v2vXvcv = Variable(value = (3,2)) * vcv
            >>> print v2vXvcv
            [[ 0.  2.]
             [ 3.  4.]
             [ 6.  6.]]
            >>> print isinstance(v2vXvcv, VectorCellVariable)
            1
            
            >>> vcvXv3v = vcv * Variable(value = (3,2,1))
            >>> print vcvXv3v
            [[ 0.  3.]
             [ 2.  4.]
             [ 2.  3.]]
            >>> isinstance(vcvXv3v, VectorCellVariable)
            1
            >>> v3vXvcv = Variable(value = (3,2,1)) * vcv 
            >>> print v3vXvcv
            [[ 0.  3.]
             [ 2.  4.]
             [ 2.  3.]]
            >>> isinstance(v3vXvcv, VectorCellVariable)
            1

            >>> vcvXv4v = vcv * Variable(value = (3,2,1,0)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v4vXvcv = Variable(value = (3,2,1,0)) * vcv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

                        
        `VectorFaceVariable` * VectorFaceVariable

            >>> vfvXvfv = vfv * vfv
            >>> print vfvXvfv
            [[  0.   1.]
             [  1.   4.]
             [  4.   9.]
             [  9.  16.]
             [  1.   9.]
             [  4.  16.]
             [  9.  25.]
             [ 36.  81.]
             [  4.  36.]
             [  1.   9.]]
            >>> isinstance(vfvXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * Scalar

            >>> vfvXs = vfv * 3
            >>> print vfvXs
            [[  0.   3.]
             [  3.   6.]
             [  6.   9.]
             [  9.  12.]
             [  3.   9.]
             [  6.  12.]
             [  9.  15.]
             [ 18.  27.]
             [  6.  18.]
             [  3.   9.]]
            >>> print isinstance(vfvXs, VectorFaceVariable)
            1
            >>> sXvfv = 3 * vfv
            >>> print sXvfv
            [[  0.   3.]
             [  3.   6.]
             [  6.   9.]
             [  9.  12.]
             [  3.   9.]
             [  6.  12.]
             [  9.  15.]
             [ 18.  27.]
             [  6.  18.]
             [  3.   9.]]
            >>> print isinstance(sXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * Vector

            >>> vfvXv2 = vfv * (3,2)
            >>> print vfvXv2
            [[  0.   2.]
             [  3.   4.]
             [  6.   6.]
             [  9.   8.]
             [  3.   6.]
             [  6.   8.]
             [  9.  10.]
             [ 18.  18.]
             [  6.  12.]
             [  3.   6.]]
            >>> print isinstance(vfvXv2, VectorFaceVariable)
            1
            >>> v2Xvfv = (3,2) * vfv
            >>> print v2Xvfv
            [[  0.   2.]
             [  3.   4.]
             [  6.   6.]
             [  9.   8.]
             [  3.   6.]
             [  6.   8.]
             [  9.  10.]
             [ 18.  18.]
             [  6.  12.]
             [  3.   6.]]
            >>> print isinstance(v2Xvfv, VectorFaceVariable)
            1
            
            >>> vfvXv3 = vfv * (2,1,0) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3Xvfv = (2,1,0) * vfv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int


            >>> vfvXv10 = vfv * (9,8,7,6,5,4,3,2,1,0)
            >>> print vfvXv10
            [[  0.   9.]
             [  8.  16.]
             [ 14.  21.]
             [ 18.  24.]
             [  5.  15.]
             [  8.  16.]
             [  9.  15.]
             [ 12.  18.]
             [  2.   6.]
             [  0.   0.]]
            >>> isinstance(vfvXv10, VectorFaceVariable)
            1
            >>> v10Xvfv = (9,8,7,6,5,4,3,2,1,0) * vfv
            >>> print v10Xvfv
            [[  0.   9.]
             [  8.  16.]
             [ 14.  21.]
             [ 18.  24.]
             [  5.  15.]
             [  8.  16.]
             [  9.  15.]
             [ 12.  18.]
             [  2.   6.]
             [  0.   0.]]
            >>> isinstance(v10Xvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * `Variable` Scalar

            >>> vfvXsv = vfv * Variable(value = 3)
            >>> print vfvXsv
            [[  0.   3.]
             [  3.   6.]
             [  6.   9.]
             [  9.  12.]
             [  3.   9.]
             [  6.  12.]
             [  9.  15.]
             [ 18.  27.]
             [  6.  18.]
             [  3.   9.]]
            >>> print isinstance(vfvXsv, VectorFaceVariable)
            1
            >>> svXvfv = Variable(value = 3) * vfv
            >>> print svXvfv
            [[  0.   3.]
             [  3.   6.]
             [  6.   9.]
             [  9.  12.]
             [  3.   9.]
             [  6.  12.]
             [  9.  15.]
             [ 18.  27.]
             [  6.  18.]
             [  3.   9.]]
            >>> print isinstance(svXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * `Variable` Vector
            
            >>> vfvXv2v = vfv * Variable(value = (3,2))
            >>> print vfvXv2v
            [[  0.   2.]
             [  3.   4.]
             [  6.   6.]
             [  9.   8.]
             [  3.   6.]
             [  6.   8.]
             [  9.  10.]
             [ 18.  18.]
             [  6.  12.]
             [  3.   6.]]
            >>> print isinstance(vfvXv2v, VectorFaceVariable)
            1
            >>> v2vXvfv = Variable(value = (3,2)) * vfv
            >>> print v2vXvfv
            [[  0.   2.]
             [  3.   4.]
             [  6.   6.]
             [  9.   8.]
             [  3.   6.]
             [  6.   8.]
             [  9.  10.]
             [ 18.  18.]
             [  6.  12.]
             [  3.   6.]]
            >>> print isinstance(v2vXvfv, VectorFaceVariable)
            1
            
            >>> vfvXv3v = vfv * Variable(value = (2,1,0)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3vXvfv = Variable(value = (2,1,0)) * vfv #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int


            >>> vfvXv10v = vfv * Variable(value = (9,8,7,6,5,4,3,2,1,0))
            >>> print vfvXv10v
            [[  0.   9.]
             [  8.  16.]
             [ 14.  21.]
             [ 18.  24.]
             [  5.  15.]
             [  8.  16.]
             [  9.  15.]
             [ 12.  18.]
             [  2.   6.]
             [  0.   0.]]
            >>> isinstance(vfvXv10v, VectorFaceVariable)
            1
            >>> v10vXvfv = Variable(value = (9,8,7,6,5,4,3,2,1,0)) * vfv
            >>> print v10vXvfv
            [[  0.   9.]
             [  8.  16.]
             [ 14.  21.]
             [ 18.  24.]
             [  5.  15.]
             [  8.  16.]
             [  9.  15.]
             [ 12.  18.]
             [  2.   6.]
             [  0.   0.]]
            >>> isinstance(v10vXvfv, VectorFaceVariable)
            1

            
            
        Scalar * `Variable` Scalar

            >>> sXsv = 3 * Variable(value = 3)
            >>> print sXsv
            9
            >>> print isinstance(sXsv, Variable)
            1
            >>> svXs = Variable(value = 3) * 3
            >>> print svXs
            9
            >>> print isinstance(svXs, Variable)
            1

        Scalar * `Variable` Vector
            
            >>> sXv2v = 3 * Variable(value = (3,2))
            >>> print sXv2v
            [ 9.  6.]
            >>> print isinstance(sXv2v, Variable)
            1
            >>> v2vXs = Variable(value = (3,2)) * 3
            >>> print v2vXs
            [ 9.  6.]
            >>> print isinstance(v2vXs, Variable)
            1
            
            
            
        Vector * `Variable` Scalar

            >>> vXsv = (3, 2) * Variable(value = 3)
            >>> print vXsv
            [ 9.  6.]
            >>> print isinstance(vXsv, Variable)
            1
            >>> svXv = Variable(value = 3) * (3, 2)
            >>> print svXv
            [ 9.  6.]
            >>> print isinstance(svXv, Variable)
            1

        Vector * `Variable` Vector
            
            >>> vXv2v = (3, 2) * Variable(value = (3,2))
            >>> print vXv2v
            [ 9.  4.]
            >>> print isinstance(vXv2v, Variable)
            1
            >>> v2vXv = Variable(value = (3,2)) * (3, 2)
            >>> print v2vXv
            [ 9.  4.]
            >>> print isinstance(v2vXv, Variable)
            1

            >>> vXv3v = (3, 2, 1) * Variable(value = (3,2)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3vXv = Variable(value = (3,2)) * (3, 2, 1) #doctest: +IGNORE_EXCEPTION_DETAIL 
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            

        `Variable` Scalar * `Variable` Scalar

            >>> svXsv = Variable(value = 3) * Variable(value = 3)
            >>> print svXsv
            9
            >>> print isinstance(svXsv, Variable)
            1

        `Variable` Scalar * `Variable` Vector
            
            >>> svXv2v = Variable(value = 3) * Variable(value = (3,2))
            >>> print svXv2v
            [ 9.  6.]
            >>> print isinstance(svXv2v, Variable)
            1
            >>> v2vXsv = Variable(value = (3,2)) * Variable(value = 3)
            >>> print v2vXsv
            [ 9.  6.]
            >>> print isinstance(v2vXsv, Variable)
            1

            
        `Variable` Vector * `Variable` Vector
            
            >>> v2vXv2v = Variable(value = (3, 2)) * Variable(value = (3,2))
            >>> print v2vXv2v
            [ 9.  4.]
            >>> print isinstance(v2vXv2v, Variable)
            1
            
            >>> v3vXv2v = Variable(value = (3, 2, 1)) * Variable(value = (3,2)) #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        Test for weird bug that was appearing in inline. Caused by the intermediate
        operators not getting marked fresh.

            >>> class Alpha(Variable):
            ...     def __init__(self, var):
            ...         Variable.__init__(self)
            ...         self.var = self._requires(var)
            ...     def _calcValue(self):
            ...         return self.var.getValue()

            >>> coeff = Variable()
            >>> alpha = Alpha(-coeff / 1)
            >>> alpha.getValue()
            -0.0
            >>> coeff.setValue(-10.0)
            >>> alpha.getValue()
            10.0
            >>> coeff.setValue(10.0)
            >>> alpha.getValue()
            -10.0

        Test to prevent divide by zero evaluation before value is
        requested.  The request is caused by the Variable requiring
        its unit to see whether it can do an inline calculation in
        _getUnaryOperatorVariable().
        
            >>> T = Variable()
            >>> from fipy import numerix
            >>> v = numerix.exp(-T / (1. *  T))
        """
        pass

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
