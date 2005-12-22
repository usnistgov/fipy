#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "variable.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 12/22/05 {4:26:51 PM} 
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


import Numeric

import fipy.tools.dimensions.physicalField

from fipy.tools import numerix

class Variable:
    
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
    
    def __init__(self, value=0., unit = None, array = None, name = '', mesh = None):
	"""
	Create a `Variable`.
	
	    >>> Variable(value = 3)
	    Variable(value = 3)
	    >>> Variable(value = 3, unit = "m")
	    Variable(value = PhysicalField(3,'m'))
	    >>> Variable(value = 3, unit = "m", array = Numeric.zeros((3,2)))
	    Variable(value = PhysicalField([[3,3,]
	     [3,3,]
	     [3,3,]],'m'))

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
	    unit = None
	    array = None
	    
	self._setValue(value = value, unit = unit, array = array)
	
	self.name = name
	self.mesh = mesh
		
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
        Attempt to convert the `Variable` to a Numeric `array` object
        
            >>> v = Variable(value = [2,3])
            >>> Numeric.array(v)
            [ 2., 3.,]
            
        It is an error to convert a dimensional `Variable` to a 
        Numeric `array`
        
            >>> v = Variable(value = [2,3], unit = "m")
            >>> Numeric.array(v)
            Traceback (most recent call last):
                ...
            TypeError: Numeric array value must be dimensionless

	"""
	return Numeric.array(self.getValue(), t)
	
    def copy(self):
	"""
	Make an duplicate of the `Variable`
	
	    >>> a = Variable(value = 3)
	    >>> b = a.copy()
	    >>> b
	    Variable(value = 3)

	The duplicate will not reflect changes made to the original
	
	    >>> a.setValue(5)
	    >>> b
	    Variable(value = 3)
	"""
	return Variable(value = self)
	
    def getUnit(self):
	"""
	Return the unit object of `self`.
	
	    >>> Variable(value = "1 m").getUnit()
	    <PhysicalUnit m>
	"""
	value = self.getValue()
	if isinstance(value, fipy.tools.dimensions.physicalField.PhysicalField):
	    return value.getUnit()
	else:
	    return "1"
	
    def inBaseUnits(self):
	"""
	Return the value of the `Variable` with all units reduced to 
	their base SI elements.
	
	    >>> e = Variable(value = "2.7 Hartree*Nav")
	    >>> print e.inBaseUnits()
	    7088849.01085 kg*m**2/s**2/mol
	"""
	value = self.getValue()
	if isinstance(value, fipy.tools.dimensions.physicalField.PhysicalField):
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
        if isinstance(value, fipy.tools.dimensions.physicalField.PhysicalField):
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
            IndexError: index out of bounds

        """
	return self.getValue()[index]

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
	
    def tostring(self, max_line_width = None, precision = None, suppress_small = None, separator = ' '):
        return numerix.tostring(self.getValue(), 
                                max_line_width = max_line_width,
                                precision = precision, 
                                suppress_small = suppress_small, 
                                separator = separator)
        
    def __setitem__(self, index, value):
	self.value[index] = value
	self._markFresh()
	
    def __call__(self):
	"""
	"Evaluate" the `Variable` and return its value
	
	    >>> a = Variable(value = 3)
	    >>> a()
	    3
	    >>> b = a + 4
	    >>> b
	    (Variable(value = 3) + 4)
	    >>> b()
	    7
	"""
	return self.getValue()
		
    def getValue(self):
	"""
	"Evaluate" the `Variable` and return its value (longhand)
	
	    >>> a = Variable(value = 3)
	    >>> a.getValue()
	    3
	    >>> b = a + 4
	    >>> b
	    (Variable(value = 3) + 4)
	    >>> b.getValue()
	    7
	"""
        if self.stale:           
            self.value = self._calcValue()
            self._markFresh()

	return self.value
 
    def _setValue(self, value, unit = None, array = None):
        self.value = self._makeValue(value = value, unit = unit, array = array)
     
    def _makeValue(self, value, unit = None, array = None):
        PF = fipy.tools.dimensions.physicalField.PhysicalField

        from fipy.tools.numerix import MA
        if not isinstance(value, PF):
            if unit is not None or type(value) in [type(''), type(()), type([])]:
                value = PF(value = value, unit = unit, array = array)
            elif array is not None:
                array[:] = value
                value = array
            elif type(value) not in (type(Numeric.array(1)), type(MA.array(1))):
                value = Numeric.array(value)
            
        if isinstance(value, PF) and value.getUnit().isDimensionless():
            value = value.getNumericValue()
            
        return value

    def setValue(self, value, unit = None, array = None):
	self._setValue(value = value, unit = unit, array = array)
	self._markFresh()
	
    def _setNumericValue(self, value):
	if isinstance(self.value, fipy.tools.dimensions.physicalField.PhysicalField):
	    self.value.value = value
	else:
	    self.value = value
	
    def _getArray(self):
	if isinstance(self.value, fipy.tools.dimensions.physicalField.PhysicalField):
	    return self.value._getArray()
	else:
	    return self.value
     
    def getNumericValue(self):
	value = self.getValue()
	if isinstance(value, fipy.tools.dimensions.physicalField.PhysicalField):
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
        return numerix.array(self._getArray()).shape
	
    def _calcValue(self):
	pass	
 
    def __markStale(self):
        import weakref
        remainingSubscribedVariables = []
        for subscriber in self.subscribedVariables:
            try:
                subscriber._markStale() 
                remainingSubscribedVariables.append(subscriber)
            except weakref.ReferenceError:
                pass
        self.subscribedVariables = remainingSubscribedVariables

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
	self.subscribedVariables.append(weakref.proxy(var))
	
    def _getVariableClass(self):
	return Variable

    def _getOperatorVariableClass(self, baseClass = None):
	if baseClass is None:
            baseClass = self._getVariableClass()
            
	class OperatorVariable(baseClass):
	    def __init__(self, op, var, mesh = None):
                mesh = mesh or var[0].getMesh() or (len(var) > 1 and var[1].getMesh())
		self.op = op
		self.var = var
                baseClass.__init__(self, value = 0, mesh = mesh)
                self.name = ''
		for aVar in self.var:
		    self._requires(aVar)

                self.old = None

            def getOld(self):
                if self.old is None:
                    oldVar = []
                    for v in self.var:
                        from fipy.variables.cellVariable import CellVariable
                        if isinstance(v, CellVariable):
                            oldVar.append(v.getOld())
                        else:
                            oldVar.append(v)
                    self.old = self.__class__(self.op, oldVar, self.getMesh())

                return self.old
            
	    def _getRepresentation(self, style = "__repr__"):
                """
                :Parameters:
                    
                  - `style`: one of `'__repr__'`, `'name'`, `'TeX'`
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
			return stack.pop()
                    elif opcode.opname[bytecode] == 'LOAD_CONST':
                        stack.append(self.op.func_code.co_consts[_popIndex()])
		    elif opcode.opname[bytecode] == 'LOAD_ATTR':
			stack.append(stack.pop() + "." + self.op.func_code.co_names[_popIndex()])
		    elif opcode.opname[bytecode] == 'COMPARE_OP':
			stack.append(stack.pop(-2) + " " + opcode.cmp_op[_popIndex()] + " " + stack.pop())
		    elif opcode.opname[bytecode] == 'LOAD_GLOBAL':
			stack.append(self.op.func_code.co_names[_popIndex()])
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
			raise SyntaxError, "Unknown bytecode: %s in %s: %s" % (`bytecode`, `[ord(byte) for byte in self.op.func_code.co_code]`)

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
		   mesh = self.getMesh())

	return OperatorVariable
	
    def _getArithmeticBaseClass(self, other = None):
        """
        Given `self` and `other`, return the desired base
        class for an operation result.
        """
        if other is None:
            return Variable
            
        if other.__class__.__name__ is self.__class__.__name__ \
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
            
    def _getUnaryOperatorVariable(self, op, baseClass = None):
	class unOp(self._getOperatorVariableClass(baseClass)):
	    def _calcValue(self):
		return self.op(self.var[0].getValue())
		
	return unOp(op, [self])
	    
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
        
        a = object._getArray()
        
        # we don't want to meddle with the contents of the actual object
        if type(a) in (type(()), type([]), type(numerix.array(1))):
            a = a.copy()
        else:
            a = numerix.array(1)
            
        if object not in valueMattersForShape:
            if a.shape == ():
                # if Numeric thinks the array is a scalar (rather than a 1x1 array)
                # it won't slice
                a = numerix.array(1)
            else:
                a[:] = 1
            
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
    
    def _verifyShape(self, op, var0, var1, var0Array, var1Array, opShape, otherClass):
        try:
            # check if applying the operation to the inputs will produce the desired shape
            if numerix.getShape(op(var0Array, var1Array)) != opShape:
                raise ValueError
        except ValueError:
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
                    
        return (var0, var1)

    def _getBinaryOperatorVariable(self, op, other, baseClass = None, opShape = None, valueMattersForShape = ()):
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
            [ 0., 1., 4.,]
            >>> print isinstance(cvXcv, CellVariable)
            1
        
        `CellVariable` * FaceVariable
        
            >>> fv = FaceVariable(mesh = mesh, value = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
            >>> fvXcv = fv * cv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> cvXfv = cv * fv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

        `CellVariable` * VectorCellVariable
        
            >>> vcv = VectorCellVariable(mesh = mesh, value = ((0,1),(1,2),(2,3)))
            >>> vcvXcv = vcv * cv
            >>> print vcvXcv
            [[ 0., 0.,]
             [ 1., 2.,]
             [ 4., 6.,]]
            >>> print isinstance(vcvXcv, VectorCellVariable)
            1
            >>> cvXvcv = cv * vcv
            >>> print cvXvcv
            [[ 0., 0.,]
             [ 1., 2.,]
             [ 4., 6.,]]
            >>> print isinstance(cvXvcv, VectorCellVariable)
            1

        `CellVariable` * VectorFaceVariable

            >>> vfv = VectorFaceVariable(mesh = mesh, value = ((0,1),(1,2),(2,3),(3,4),(1,3),(2,4),(3,5),(6,9),(2,6),(1,3)))
            >>> vfvXcv = vfv * cv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> cvXvfv = cv * vfv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

        `CellVariable` * Scalar
        
            >>> cvXs = cv * 3
            >>> print cvXs
            [ 0., 3., 6.,]
            >>> print isinstance(cvXs, CellVariable)
            1
            >>> sXcv = 3 * cv
            >>> print sXcv
            [ 0., 3., 6.,]
            >>> print isinstance(sXcv, CellVariable)
            1

        `CellVariable` * Vector
        
            >>> cvXv2 = cv * (3,2)
            >>> print cvXv2
            [[ 0., 0.,]
             [ 3., 2.,]
             [ 6., 4.,]]
            >>> print isinstance(cvXv2, VectorCellVariable)
            1
            >>> v2Xcv = (3,2) * cv
            >>> print v2Xcv
            [[ 0., 0.,]
             [ 3., 2.,]
             [ 6., 4.,]]
            >>> print isinstance(v2Xcv, VectorCellVariable)
            1
            
            >>> cvXv3 = cv * (3,2,1)
            >>> print cvXv3
            [ 0., 2., 2.,]
            >>> print isinstance(cvXv3, CellVariable)
            1
            >>> v3Xcv = (3,2,1) * cv
            >>> print v3Xcv
            [ 0., 2., 2.,]
            >>> print isinstance(v3Xcv, CellVariable)
            1
            
            >>> cvXv4 = cv * (3,2,1,0)
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int
            >>> v4Xcv = (3,2,1,0) * cv
            Traceback (most recent call last):
                ...
            TypeError: can't multiply sequence to non-int


        `CellVariable` * `Variable` Scalar
        
            >>> cvXsv = cv * Variable(value = 3)
            >>> print cvXsv
            [ 0., 3., 6.,]
            >>> print isinstance(cvXsv, CellVariable)
            1
            >>> svXcv = Variable(value = 3) * cv
            >>> print svXcv
            [ 0., 3., 6.,]
            >>> print isinstance(svXcv, CellVariable)
            1
        
        `CellVariable` * `Variable` Vector
            
            >>> cvXv2v = cv * Variable(value = (3,2))
            >>> print cvXv2v
            [[ 0., 0.,]
             [ 3., 2.,]
             [ 6., 4.,]]
            >>> print isinstance(cvXv2v, VectorCellVariable)
            1
            >>> v2vXcv = Variable(value = (3,2)) * cv
            >>> print v2vXcv
            [[ 0., 0.,]
             [ 3., 2.,]
             [ 6., 4.,]]
            >>> print isinstance(v2vXcv, VectorCellVariable)
            1
            
            >>> cvXv3v = cv * Variable(value = (3,2,1))
            >>> print cvXv3v
            [ 0., 2., 2.,]
            >>> print isinstance(cvXv3v, CellVariable)
            1
            >>> v3vXcv = Variable(value = (3,2,1)) * cv
            >>> print v3vXcv
            [ 0., 2., 2.,]
            >>> print isinstance(v3vXcv, CellVariable)
            1

            >>> cvXv4v = cv * Variable(value = (3,2,1,0))
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> v4vXcv = Variable(value = (3,2,1,0)) * cv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            

        `CellVariable` * CellGradVariable
        
            >>> cvXcgv = cv * cv.getGrad()
            >>> print cvXcgv
            [[ 0., 0.,]
             [ 1., 0.,]
             [ 1., 0.,]]
            >>> print isinstance(cvXcgv, VectorCellVariable)
            1
            
        `FaceVariable` * FaceVariable

            >>> fvXfv = fv * fv
            >>> print fvXfv
            [  0.,  1.,  4.,  9., 16., 25., 36., 49., 64., 81.,]
            >>> print isinstance(fvXfv, FaceVariable)
            1

        `FaceVariable` * VectorCellVariable

            >>> vcvXfv = vcv * fv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> fvXvcv = fv * vcv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

        `FaceVariable` * VectorFaceVariable

            >>> vfvXfv = vfv * fv
            >>> print vfvXfv
            [[  0.,  0.,]
             [  1.,  2.,]
             [  4.,  6.,]
             [  9., 12.,]
             [  4., 12.,]
             [ 10., 20.,]
             [ 18., 30.,]
             [ 42., 63.,]
             [ 16., 48.,]
             [  9., 27.,]]
            >>> print isinstance(vfvXfv, VectorFaceVariable)
            1
            >>> fvXvfv = fv * vfv
            >>> print fvXvfv
            [[  0.,  0.,]
             [  1.,  2.,]
             [  4.,  6.,]
             [  9., 12.,]
             [  4., 12.,]
             [ 10., 20.,]
             [ 18., 30.,]
             [ 42., 63.,]
             [ 16., 48.,]
             [  9., 27.,]]
            >>> print isinstance(fvXvfv, VectorFaceVariable)
            1

        `FaceVariable` * Scalar

            >>> fvXs = fv * 3
            >>> print fvXs
            [  0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.,]
            >>> print isinstance(fvXs, FaceVariable)
            1
            >>> sXfv = 3 * fv
            >>> print sXfv
            [  0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.,]
            >>> print isinstance(sXfv, FaceVariable)
            1

        `FaceVariable` * Vector

            >>> fvXv2 = fv * (3,2)
            >>> print fvXv2
            [[  0.,  0.,]
             [  3.,  2.,]
             [  6.,  4.,]
             [  9.,  6.,]
             [ 12.,  8.,]
             [ 15., 10.,]
             [ 18., 12.,]
             [ 21., 14.,]
             [ 24., 16.,]
             [ 27., 18.,]]
            >>> print isinstance(fvXv2, VectorFaceVariable)
            1
            >>> v2Xfv = (3,2) * fv
            >>> print v2Xfv
            [[  0.,  0.,]
             [  3.,  2.,]
             [  6.,  4.,]
             [  9.,  6.,]
             [ 12.,  8.,]
             [ 15., 10.,]
             [ 18., 12.,]
             [ 21., 14.,]
             [ 24., 16.,]
             [ 27., 18.,]]
            >>> print isinstance(v2Xfv, VectorFaceVariable)
            1
            
            >>> fvXv3 = fv * (3,2,1)
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3Xfv = (3,2,1) * fv 
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

            >>> fvXv10 = fv * (9,8,7,6,5,4,3,2,1,0)
            >>> print fvXv10
            [  0.,  8., 14., 18., 20., 20., 18., 14.,  8.,  0.,]
            >>> print isinstance(fvXv10, FaceVariable)
            1
            >>> v10Xfv = (9,8,7,6,5,4,3,2,1,0) * fv
            >>> print v10Xfv
            [  0.,  8., 14., 18., 20., 20., 18., 14.,  8.,  0.,]
            >>> print isinstance(v10Xfv, FaceVariable)
            1

        `FaceVariable` * `Variable` Scalar

            >>> fvXsv = fv * Variable(value = 3)
            >>> print fvXsv
            [  0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.,]
            >>> print isinstance(fvXsv, FaceVariable)
            1
            >>> svXfv = Variable(value = 3) * fv
            >>> print svXfv
            [  0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27.,]
            >>> print isinstance(svXfv, FaceVariable)
            1

        `FaceVariable` * `Variable` Vector
            
            >>> fvXv2v = fv * Variable(value = (3,2))
            >>> print fvXv2v
            [[  0.,  0.,]
             [  3.,  2.,]
             [  6.,  4.,]
             [  9.,  6.,]
             [ 12.,  8.,]
             [ 15., 10.,]
             [ 18., 12.,]
             [ 21., 14.,]
             [ 24., 16.,]
             [ 27., 18.,]]
            >>> print isinstance(fvXv2v, VectorFaceVariable)
            1
            >>> v2vXfv = Variable(value = (3,2)) * fv
            >>> print v2vXfv
            [[  0.,  0.,]
             [  3.,  2.,]
             [  6.,  4.,]
             [  9.,  6.,]
             [ 12.,  8.,]
             [ 15., 10.,]
             [ 18., 12.,]
             [ 21., 14.,]
             [ 24., 16.,]
             [ 27., 18.,]]
            >>> print isinstance(v2vXfv, VectorFaceVariable)
            1
            
            >>> fvXv3v = fv * Variable(value = (3,2,1))
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> v3vXfv = Variable(value = (3,2,1)) * fv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

            >>> fvXv10v = fv * Variable(value = (9,8,7,6,5,4,3,2,1,0))
            >>> print fvXv10v
            [  0.,  8., 14., 18., 20., 20., 18., 14.,  8.,  0.,]
            >>> print isinstance(fvXv10v, FaceVariable)
            1
            >>> v10vXfv = Variable(value = (9,8,7,6,5,4,3,2,1,0)) * fv
            >>> print v10vXfv
            [  0.,  8., 14., 18., 20., 20., 18., 14.,  8.,  0.,]
            >>> print isinstance(v10vXfv, FaceVariable)
            1

            
            
        `VectorCellVariable` * VectorCellVariable

            >>> vcvXvcv = vcv * vcv
            >>> print vcvXvcv
            [[ 0., 1.,]
             [ 1., 4.,]
             [ 4., 9.,]]
            >>> print isinstance(vcvXvcv, VectorCellVariable)
            1

        `VectorCellVariable` * VectorFaceVariable

            >>> vfvXvcv = vfv * vcv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> vcvXvfv = vcv * vfv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

        `VectorCellVariable` * Scalar

            >>> vcvXs = vcv * 3
            >>> print vcvXs
            [[ 0., 3.,]
             [ 3., 6.,]
             [ 6., 9.,]]
            >>> print isinstance(vcvXs, VectorCellVariable)
            1
            >>> sXvcv = 3 * vcv
            >>> print sXvcv
            [[ 0., 3.,]
             [ 3., 6.,]
             [ 6., 9.,]]
            >>> print isinstance(vcvXs, VectorCellVariable)
            1

        `VectorCellVariable` * Vector

            >>> vcvXv2 = vcv * (3,2)
            >>> print vcvXv2
            [[ 0., 2.,]
             [ 3., 4.,]
             [ 6., 6.,]]
            >>> print isinstance(vcvXv2, VectorCellVariable)
            1
            >>> v2Xvcv = (3,2) * vcv
            >>> print v2Xvcv
            [[ 0., 2.,]
             [ 3., 4.,]
             [ 6., 6.,]]
            >>> print isinstance(v2Xvcv, VectorCellVariable)
            1
            
            >>> vcvXv3 = vcv * (3,2,1)
            >>> print vcvXv3
            [[ 0., 3.,]
             [ 2., 4.,]
             [ 2., 3.,]]
            >>> isinstance(vcvXv3, VectorCellVariable)
            1
            >>> v3Xvcv = (3,2,1) * vcv 
            >>> print v3Xvcv
            [[ 0., 3.,]
             [ 2., 4.,]
             [ 2., 3.,]]
            >>> isinstance(v3Xvcv, VectorCellVariable)
            1

            >>> vcvXv4 = vcv * (3,2,1,0)
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v4Xvcv = (3,2,1,0) * vcv
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int

        `VectorCellVariable` * `Variable` Scalar

            >>> vcvXsv = vcv * Variable(value = 3)
            >>> print vcvXsv
            [[ 0., 3.,]
             [ 3., 6.,]
             [ 6., 9.,]]
            >>> print isinstance(vcvXsv, VectorCellVariable)
            1
            >>> svXvcv = Variable(value = 3) * vcv
            >>> print svXvcv
            [[ 0., 3.,]
             [ 3., 6.,]
             [ 6., 9.,]]
            >>> print isinstance(svXvcv, VectorCellVariable)
            1

        `VectorCellVariable` * `Variable` Vector
            
            >>> vcvXv2v = vcv * Variable(value = (3,2))
            >>> print vcvXv2v
            [[ 0., 2.,]
             [ 3., 4.,]
             [ 6., 6.,]]
            >>> print isinstance(vcvXv2v, VectorCellVariable)
            1
            >>> v2vXvcv = Variable(value = (3,2)) * vcv
            >>> print v2vXvcv
            [[ 0., 2.,]
             [ 3., 4.,]
             [ 6., 6.,]]
            >>> print isinstance(v2vXvcv, VectorCellVariable)
            1
            
            >>> vcvXv3v = vcv * Variable(value = (3,2,1))
            >>> print vcvXv3v
            [[ 0., 3.,]
             [ 2., 4.,]
             [ 2., 3.,]]
            >>> isinstance(vcvXv3v, VectorCellVariable)
            1
            >>> v3vXvcv = Variable(value = (3,2,1)) * vcv 
            >>> print v3vXvcv
            [[ 0., 3.,]
             [ 2., 4.,]
             [ 2., 3.,]]
            >>> isinstance(v3vXvcv, VectorCellVariable)
            1

            >>> vcvXv4v = vcv * Variable(value = (3,2,1,0))
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> v4vXvcv = Variable(value = (3,2,1,0)) * vcv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'

            
            
            
            
            
        `VectorFaceVariable` * VectorFaceVariable

            >>> vfvXvfv = vfv * vfv
            >>> print vfvXvfv
            [[  0.,  1.,]
             [  1.,  4.,]
             [  4.,  9.,]
             [  9., 16.,]
             [  1.,  9.,]
             [  4., 16.,]
             [  9., 25.,]
             [ 36., 81.,]
             [  4., 36.,]
             [  1.,  9.,]]
            >>> isinstance(vfvXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * Scalar

            >>> vfvXs = vfv * 3
            >>> print vfvXs
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print isinstance(vfvXs, VectorFaceVariable)
            1
            >>> sXvfv = 3 * vfv
            >>> print sXvfv
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print isinstance(sXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * Vector

            >>> vfvXv2 = vfv * (3,2)
            >>> print vfvXv2
            [[  0.,  2.,]
             [  3.,  4.,]
             [  6.,  6.,]
             [  9.,  8.,]
             [  3.,  6.,]
             [  6.,  8.,]
             [  9., 10.,]
             [ 18., 18.,]
             [  6., 12.,]
             [  3.,  6.,]]
            >>> print isinstance(vfvXv2, VectorFaceVariable)
            1
            >>> v2Xvfv = (3,2) * vfv
            >>> print v2Xvfv
            [[  0.,  2.,]
             [  3.,  4.,]
             [  6.,  6.,]
             [  9.,  8.,]
             [  3.,  6.,]
             [  6.,  8.,]
             [  9., 10.,]
             [ 18., 18.,]
             [  6., 12.,]
             [  3.,  6.,]]
            >>> print isinstance(v2Xvfv, VectorFaceVariable)
            1
            
            >>> vfvXv3 = vfv * (2,1,0)
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3Xvfv = (2,1,0) * vfv
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int


            >>> vfvXv10 = vfv * (9,8,7,6,5,4,3,2,1,0)
            >>> print vfvXv10
            [[  0.,  9.,]
             [  8., 16.,]
             [ 14., 21.,]
             [ 18., 24.,]
             [  5., 15.,]
             [  8., 16.,]
             [  9., 15.,]
             [ 12., 18.,]
             [  2.,  6.,]
             [  0.,  0.,]]
            >>> isinstance(vfvXv10, VectorFaceVariable)
            1
            >>> v10Xvfv = (9,8,7,6,5,4,3,2,1,0) * vfv
            >>> print v10Xvfv
            [[  0.,  9.,]
             [  8., 16.,]
             [ 14., 21.,]
             [ 18., 24.,]
             [  5., 15.,]
             [  8., 16.,]
             [  9., 15.,]
             [ 12., 18.,]
             [  2.,  6.,]
             [  0.,  0.,]]
            >>> isinstance(v10Xvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * `Variable` Scalar

            >>> vfvXsv = vfv * Variable(value = 3)
            >>> print vfvXsv
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print isinstance(vfvXsv, VectorFaceVariable)
            1
            >>> svXvfv = Variable(value = 3) * vfv
            >>> print svXvfv
            [[  0.,  3.,]
             [  3.,  6.,]
             [  6.,  9.,]
             [  9., 12.,]
             [  3.,  9.,]
             [  6., 12.,]
             [  9., 15.,]
             [ 18., 27.,]
             [  6., 18.,]
             [  3.,  9.,]]
            >>> print isinstance(svXvfv, VectorFaceVariable)
            1

        `VectorFaceVariable` * `Variable` Vector
            
            >>> vfvXv2v = vfv * Variable(value = (3,2))
            >>> print vfvXv2v
            [[  0.,  2.,]
             [  3.,  4.,]
             [  6.,  6.,]
             [  9.,  8.,]
             [  3.,  6.,]
             [  6.,  8.,]
             [  9., 10.,]
             [ 18., 18.,]
             [  6., 12.,]
             [  3.,  6.,]]
            >>> print isinstance(vfvXv2v, VectorFaceVariable)
            1
            >>> v2vXvfv = Variable(value = (3,2)) * vfv
            >>> print v2vXvfv
            [[  0.,  2.,]
             [  3.,  4.,]
             [  6.,  6.,]
             [  9.,  8.,]
             [  3.,  6.,]
             [  6.,  8.,]
             [  9., 10.,]
             [ 18., 18.,]
             [  6., 12.,]
             [  3.,  6.,]]
            >>> print isinstance(v2vXvfv, VectorFaceVariable)
            1
            
            >>> vfvXv3v = vfv * Variable(value = (2,1,0))
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            >>> v3vXvfv = Variable(value = (2,1,0)) * vfv
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'


            >>> vfvXv10v = vfv * Variable(value = (9,8,7,6,5,4,3,2,1,0))
            >>> print vfvXv10v
            [[  0.,  9.,]
             [  8., 16.,]
             [ 14., 21.,]
             [ 18., 24.,]
             [  5., 15.,]
             [  8., 16.,]
             [  9., 15.,]
             [ 12., 18.,]
             [  2.,  6.,]
             [  0.,  0.,]]
            >>> isinstance(vfvXv10v, VectorFaceVariable)
            1
            >>> v10vXvfv = Variable(value = (9,8,7,6,5,4,3,2,1,0)) * vfv
            >>> print v10vXvfv
            [[  0.,  9.,]
             [  8., 16.,]
             [ 14., 21.,]
             [ 18., 24.,]
             [  5., 15.,]
             [  8., 16.,]
             [  9., 15.,]
             [ 12., 18.,]
             [  2.,  6.,]
             [  0.,  0.,]]
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
            [ 9., 6.,]
            >>> print isinstance(sXv2v, Variable)
            1
            >>> v2vXs = Variable(value = (3,2)) * 3
            >>> print v2vXs
            [ 9., 6.,]
            >>> print isinstance(v2vXs, Variable)
            1
            
            
            
        Vector * `Variable` Scalar

            >>> vXsv = (3, 2) * Variable(value = 3)
            >>> print vXsv
            [ 9., 6.,]
            >>> print isinstance(vXsv, Variable)
            1
            >>> svXv = Variable(value = 3) * (3, 2)
            >>> print svXv
            [ 9., 6.,]
            >>> print isinstance(svXv, Variable)
            1

        Vector * `Variable` Vector
            
            >>> vXv2v = (3, 2) * Variable(value = (3,2))
            >>> print vXv2v
            [ 9., 4.,]
            >>> print isinstance(vXv2v, Variable)
            1
            >>> v2vXv = Variable(value = (3,2)) * (3, 2)
            >>> print v2vXv
            [ 9., 4.,]
            >>> print isinstance(v2vXv, Variable)
            1

            >>> vXv3v = (3, 2, 1) * Variable(value = (3,2))
            Traceback (most recent call last):
                  ...
            TypeError: can't multiply sequence to non-int
            >>> v3vXv = Variable(value = (3,2)) * (3, 2, 1) 
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
            [ 9., 6.,]
            >>> print isinstance(svXv2v, Variable)
            1
            >>> v2vXsv = Variable(value = (3,2)) * Variable(value = 3)
            >>> print v2vXsv
            [ 9., 6.,]
            >>> print isinstance(v2vXsv, Variable)
            1

            
        `Variable` Vector * `Variable` Vector
            
            >>> v2vXv2v = Variable(value = (3, 2)) * Variable(value = (3,2))
            >>> print v2vXv2v
            [ 9., 4.,]
            >>> print isinstance(v2vXv2v, Variable)
            1
            
            >>> v3vXv2v = Variable(value = (3, 2, 1)) * Variable(value = (3,2))
            Traceback (most recent call last):
                  ...
            TypeError: unsupported operand type(s) for *: 'instance' and 'instance'
            
        :Parameters:
          - `op`: the operator function to apply (takes two arguments for `self` and `other`)
          - `other`: the quantity to be operated with
          - `baseClass`: the `Variable` class that the binary operator should inherit from 
          - `opShape`: the shape that should result from the operation
          - `valueMattersForShape`: tuple of elements that must have a particular value for the operation to succeed.
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
        
        selfArray = self._getArrayAsOnes(self, valueMattersForShape)
        otherArray = self._getArrayAsOnes(other, valueMattersForShape)

        try:
            var0, var1 = self._verifyShape(op, var0, var1, selfArray, otherArray, opShape, otherClass)
        except SyntaxError:
            return NotImplemented
            
        
        # obtain a general operator class with the desired base class
	operatorClass = self._getOperatorVariableClass(baseClass)
        
        # declare a binary operator class with the desired base class
	class binOp(operatorClass):
	    def _calcValue(self):
		if isinstance(self.var[1], Variable):
		    val1 = self.var[1].getValue()
		else:
		    if type(self.var[1]) is type(''):
			self.var[1] = fipy.tools.dimensions.physicalField.PhysicalField(value = self.var[1])
		    val1 = self.var[1]
		    
		return self.op(self.var[0].getValue(), val1)
	
            def _getRepresentation(self, style = "__repr__"):
                return "(" + operatorClass._getRepresentation(self, style = style) + ")"
		
        # return the binary operator variable instance
	return binOp(op, [var0, var1])
	
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
	return self._getBinaryOperatorVariable(lambda a,b: a**b, other)
	    
    def __rpow__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: b**a, other)
	    
    def __div__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: a/b, other)
	
    def __rdiv__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: b/a, other)
	    
    def __neg__(self):
	return self._getUnaryOperatorVariable(lambda a: -a)
	
    def __pos__(self):
	return self
	
    def __abs__(self):
	return self._getUnaryOperatorVariable(lambda a: abs(a))

    def __lt__(self,other):
	"""
	Test if a `Variable` is less than another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a < 4)
	    >>> b
	    (Variable(value = 3) < 4)
	    >>> b()
	    1
	    >>> a.setValue(4)
	    >>> b()
	    0
	    
	Python automatically reverses the arguments when necessary
	
	    >>> 4 > Variable(value = 3)
	    (Variable(value = 3) < 4)
	"""
	return self._getBinaryOperatorVariable(lambda a,b: a<b, other)

    def __le__(self,other):
	"""
	Test if a `Variable` is less than or equal to another quantity
	
	    >>> a = Variable(value = 3)
	    >>> b = (a <= 4)
	    >>> b
	    (Variable(value = 3) <= 4)
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
	    (Variable(value = 3) == 4)
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
	    (Variable(value = 3) != 4)
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
	    (Variable(value = 3) > 4)
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
	    (Variable(value = 3) >= 4)
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
        return self._getBinaryOperatorVariable(lambda a,b: a & b, other)
        
    def __len__(self):
	return len(self.value)
	
    def __float__(self):
	return float(self.value)

    def arccos(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arccos(a))

    def arccosh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arccosh(a))

    def arcsin(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arcsin(a))

    def arcsinh(self):
        return self._getUnaryOperatorVariable(lambda a: numerix.arcsinh(a))

    def sqrt(self):
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
        return self._getUnaryOperatorVariable(lambda a: numerix.conjugate(a))

    def arctan2(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.arctan2(a,b), other)
		
    def dot(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: numerix.dot(a,b), other)
        
    def reshape(self, shape):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.reshape(a,b), shape, valueMattersForShape = (shape,))
        
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
           [[ 2., 1.,]
            [ 3., 2.,]
            [ 4., 3.,]
            [ 5., 4.,]]
           >>> var[3, 0] = 1
           >>> print v10
           [[ 2., 1.,]
            [ 3., 2.,]
            [ 4., 3.,]
            [ 5., 1.,]]
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
            return self._getUnaryOperatorVariable(lambda a: numerix.take(a, ids, axis = axis))
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
           1
           >>> print var.allclose((1,1,0))
           0
           
        """
        


        return self._getBinaryOperatorVariable(lambda a,b: numerix.allclose(a, b, atol = atol, rtol = rtol), 
                                               other, 
                                               baseClass = Variable,
                                               opShape = "number")
        
    def allequal(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.allequal(a,b), 
                                               other,
                                               baseClass = Variable,
                                               opShape = "number")

    def getMag(self):
        if self.mag is None:
	    self.mag = self.dot(self).sqrt()
	    
	return self.mag

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
