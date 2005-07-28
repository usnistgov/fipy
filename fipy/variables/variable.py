#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "variable.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 7/13/05 {4:52:33 PM} 
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
    
    Using a Variable in a mathematical expression will create an automatic
    dependency Variable, e.g.,
    
	>>> a = Variable(value = 3)
	>>> b = 4 * a
	>>> b
	(Variable(value = 3) * 4)
	>>> b()
	12
	
    Changes to the value of a Variable will automatically trigger changes in any dependent Variables
    
	>>> a.setValue(5)
	>>> b
	(Variable(value = 5) * 4)
	>>> b()
	20
	
    """
    
    def __init__(self, value=0., unit = None, array = None, name = '', mesh = None):
	"""
	Create a Variable.
	
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
	  - `unit`: the physical units of the variable
	  - `array`: the storage array for the variable
	  - `name`: the user-readable name of the variable
	  - `mesh`: the mesh that defines the geometry of this variable
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
        
	self.transposeVar = None
	self.sumVar = {}
	self.faceDifferences = {}
	self.laplacian = {}
        self.mag = None
    
    def getMesh(self):
	return self.mesh
	
    def __array__(self, t = None):
	"""
        Attempt to convert the Variable to a Numeric `array` object
        
            >>> v = Variable(value = [2,3])
            >>> Numeric.array(v)
            [ 2., 3.,]
            
        It is an error to convert a dimensional Variable to a 
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
	Make an duplicate of the Variable
	
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
	Return the value of the Variable with all units reduced to 
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

    def __getitem__(self, index):
        """    
        "Evaluate" the variable and return the specified element
        
            >>> a = Variable(value = ((3.,4.),(5.,6.)), unit = "m") + "4 m"
            >>> print a[1,1]
            10.0 m

        It is an error to slice a Variable whose `value` is not sliceable

            >>> Variable(value = 3)[2]
            Traceback (most recent call last):
                  ...
            TypeError: unsubscriptable object

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
	
    def __setitem__(self, index, value):
	self.value[index] = value
	self._markFresh()
	
    def __call__(self):
	"""
	"Evaluate" the Variable and return its value
	
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
	"Evaluate" the Variable and return its value (longhand)
	
	    >>> a = Variable(value = 3)
	    >>> a.getValue()
	    3
	    >>> b = a + 4
	    >>> b
	    (Variable(value = 3) + 4)
	    >>> b.getValue()
	    7
	"""
	self._refresh()
	return self.value

    def _setValue(self, value, unit = None, array = None):
	PF = fipy.tools.dimensions.physicalField.PhysicalField
	if not isinstance(value, PF) \
        and (unit is not None or type(value) in [type(''), type(()), type([])]):
	    self.value = PF(value = value, unit = unit, array = array)
	elif array is not None:
	    array[:] = value
	    self.value = array
	else:
	    self.value = value
	    
	if isinstance(self.value, PF) and self.value.getUnit().isDimensionless():
	    self.value = self.value.getNumericValue()

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
	
    def _refresh(self):
	if self.stale:           
	    for required in self.requiredVariables:
		required._refresh()
	    self._calcValue()
	    self._markFresh()
		    
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

    def _getOperatorVariableClass(self, parentClass = None):
	if parentClass is None:
	    parentClass = self._getVariableClass()

	class OperatorVariable(parentClass):
	    def __init__(self, op, var, mesh = None):
		if mesh is None:
		    mesh = var[0].getMesh()
		self.op = op
		self.var = var
		parentClass.__init__(self, value = var[0], mesh = mesh)
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
                name = parentClass.getName(self)
                if len(name) == 0:
                    name = self._getRepresentation(style = "name")
                return name
                
	    def copy(self):
	       return self.__class__(
		   op = self.op,
		   var = self.var,
		   mesh = self.getMesh())

	return OperatorVariable
	
    def _getUnaryOperatorVariable(self, op, parentClass = None):
	class unOp(self._getOperatorVariableClass(parentClass)):
	    def _calcValue(self):
		self._setValue(value = self.op(self.var[0].getValue())) 
		
	return unOp(op, [self])
	    
    def _getBinaryOperatorVariable(self, op, other, parentClass = None):
	operatorClass = self._getOperatorVariableClass(parentClass)
	
	class binOp(operatorClass):
	    def _calcValue(self):
		if isinstance(self.var[1], Variable):
		    val1 = self.var[1].getValue()
		else:
		    if type(self.var[1]) is type(''):
			self.var[1] = fipy.tools.dimensions.physicalField.PhysicalField(value = self.var[1])
		    val1 = self.var[1]
		    
		self._setValue(value = self.op(self.var[0].getValue(), val1)) 
	
            def _getRepresentation(self, style = "__repr__"):
                return "(" + operatorClass._getRepresentation(self, style = style) + ")"
		
	return binOp(op, [self, other])
	
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
	Test if a Variable is less than another quantity
	
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
	Test if a Variable is less than or equal to another quantity
	
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
	Test if a Variable is equal to another quantity
	
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
	Test if a Variable is not equal to another quantity
	
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
	Test if a Variable is greater than another quantity
	
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
	Test if a Variable is greater than or equal to another quantity
	
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
        return self._getBinaryOperatorVariable(lambda a,b: numerix.reshape(a,b), shape)
        
    def transpose(self):
	if self.transposeVar is None:
	    from transposeVariable import _TransposeVariable
	    self.transposeVar = _TransposeVariable(self)
	
	return self.transposeVar

    def sum(self, index = 0):
	if not self.sumVar.has_key(index):
	    from sumVariable import _SumVariable
	    self.sumVar[index] = _SumVariable(self, index)
	
	return self.sumVar[index]
	
    def take(self, ids):
	return numerix.take(self.getValue(), ids)
	
##     def allclose(self, other, atol = 1.e-5, rtol = 1.e-8):
## 	return numerix.allclose(first = self.getValue(), second = other, atol = atol, rtol = rtol)
        
    def allclose(self, other, rtol = 1.e-10, atol = 1.e-10):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.allclose(a, b, atol = atol, rtol = rtol), other)
        
    def allequal(self, other):
        return self._getBinaryOperatorVariable(lambda a,b: numerix.allequal(a,b), other)

    def getMag(self):
        if self.mag is None:
	    self.mag = self.dot(self).sqrt()
	    
## 	    from magVariable import _MagVariable
## 	    self.mag = _MagVariable(self)
	    
	return self.mag

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
