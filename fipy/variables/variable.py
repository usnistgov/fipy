#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "variable.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 2/25/05 {5:47:25 PM} 
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

import fipy.tools.array as array

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
    
    __variable__ = True
    
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
	self.markFresh()
        
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
	    [2,3,]
	    
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
	self.markFresh()
	
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
	self.refresh()
	return self.value

    def _setValue(self, value, unit = None, array = None):
	PF = fipy.tools.dimensions.physicalField.PhysicalField
	if not isinstance(value, PF) and (unit is not None or type(value) is type('')):
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
	self.markFresh()
	
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
	
    def refresh(self):
	if self.stale:           
	    for required in self.requiredVariables:
		required.refresh()
	    self._calcValue()
	    self.markFresh()
		    
    def _calcValue(self):
	pass
	
    def _markStale(self):
        import weakref
        remainingSubscribedVariables = []
        for subscriber in self.subscribedVariables:
            try:
                subscriber.markStale() 
                remainingSubscribedVariables.append(subscriber)
            except weakref.ReferenceError:
                pass
        self.subscribedVariables = remainingSubscribedVariables

    def markFresh(self):
	self.stale = 0
        self._markStale()

    def markStale(self):
	if not self.stale:
	    self.stale = 1
            self._markStale()
	    
    def requires(self, var):
	if isinstance(var, Variable):
	    self.requiredVariables.append(var)
	    var.requiredBy(self)
	    self.markStale()
	return var
	    
    def requiredBy(self, var):
	assert isinstance(var, Variable)
        
        # we retain a weak reference to avoid a memory leak 
        # due to circular references between the subscriber
        # and the subscribee
        import weakref
	self.subscribedVariables.append(weakref.proxy(var))
	
    def getVariableClass(self):
	return Variable
	
    def getOperatorVariableClass(self, parentClass = None):
	if parentClass is None:
	    parentClass = self.getVariableClass()

	class OperatorVariable(parentClass):
	    def __init__(self, op, var, mesh = None):
		if mesh is None:
		    mesh = var[0].getMesh()
		self.op = op
		self.var = var
		parentClass.__init__(self, value = var[0], mesh = mesh)
		for aVar in self.var:
		    self.requires(aVar)
		    
	    def __repr__(self):
		bytecodes = [ord(byte) for byte in self.op.func_code.co_code]
		
		def _getIndex():
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
		    
		    if bytecode == 13:
			# UNARY_CONVERT
			stack.append("`" + stack.pop() + "`")
		    elif bytecode == 25:	
			# BINARY_SUBSCR
			stack.append(stack.pop(-2) + "[" + stack.pop() + "]")
		    elif bytecode == 83:	
			# RETURN_VALUE
			return stack.pop()
		    elif bytecode == 105:
			# LOAD_ATTR
			stack.append(stack.pop() + "." + self.op.func_code.co_names[_getIndex()])
		    elif bytecode == 106:
			# COMPARE_OP
			import dis
			stack.append(stack.pop(-2) + " " + dis.cmp_op[_getIndex()] + " " + stack.pop())
		    elif bytecode == 116:
			# LOAD_GLOBAL
			stack.append(self.op.func_code.co_names[_getIndex()])
		    elif bytecode == 124:	
			# LOAD_FAST
			stack.append(repr(self.var[_getIndex()]))
		    elif bytecode == 131:
			# CALL_FUNCTION
			args = []
			for j in range(bytecodes.pop(1)):
			    # keyword parameters
			    args.insert(0, stack.pop(-2) + " = " + stack.pop())
			for j in range(bytecodes.pop(0)):
			    # positional parameters
			    args.insert(0, stack.pop())
			stack.append(stack.pop() + "(" + ", ".join(args) + ")")
		    elif unop.has_key(bytecode):
			stack.append(unop[bytecode] + stack.pop())
		    elif binop.has_key(bytecode):
			stack.append(stack.pop(-2) + " " + binop[bytecode] + " " + stack.pop())
		    else:
			raise SyntaxError, "Unknown bytecode: " + `bytecode` + " in " + `[ord(byte) for byte in self.op.func_code.co_code]`

	    def copy(self):
	       return self.__class__(
		   op = self.op,
		   var = self.var,
		   mesh = self.getMesh())

	return OperatorVariable
	
    def getUnaryOperatorVariable(self, op, parentClass = None):
	class unOp(self.getOperatorVariableClass(parentClass)):
	    def _calcValue(self):
		self._setValue(value = self.op(self.var[0].getValue())) 
		
	return unOp(op, [self])
	    
    def getBinaryOperatorVariable(self, op, other, parentClass = None):
	operatorClass = self.getOperatorVariableClass(parentClass)
	
	class binOp(operatorClass):
	    def _calcValue(self):
		if isinstance(self.var[1], Variable):
		    val1 = self.var[1].getValue()
		else:
		    if type(self.var[1]) is type(''):
			self.var[1] = fipy.tools.dimensions.physicalField.PhysicalField(value = self.var[1])
		    val1 = self.var[1]
		    
		self._setValue(value = self.op(self.var[0].getValue(), val1)) 
	
	    def __repr__(self):
		return "(" + operatorClass.__repr__(self) + ")"
		
	return binOp(op, [self, other])
	
    def __add__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return other + self
        else:
            return self.getBinaryOperatorVariable(lambda a,b: a+b, other)
	
    __radd__ = __add__

    def __sub__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return -other + self
        else:
            return self.getBinaryOperatorVariable(lambda a,b: a-b, other)
	
    def __rsub__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: b-a, other)
	    
    def __mul__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: a*b, other)
	
    __rmul__ = __mul__
	    
    def __mod__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: a%b, other)
	    
    def __pow__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: a**b, other)
	    
    def __rpow__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: b**a, other)
	    
    def __div__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: a/b, other)
	
    def __rdiv__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: b/a, other)
	    
    def __neg__(self):
	return self.getUnaryOperatorVariable(lambda a: -a)
	
    def __pos__(self):
	return self
	
    def __abs__(self):
	return self.getUnaryOperatorVariable(lambda a: abs(a))

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
	return self.getBinaryOperatorVariable(lambda a,b: a<b, other)

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
	return self.getBinaryOperatorVariable(lambda a,b: a<=b, other)
	
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
	return self.getBinaryOperatorVariable(lambda a,b: a==b, other)
	
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
	return self.getBinaryOperatorVariable(lambda a,b: a!=b, other)
	
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
	return self.getBinaryOperatorVariable(lambda a,b: a>b, other)
	
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
	return self.getBinaryOperatorVariable(lambda a,b: a>=b, other)

    def __len__(self):
	return len(self.value)
	
    def __float__(self):
	return float(self.value)
	
    def sqrt(self):
	return self.getUnaryOperatorVariable(lambda a: array.sqrt(a))
	
    def tan(self):
	return self.getUnaryOperatorVariable(lambda a: array.tan(a))

    def arctan(self):
	return self.getUnaryOperatorVariable(lambda a: array.arctan(a))

    def exp(self):
	return self.getUnaryOperatorVariable(lambda a: array.exp(a))

    def sin(self):
	return self.getUnaryOperatorVariable(lambda a: array.sin(a))
		
    def cos(self):
	return self.getUnaryOperatorVariable(lambda a: array.cos(a))

    def arctan2(self, other):
        return self.getBinaryOperatorVariable(lambda a,b: array.arctan2(a,b), other)
		
    def dot(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: array.dot(a,b), other)

    def min(self):
        return self.getUnaryOperatorVariable(lambda a: array._min(a), parentClass = Variable)
        
    def max(self):
        return self.getUnaryOperatorVariable(lambda a: array._max(a), parentClass = Variable)
        
    def transpose(self):
	if self.transposeVar is None:
	    from transposeVariable import TransposeVariable
	    self.transposeVar = TransposeVariable(self)
	
	return self.transposeVar

    def sum(self, index = 0):
	if not self.sumVar.has_key(index):
	    from sumVariable import SumVariable
	    self.sumVar[index] = SumVariable(self, index)
	
	return self.sumVar[index]
	
    def take(self, ids):
	return array.take(self.getValue(), ids)
	
    def allclose(self, other, atol = 1.e-5, rtol = 1.e-8):
	return array.allclose(first = self.getValue(), second = other, atol = atol, rtol = rtol)

    def getMag(self):
        if self.mag is None:
	    self.mag = self.dot(self).sqrt()
	    
## 	    from magVariable import MagVariable
## 	    self.mag = MagVariable(self)
	    
	return self.mag

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
