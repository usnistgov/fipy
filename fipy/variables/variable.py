#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "variable.py"
 #                                    created: 11/10/03 {3:15:38 PM} 
 #                                last update: 3/9/04 {11:12:13 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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


import Numeric

import fivol.tools.dimensions.physicalField

import fivol.tools.array as array

class Variable:
    
    """Lazily evaluated physical field or quantity with units

    Constructor:

    - Variable(|mesh|, |value|, |unit|), where |value| is a number of
      arbitrary type and |unit| is a string containing the unit name.

    - PhysicalField(|mesh|, |string|), where |string| contains both the value
      and the unit. This form is provided to make interactive use more
      convenient.

    Variable instances allow addition, subtraction,
    multiplication, and division with each other as well as
    multiplication, division, and exponentiation with numbers.
    Addition and subtraction check that the units of the two operands
    are compatible and return the result in the units of the first
    operand. A limited set of mathematical functions (from module
    Numeric) is applicable as well:

    sqrt -- equivalent to exponentiation with 0.5.

    sin, cos, tan -- applicable only to objects whose unit is compatible
		     with 'rad'.
    """
    
    def __init__(self, value=0., unit = None, array = None, name = '', mesh = None):
	self.requiredVariables = []
	self.subscribedVariables = []

	self.value = self.getPhysicalFieldClass()(value = value, unit = unit, array = array)
## 	self.value = value
## 	if array is not None:
## 	    array[:] = self.value
## 	    self.value = array
		
	self.mesh = mesh
	self.name = name
		
	self.stale = 1
	self.markFresh()
        
	self.transposeVar = None
	self.sumVar = {}
        self.mag = None
	
    def getPhysicalFieldClass(self):
	return fivol.tools.dimensions.physicalField.PhysicalField
	
    def copy(self):
	return Variable(self)
	
    def getMesh(self):
	return self.mesh
	
    def getUnit(self):
	value = self.getValue()
	if isinstance(value,fivol.tools.dimensions.physicalField.PhysicalField):
	    return value.getUnit()
	else:
	    return "1"
	
    def inBaseUnits(self):
	value = self.getValue()
	if isinstance(value,fivol.tools.dimensions.physicalField.PhysicalField):
	    return value.inBaseUnits()
	else:
	    return value

    def __getitem__(self, index): 
	return self.getValue()[index]
	
    def __str__(self):
	return str(self.getValue())
	    
    def __repr__(self):
	return (self.__class__.__name__ + '("' + self.name + '",' + `self.getValue()` + ')')
	
    def __setitem__(self, index, value):
	self.value[index] = value
	self.markFresh()
		
    def getValue(self):
	self.refresh()
	return self.value
	
    def setValue(self,value):
	self.value = value
	self.markFresh()
	    
    def getNumericValue(self):
	value = self.getValue()
	if isinstance(value,fivol.tools.dimensions.physicalField.PhysicalField):
	    return value.getNumericValue()
	else:
	    return value
	
    def refresh(self):
	if self.stale:           
	    for required in self.requiredVariables:
		required.refresh()
	    self.calcValue()
	    self.markFresh()
		    
    def calcValue(self):
	pass
	
    def markFresh(self):
	self.stale = 0
	for subscriber in self.subscribedVariables:
	    subscriber.markStale() 

    def markStale(self):
	if not self.stale:
	    self.stale = 1
	    for subscriber in self.subscribedVariables:
		subscriber.markStale()
	    
    def requires(self, var):
	if isinstance(var, Variable):
	    self.requiredVariables.append(var)
	    var.requiredBy(self)
	    self.markStale()
	return var
	    
    def requiredBy(self, var):
	assert isinstance(var, Variable)
	self.subscribedVariables.append(var)
	
    def getVariableClass(self):
	return Variable
	
    def getUnaryOperatorVariable(self, op):
        parentClass = self.getVariableClass()

	class unOp(parentClass):
	    def __init__(self, op, var, mesh = None):
		if mesh is None:
		    mesh = var.getMesh()
		self.op = op
		self.var = var
		# this horrendous hack is necessary because older Pythons
		# (2.1) don't know the value of 'parentClass' at this
		# point.  Since we're good and dog-fearing people, we don't
		# ever, ever, ever do multiple inheritance, so we know
		# there is only one base class.
		self.__class__.__bases__[0].__init__(self, mesh = mesh)
		self.requires(self.var)
		
	    def calcValue(self):
		self.value = self.getPhysicalFieldClass()(self.op(self.var.getValue()))
		
	    def __repr__(self):
		return ("\n" + `self.op` + "(" + `self.var` + ") = " + `self.value`)
		
	return unOp(op, self)
	    
    def getBinaryOperatorVariable(self, op, var2):
	parentClass = self.getVariableClass()
	
	class binOp(parentClass):
	    def __init__(self, op, var1, var2, mesh = None):
		if mesh is None:
		    mesh = var1.getMesh()
		self.op = op
		self.var1 = var1
		self.var2 = var2
		# this horrendous hack is necessary because older Pythons
		# (2.1) don't know the value of 'parentClass' at this
		# point.  Since we're good and dog-fearing people, we don't
		# ever, ever, ever do multiple inheritance, so we know
		# there is only one base class.
		self.__class__.__bases__[0].__init__(self, mesh = mesh)
		self.requires(self.var1)
		self.requires(self.var2)

	    def calcValue(self):
		if isinstance(self.var2, Variable):
		    val2 = self.var2.getValue()
		else:
		    val2 = self.var2
		    
		self.value = self.getPhysicalFieldClass()(self.op(self.var1.getValue(), val2))
		
	    def __repr__(self):
		return ("\n" + `self.op` + "(\n\t" + `self.var1` + ",\n\t" + `self.var2` + "\n) = " + `self.getValue()`)
		
	return binOp(op, self, var2)
	
    def __add__(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: a+b, other)
	
    __radd__ = __add__

    def __sub__(self, other):
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
	return self.getBinaryOperatorVariable(lambda a,b: a<b, other)

    def __le__(self,other):
	return self.getBinaryOperatorVariable(lambda a,b: a<=b, other)
	
    def __eq__(self,other):
	return self.getBinaryOperatorVariable(lambda a,b: a==b, other)
	
    def __ne__(self,other):
	return self.getBinaryOperatorVariable(lambda a,b: a!=b, other)
	
    def __gt__(self,other):
	return self.getBinaryOperatorVariable(lambda a,b: a>b, other)
	
    def __ge__(self,other):
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

    def sin(self):
	return self.getUnaryOperatorVariable(lambda a: array.sin(a))
		
    def cos(self):
	return self.getUnaryOperatorVariable(lambda a: array.cos(a))
		
    def dot(self, other):
	return self.getBinaryOperatorVariable(lambda a,b: array.dot(a,b))
	    
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

    def getMag(self):
        if self.mag is None:
	    from magVariable import MagVariable
	    self.mag = MagVariable(self)
	    
	return self.mag

        

