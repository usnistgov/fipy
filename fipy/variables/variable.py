"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "variable.py"
                                   created: 11/10/03 {3:15:38 PM} 
                               last update: 12/22/03 {3:12:10 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-10 JEG 1.0 original
###################################################################
"""

from tools.dimensionalization import PhysicalField
from Scientific.Physics.PhysicalQuantities import isPhysicalQuantity
# from binaryOperatorVariable import BinaryOperatorVariable
import Numeric

class Variable:
# class Variable(PhysicalField):
    
    def __init__(self, mesh, name = '', value=0., array = None, scaling = None, unit = None):
	self.mesh = mesh
	self.name = name
	
	self.requiredVariables = []
	self.subscribedVariables = []

# 	if type(value) not in [type(1),type(1.),type(Numeric.array((1)))] and not isPhysicalQuantity(value):
# 	    value = PhysicalField(value)
	    
	self.scaling = scaling
# 	if scaling is not None:
# 	    self.scaling = PhysicalField(scaling)
# 	    if unit is not None:
# 		self.scaling = scaling.inUnitsOf(unit)
# 	    else:
# 		unit = self.scaling.unit
# 	elif unit is not None:
# 	    self.scaling = PhysicalField(1., unit)
# 	else:
# 	    self.scaling = 1
# 	if isPhysicalQuantity(value):
# 	    unit = value.unit
# 	    value = value.value
# 	else:
# 	    unit = "m/m"
	    
	if array is None:
	    array = Numeric.array(value)
	else:
	    array[:] = value

	self.value = array
# 	PhysicalField.__init__(self, array, unit)
		
	self.markFresh()
	
	self.transposeVar = None
	self.sumVar = {}
	
    def __getitem__(self, index): 
	return self.getValue()[index]
	
    def __repr__(self):
	return (self.__class__.__name__ + '(' + `self.getValue()` + ')')
	
    def __setitem__(self, index, value): 
	self.value[index] = value
# 	PhysicalField.__setitem__(self, index, value)
	self.markFresh()
		
    def getMesh(self):
	return self.mesh

    def getValue(self):
	self.refresh()
        return self.value
	
    def __float__(self):
	self.refresh()
	return self.value
	
    def refresh(self):
	if self.stale:
	    for required in self.requiredVariables:
		required.refresh()
	    self.calcValue()
	    self.markFresh()
		    
    def calcValue(self):
	pass
	
    def markFresh(self):
# 	print self, "is fresh"
	self.stale = 0
	for subscriber in self.subscribedVariables:
	    subscriber.markStale() 

    def markStale(self):
# 	print self, "is stale"
	self.stale = 1
	for subscriber in self.subscribedVariables:
	    subscriber.markStale()
	    
    def requires(self, var):
	if isinstance(var, Variable):
	    self.requiredVariables.append(var)
	    var.requiredBy(self)
# 	    print self, "requires", self.requiredVariables
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
		# this horrendous hack is necessary because older Python's
		# (2.1) don't know the value of 'parentClass' at this
		# point.  Since we're good and dog-fearing people, we don't
		# ever, ever, ever do multiple inheritance, so we know
		# there is only one base class.
		self.__class__.__bases__[0].__init__(self, mesh = mesh)
		self.var = self.requires(var)
		
	    def calcValue(self):
		self.value = self.op(self.var.getValue())
		
	    def __repr__(self):
		return ("\n" + `self.op`)
# 		return ("\n" + `self.op` + "(" + `self.var` + ") = " + `self.value`)
		
	return unOp(op, self)
	    
    def getBinaryOperatorVariable(self, op, var2):
	parentClass = self.getVariableClass()

	class binOp(parentClass):
	    def __init__(self, op, var1, var2, mesh = None, parentClass = None):
		if mesh is None:
		    mesh = var1.getMesh()
		self.op = op
		# this horrendous hack is necessary because older Python's
		# (2.1) don't know the value of 'parentClass' at this
		# point.  Since we're good and dog-fearing people, we don't
		# ever, ever, ever do multiple inheritance, so we know
		# there is only one base class.
		self.__class__.__bases__[0].__init__(self, mesh = mesh)
		self.var1 = self.requires(var1)
		self.var2 = self.requires(var2)
		
	    def calcValue(self):
		if isinstance(self.var2, Variable):
		    val2 = self.var2.getValue()
		else:
		    val2 = self.var2
		    
		self.value = self.op(self.var1.getValue(), val2)
		
	    def __repr__(self):
		return ("\n" + `self.op`)
# 		return ("\n" + `self.op` + "(" + `self.var1` + "," + `self.var2` + ") = " + `self.value`)
		
	return binOp(op, self, var2, parentClass = parentClass)
	
    def __add__(self, other):
	return self.getBinaryOperatorVariable(Numeric.add, other)
	
    def __radd__(self, other):
	return self.__add__(other)

    def __sub__(self, other):
	return self.getBinaryOperatorVariable(Numeric.subtract, other)
	
    def __rsub__(self, other):
	return -self + other
	    
    def __mul__(self, other):
	return self.getBinaryOperatorVariable(Numeric.multiply, other)
	
    def __rmul__(self, other):
	return self.__mul__(other)
	    
    def __mod__(self, other):
	return self.getBinaryOperatorVariable(Numeric.fmod, other)
	    
#     def __rmod__(self, other):
# 	return self.__mod__(other)
	    
    def __pow__(self, other):
	return self.getBinaryOperatorVariable(Numeric.power, other)
	    
    def __rpow__(self, other):
	return self.__pow__(other)
	    
    def __div__(self, other):
	return self.getBinaryOperatorVariable(Numeric.divide, other)
	
    def __rdiv__(self, other):
	return self**-1 * other
	    
    def __neg__(self):
	return -1 * self
	
    def __pos__(self):
	return self
	
    def __abs__(self):
	return self.getUnaryOperatorVariable(Numeric.fabs)
	
    def __lt__(self,other):
	return self.getBinaryOperatorVariable(Numeric.less, other)

    def __le__(self,other):
        return self.getBinaryOperatorVariable(Numeric.less_equal, other)
	
    def __eq__(self,other):
        return self.getBinaryOperatorVariable(Numeric.equal, other)
	
    def __ne__(self,other):
        return self.getBinaryOperatorVariable(Numeric.not_equal, other)
	
    def __gt__(self,other):
        return self.getBinaryOperatorVariable(Numeric.greater, other)
	
    def __ge__(self,other):
        return self.getBinaryOperatorVariable(Numeric.greater_equal, other)

    def tan(self):
	return self.getUnaryOperatorVariable(Numeric.tan)
	
    def transpose(self):
	if self.transposeVar is None:
	    from transposeVariable import TransposeVariable
	    self.transposeVar = TransposeVariable(self)
	
	return self.transposeVar

    def sum(self, index):
	if not self.sumVar.has_key(index):
	    from sumVariable import SumVariable
	    self.sumVar[index] = SumVariable(self, index)
	
	return self.sumVar[index]
