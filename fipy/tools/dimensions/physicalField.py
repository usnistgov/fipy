#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "physicalField.py"
 #                                    created: 12/28/03 {10:56:55 PM} 
 #                                last update: 1/28/04 {4:21:15 PM} 
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
 #  Copyright 1997-2004 by Konrad Hinsen, except as noted below.
 # 
 #  Permission to use, copy, modify, and distribute this software and its
 #  documentation for any purpose and without fee is hereby granted,
 #  provided that the above copyright notice appear in all copies and that
 #  both that copyright notice and this permission notice appear in
 #  supporting documentation.
 # 
 #  THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 #  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
 #  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 #  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
 #  USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 #  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 #  PERFORMANCE OF THIS SOFTWARE.
 #  
 #  Description: 
 # 
 # Physical fields or quantities with units
 #
 # Based on PhysicalQuantities of the Scientific package, written by Konrad
 # Hinsen <hinsen@cnrs-orleans.fr>
 #
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-12-28 JEG 1.0 original
 #  1998/09/29 GPW     now supports conversions with offset 
 #                     (for temperature units)
 #  1998/09/28 GPW     now removes __args__ from local dict after eval
 # ###################################################################
 ##


"""Physical quantities with units.

This module provides a data type that represents a physical
quantity together with its unit. It is possible to add and
subtract these quantities if the units are compatible, and
a quantity can be converted to another compatible unit.
Multiplication, subtraction, and raising to integer powers
is allowed without restriction, and the result will have
the correct unit. A quantity can be raised to a non-integer
power only if the result can be represented by integer powers
of the base units.

The values of physical constants are taken from the 1986
recommended values from CODATA. Other conversion factors
(e.g. for British units) come from various sources. I can't
guarantee for the correctness of all entries in the unit
table, so use this at your own risk!
"""

import re, string, umath

import Numeric

import fivol.variables.variable

from NumberDict import NumberDict

# Class definitions

class PhysicalField:

    """Physical field or quantity with units

    Constructor:

    - PhysicalField(|value|, |unit|), where |value| is a number of
      arbitrary type and |unit| is a string containing the unit name.

    - PhysicalField(|string|), where |string| contains both the value
      and the unit. This form is provided to make interactive use more
      convenient.

    PhysicalField instances allow addition, subtraction,
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

    def __init__(self, value, unit = None, array = None):
	if isinstance(value, PhysicalField):
	    unit = value.unit
	    value = value.value
	elif unit is not None:
	    unit = _findUnit(unit)
	elif type(value) is type(''):
	    s = string.strip(value)
	    match = PhysicalField._number.match(s)
	    if match is None:
		value = 1
		unit = _findUnit(s)
# 		    raise TypeError, 'No number found'
	    else:
		value = string.atof(match.group(0))
		unit = _findUnit(s[len(match.group(0)):])
	    
	if type(value) in [type([]),type(())]:
	    value = [PhysicalField(item,unit) for item in value]
	    if unit is None:
		unit = value[0].unit
	    value = [item.inUnitsOf(unit) / PhysicalField(1,unit) for item in value]
	    value = Numeric.array(value)
	    
	if unit is None:
	    unit = _findUnit("")

	self.value = value
	self.unit = unit
	if array is not None:
	    array[:] = self.value
	    self.value = array
# 	self.value = Numeric.array(self.value)

    _number = re.compile('[+-]?[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?')

    def copy(self):
	return PhysicalField(self)
	
    def __str__(self):
	return str(self.value) + ' ' + self.unit.name()

    def __repr__(self):
 	return (self.__class__.__name__ + '(' + `self.value` + ',' + 
 		`self.unit.name()` + ')')

    def _sum(self, other, sign1 = lambda a: a, sign2 = lambda b: b):
	if _isVariable(other):
	    return sign2(other) + self.__class__(value = sign1(self.value), unit = self.unit)
	if type(other) is type(''):
	    other = PhysicalField(value = other)
	if not isinstance(other,PhysicalField):
	    if Numeric.alltrue(other == 0):
		new_value = sign1(self.value)
	    elif self.unit.isDimensionlessOrAngle():
		# stupid Numeric bug
		# it's smart enough to negate a boolean, but not 
		# smart enough to know that it's negative when it
		# gets done.
		if type(self.value) is type(Numeric.array((0))) and self.value.typecode() is 'b':
		    self.value = 1. * other
		if type(other) is type(Numeric.array((0))) and other.typecode() is 'b':
		    other = 1. * other
		new_value = sign1(self.value) + sign2(other)
	    else:
		raise TypeError, str(self) + ' and ' + str(other) + ' are incompatible.'
	else:
	    new_value = sign1(self.value) + \
			sign2(other.value)*other.unit.conversionFactorTo(self.unit)
	return self.__class__(value = new_value, unit = self.unit)

    def __add__(self, other):
	return self._sum(other)

    __radd__ = __add__

    def __sub__(self, other):
	return self._sum(other, sign2 = lambda b: -b)

    def __rsub__(self, other):
	return self._sum(other, sign1 = lambda a: -a)

    def __mul__(self, other):
	if _isVariable(other):
	    return other.__mul__(self)
	if type(other) is type(''):
	    other = PhysicalField(value = other)
	if not isinstance(other,PhysicalField):
	    return self.__class__(value = self.value*other, unit = self.unit)
	value = self.value*other.value
	unit = self.unit*other.unit
	if unit.isDimensionless():
	    return value*unit.factor
	else:
	    return self.__class__(value = value, unit = unit)

    __rmul__ = __mul__

    def __div__(self, other):
	if _isVariable(other):
	    return other.__rdiv__(self)
	if type(other) is type(''):
	    other = self.__class__(value = other)
	if not isinstance(other,PhysicalField):
	    value = self.value/other 
	    unit = self.unit
	else:
	    value = self.value/other.value
	    unit = self.unit/other.unit
	if unit.isDimensionless():
	    return value*unit.factor
	else:
	    return self.__class__(value = value, unit = unit)

    def __rdiv__(self, other):
	if _isVariable(other):
	    return other.__div__(self)
	if type(other) is type(''):
	    other = PhysicalField(value = other)
	if not isinstance(other,PhysicalField):
	    value = other/self.value
	    unit = pow(self.unit, -1)
	else:
	    value = other.value/self.value
	    unit = other.unit/self.unit
	if unit.isDimensionless():
	    return value*unit.factor
	else:
	    return self.__class__(value = value, unit = unit)

    def __mod__(self, other):
	if _isVariable(other):
	    return other.__mod__(self)
	if type(other) is type(''):
	    other = self.__class__(value = other)
	if not isinstance(other,PhysicalField):
	    value = self.value % other 
	    unit = self.unit
	else:
	    value = self.value % other.value
	    unit = self.unit/other.unit
	if unit.isDimensionless():
	    return value*unit.factor
	else:
	    return self.__class__(value = value, unit = unit)
	    
    def __pow__(self, other):
	if type(other) is type(''):
	    other = PhysicalField(value = other)
## 	if isinstance(other,PhysicalField):
## 	    raise TypeError, 'Exponents must be dimensionless'
	return self.__class__(value = pow(self.value, float(other)), unit = pow(self.unit, other))

    def __rpow__(self, other):
	return pow(other, float(self))
## 	raise TypeError, 'Exponents must be dimensionless'

    def __abs__(self):
	return self.__class__(value = abs(self.value), unit = self.unit)

    def __pos__(self):
	return self

    def __neg__(self):
	return self.__class__(value = -self.value, unit = self.unit)

    def __nonzero__(self):
	return self.value != 0
	
##     def __float__(self):
## 	if self.unit.isDimensionless():
## 	    return float(self.value)
## 	elif self.value == 0:
## 	    return self.value
## 	else:
## 	    raise TypeError, 'Quantity has dimensions ' + str(self.unit) + '. Cannot be converted to type(float).'

    def _inMyUnits(self, other):
	if not isinstance(other,PhysicalField):
	    if type(other) is type(''):
		other = PhysicalField(other)
	    elif other == 0 or self.unit.isDimensionlessOrAngle():
		other = PhysicalField(value = other, unit = self.unit)
	    else:
		raise TypeError, 'Incompatible types'
	return other.inUnitsOf(self.unit)
	
    def __getitem__(self, index): 
	return PhysicalField(self.value[index],self.unit)
	
    def __setitem__(self, index, value):
	if type(value) is type(''):
	    value = PhysicalField(value)
	if isinstance(value,PhysicalField):
	    value = self._inMyUnits(value)
	    self.value[index] = value.value
	else:
	    self.value[index] = value
	    
    def __lt__(self,other):
	other = self._inMyUnits(other)
	return self.value < other.value

    def __le__(self,other):
	other = self._inMyUnits(other)
	return self.value <= other.value
	
    def __eq__(self,other):
	other = self._inMyUnits(other)
	return self.value == other.value
	
    def __ne__(self,other):
	other = self._inMyUnits(other)
	return self.value != other.value
	
    def __gt__(self,other):
	other = self._inMyUnits(other)
	return self.value > other.value
	
    def __ge__(self,other):
	other = self._inMyUnits(other)
	return self.value >= other.value
	    
    def __len__(self):
	return len(self.value)
	
    def convertToUnit(self, unit):
        """Changes the unit to |unit| and adjusts the value such that
        the combination is equivalent. The new unit is by a string containing
        its name. The new unit must be compatible with the previous unit
        of the object."""
	unit = _findUnit(unit)
	self.value = _convertValue (self.value, self.unit, unit)
	self.unit = unit

    def inUnitsOf(self, *units):
        """Returns one or more PhysicalField objects that express
        the same physical quantity in different units. The units are
        specified by strings containing their names. The units must be
        compatible with the unit of the object. If one unit is
        specified, the return value is a single PhysicalObject. If
        several units are specified, the return value is a tuple of
        PhysicalObject instances with with one element per unit such
        that the sum of all quantities in the tuple equals the the
        original quantity and all the values except for the last one
        are integers. This is used to convert to irregular unit
        systems like hour/minute/second. The original object will not
        be changed.
        """
	units = map(_findUnit, units)
	if len(units) == 1:
	    unit = units[0]
	    value = _convertValue (self.value, self.unit, unit)
	    return self.__class__(value = value, unit = unit)
	else:
	    units.sort()
	    result = []
	    value = self.value
	    unit = self.unit
	    for i in range(len(units)-1,-1,-1):
		value = value*unit.conversionFactorTo(units[i])
		if i == 0:
		    rounded = value
		else:
		    rounded = _round(value)
		result.append(self.__class__(value = rounded, unit = units[i]))
		value = value - rounded
		unit = units[i]
	    return tuple(result)

    def getUnit(self):
	return self.unit
	
    def getNumericValue(self):
	return self.value
	
    # Contributed by Berthold Hoellmann
    def inBaseUnits(self):
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in xrange(9):
            unit = _base_names[i]
            power = self.unit.powers[i]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return self.__class__(value = new_value, unit = num + denom)

    def isCompatible (self, unit):
        unit = _findUnit (unit)
        return self.unit.isCompatible (unit)

    def sqrt(self):
	return pow(self, 0.5)

    def sin(self):
	if self.unit.isAngle():
	    return umath.sin(self.value * \
			     self.unit.conversionFactorTo(_unit_table['rad']))
	else:
	    raise TypeError, 'Argument of sin must be an angle'

    def cos(self):
	if self.unit.isAngle():
	    return umath.cos(self.value * \
			     self.unit.conversionFactorTo(_unit_table['rad']))
	else:
	    raise TypeError, 'Argument of cos must be an angle'

    def tan(self):
	if self.unit.isAngle():
	    return umath.tan(self.value * \
			     self.unit.conversionFactorTo(_unit_table['rad']))
	else:
	    raise TypeError, 'Argument of tan must be an angle'
	    
    def arctan2(self,other):
	other = self._inMyUnits(other)
	return PhysicalField(value = umath.arctan2(self.value, other.value), unit = "rad")
	    
    def arctan(self):
	return PhysicalField(value = umath.arctan(self.value), unit = "rad")
	    
    def dot(self, other):
	if not isinstance(other,PhysicalField):
	    return self.__class__(value = Numeric.dot(self.value,other), unit = self.unit)
	value = Numeric.dot(self.value,other.value)
	unit = self.unit*other.unit
	if unit.isDimensionless():
	    return value*unit.factor
	else:
	    return self.__class__(value = value, unit = unit)
	    
    def take(self, ids):
	return self.__class__(value = Numeric.take(self.value, ids), unit = self.unit)
	
    def reshape(self, shape):
	return self.__class__(value = Numeric.reshape(self.value, shape), unit = self.unit)
	    
    def sum(self, index = 0):
	return self.__class__(value = Numeric.sum(self.value, index), unit = self.unit)

class PhysicalUnit:

    def __init__(self, names, factor, powers, offset=0):
	if type(names) == type(''):
	    self.names = NumberDict()
	    self.names[names] = 1
	else:
	    self.names = names
	self.factor = factor
	self.offset = offset
	self.powers = powers

    def __repr__(self):
	return '<PhysicalUnit ' + self.name() + '>'

    __str__ = __repr__

    def __cmp__(self, other):
	if self.powers != other.powers:
	    raise TypeError, 'Incompatible units'
	return cmp(self.factor, other.factor)

    def __mul__(self, other):
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot multiply units with non-zero offset"
	if isinstance(other,PhysicalUnit):
	    return PhysicalUnit(self.names+other.names,
				self.factor*other.factor,
				map(lambda a,b: a+b,
				    self.powers, other.powers))
	else:
	    return PhysicalUnit(self.names+{str(other): 1},
				self.factor*other,
				self.powers,
				self.offset * other)

    __rmul__ = __mul__

    def __div__(self, other):
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot divide units with non-zero offset"
	if isinstance(other,PhysicalUnit):
	    return PhysicalUnit(self.names-other.names,
				self.factor/other.factor,
				map(lambda a,b: a-b,
				    self.powers, other.powers))
	else:
	    return PhysicalUnit(self.names+{str(other): -1},
				self.factor/other, self.powers)

    def __rdiv__(self, other):
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot divide units with non-zero offset"
	if isinstance(other,PhysicalUnit):
	    return PhysicalUnit(other.names-self.names,
				other.factor/self.factor,
				map(lambda a,b: a-b,
				    other.powers, self.powers))
	else:
	    return PhysicalUnit({str(other): 1}-self.names,
				other/self.factor,
				map(lambda x: -x, self.powers))

    def __pow__(self, other):
        if self.offset != 0:
            raise TypeError, "cannot exponentiate units with non-zero offset"
	if type(other) == type(0):
	    return PhysicalUnit(other*self.names, pow(self.factor, other),
				map(lambda x,p=other: x*p, self.powers))
	if type(other) == type(0.):
	    inv_exp = 1./other
	    rounded = int(umath.floor(inv_exp+0.5))
	    if abs(inv_exp-rounded) < 1.e-10:
		if reduce(lambda a, b: a and b,
			  map(lambda x, e=rounded: x%e == 0, self.powers)):
		    f = pow(self.factor, other)
		    p = map(lambda x,p=rounded: x/p, self.powers)
		    if reduce(lambda a, b: a and b,
			      map(lambda x, e=rounded: x%e == 0,
				  self.names.values())):
			names = self.names/rounded
		    else:
			names = NumberDict()
			if f != 1.:
			    names[str(f)] = 1
			for i in range(len(p)):
			    names[_base_names[i]] = p[i]
		    return PhysicalUnit(names, f, p)
		else:
		    raise TypeError, 'Illegal exponent'
	raise TypeError, 'Only integer and inverse integer exponents allowed'

    def conversionFactorTo(self, other):
	if self.powers != other.powers:
	    if self.isDimensionlessOrAngle() and other.isDimensionlessOrAngle():
		return self.factor/other.factor
	    else:
		raise TypeError, 'Incompatible units'
        if self.offset != other.offset and self.factor != other.factor:
	    raise TypeError, \
                  ('Unit conversion (%s to %s) cannot be expressed ' +
                   'as a simple multiplicative factor') % \
                  (self.name(), other.name())
	return self.factor/other.factor

    def conversionTupleTo(self, other): # added 1998/09/29 GPW
	if self.powers != other.powers:
	    raise TypeError, 'Incompatible units'

	# let (s1,d1) be the conversion tuple from 'self' to base units
	#   (ie. (x+d1)*s1 converts a value x from 'self' to base units,
	#   and (x/s1)-d1 converts x from base to 'self' units)
	# and (s2,d2) be the conversion tuple from 'other' to base units
	# then we want to compute the conversion tuple (S,D) from
	#   'self' to 'other' such that (x+D)*S converts x from 'self'
	#   units to 'other' units
	# the formula to convert x from 'self' to 'other' units via the
	#   base units is (by definition of the conversion tuples):
	#     ( ((x+d1)*s1) / s2 ) - d2
	#   = ( (x+d1) * s1/s2) - d2
	#   = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
	#   = ( (x+d1) - (d1*s2/s1) ) * s1/s2
	#   = (x + d1 - d2*s2/s1) * s1/s2
	# thus, D = d1 - d2*s2/s1 and S = s1/s2
	factor = self.factor / other.factor
	offset = self.offset - (other.offset * other.factor / self.factor)
	return (factor, offset)

    def isCompatible (self, other):     # added 1998/10/01 GPW
        return self.powers == other.powers

    def isDimensionless(self):
	return not reduce(lambda a,b: a or b, self.powers)

    def isAngle(self):
	return self.powers[7] == 1 and \
	       reduce(lambda a,b: a + b, self.powers) == 1
	       
    def isDimensionlessOrAngle(self):
	return self.isDimensionless() or self.isAngle()
		   
    def setName(self, name):
	self.names = NumberDict()
	self.names[name] = 1

    def name(self):
	num = ''
	denom = ''
	for unit in self.names.keys():
	    power = self.names[unit]
	    if power < 0:
		denom = denom + '/' + unit
		if power < -1:
		    denom = denom + '**' + str(-power)
	    elif power > 0:
		num = num + '*' + unit
		if power > 1:
		    num = num + '**' + str(power)
	if len(num) == 0:
	    num = '1'
	else:
	    num = num[1:]
	return num + denom

# Helper functions

def _findUnit(unit):
    if type(unit) == type(''):
	name = string.strip(unit)
	if len(name) == 0:
	    name = "m/m" # allow dimensionless quantities
	unit = eval(name, _unit_table)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass

    if not isinstance(unit,PhysicalUnit):
	if unit == 1:
	    unit = eval("m/m", _unit_table)
	else:
	    raise TypeError, str(unit) + ' is not a unit'
    return unit

def _round(x):
    if umath.greater(x, 0.):
	return umath.floor(x)
    else:
	return umath.ceil(x)


def _convertValue (value, src_unit, target_unit):
    (factor, offset) = src_unit.conversionTupleTo(target_unit)
    return (value + offset) * factor

def Scale(quantity, scaling):
    """ normalize 'quantity' by 'scaling'.
    
    'quantity' can be a PhysicalQuantity, a value-unit string convertable to
    a PhysicalQuantity, or a dimensionless number. A dimensionless number
    is left alone. It is an error for the result to have dimensions.
    """
    
    quantity = AnyQuantity(quantity)

    if isinstance(quantity,PhysicalField) and not quantity.unit.isDimensionless():
	scaling = AnyQuantity(scaling)
	# normalize quantity to scaling
	# error will be thrown if incompatible
	dimensionless = quantity / scaling
    else:
	# Assume quantity is a dimensionless number and return it.
	# Automatically throws an error if it's not a number.
	dimensionless = quantity
		
    if isinstance(dimensionless,PhysicalField) and not dimensionless.unit.isDimensionless():
	raise TypeError, `quantity.inBaseUnits().unit` + ' and ' \
	+ `scaling.inBaseUnits().unit` \
	+ ' are incompatible'
	
    return dimensionless
 

def NonDimOrUnits(quantity, units):
    quantity = AnyQuantity(quantity)
    if isinstance(quantity,PhysicalField) and not quantity.unit.isDimensionless():
	quantity.convertToUnit(units)
    return quantity
	
def AnyQuantity(quantity):
    """ normalize 'quantity'.
    
    'quantity' can be a PhysicalField, a value-unit string convertable to
    a PhysicalField, or a dimensionless number.
    """
    
    if type(quantity) == type(''):
	# input is string, so construct a PhysicalQuantity
	quantity = PhysicalField(quantity)
	
    if not isinstance(quantity,PhysicalField):
	# Assume quantity is a dimensionless number and return it.
	# Automatically throws an error if it's not a number.
	quantity = float(quantity)
	
    return quantity

def AddConstant(name, constant):
##     if _unit_table.has_key(name):
## 	raise KeyError, 'Constant ' + name + ' already defined'
    if isinstance(constant,PhysicalField):
	constant = PhysicalUnit(constant.unit.names, constant.value, constant.unit.powers, constant.unit.offset)
    elif type(constant) == type(''):
	constant = eval(constant, _unit_table)
	for cruft in ['__builtins__', '__args__']:
	    try: del _unit_table[cruft]
	    except: pass
    if isinstance(constant,PhysicalUnit):
	constant.setName(name)
    _unit_table[name] = constant

def _isVariable(var):
    return isinstance(var,fivol.variables.variable.Variable)
    
# SI unit definitions

_base_names = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [('m',   PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0])),
	       ('g',   PhysicalUnit('g',   0.001, [0,1,0,0,0,0,0,0,0])),
	       ('s',   PhysicalUnit('s',   1.,    [0,0,1,0,0,0,0,0,0])),
	       ('A',   PhysicalUnit('A',   1.,    [0,0,0,1,0,0,0,0,0])),
	       ('K',   PhysicalUnit('K',   1.,    [0,0,0,0,1,0,0,0,0])),
	       ('mol', PhysicalUnit('mol', 1.,    [0,0,0,0,0,1,0,0,0])),
	       ('cd',  PhysicalUnit('cd',  1.,    [0,0,0,0,0,0,1,0,0])),
	       ('rad', PhysicalUnit('rad', 1.,    [0,0,0,0,0,0,0,1,0])),
	       ('sr',  PhysicalUnit('sr',  1.,    [0,0,0,0,0,0,0,0,1])),
	       ]

_prefixes = [('Y',  1.e24),
             ('Z',  1.e21),
             ('E',  1.e18),
	     ('P',  1.e15),
	     ('T',  1.e12),
	     ('G',  1.e9),
	     ('M',  1.e6),
	     ('k',  1.e3),
	     ('h',  1.e2),
	     ('da', 1.e1),
	     ('d',  1.e-1),
	     ('c',  1.e-2),
	     ('m',  1.e-3),
	     ('mu', 1.e-6),
	     ('n',  1.e-9),
	     ('p',  1.e-12),
	     ('f',  1.e-15),
	     ('a',  1.e-18),
	     ('z',  1.e-21),
	     ('y',  1.e-24),
	     ]

_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]

def _addUnit(name, unit):
    if _unit_table.has_key(name):
	raise KeyError, 'Unit ' + name + ' already defined'
    if type(unit) == type(''):
	unit = eval(unit, _unit_table)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass
    unit.setName(name)
    _unit_table[name] = unit

def _addPrefixed(unit):
    for prefix in _prefixes:
	name = prefix[0] + unit
	_addUnit(name, prefix[1]*_unit_table[unit])


# SI derived units; these automatically get prefixes

_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])

_addUnit('Hz', '1/s')                # Hertz
_addUnit('N', 'm*kg/s**2')           # Newton
_addUnit('Pa', 'N/m**2')             # Pascal
_addUnit('J', 'N*m')                 # Joule
_addUnit('W', 'J/s')                 # Watt
_addUnit('C', 's*A')                 # Coulomb
_addUnit('V', 'W/A')                 # Volt
_addUnit('F', 'C/V')                 # Farad
_addUnit('ohm', 'V/A')               # Ohm
_addUnit('S', 'A/V')                 # Siemens
_addUnit('Wb', 'V*s')                # Weber
_addUnit('T', 'Wb/m**2')             # Tesla
_addUnit('H', 'Wb/A')                # Henry
_addUnit('lm', 'cd*sr')              # Lumen
_addUnit('lx', 'lm/m**2')            # Lux
_addUnit('Bq', '1/s')                # Becquerel
_addUnit('Gy', 'J/kg')               # Gray
_addUnit('Sv', 'J/kg')               # Sievert

del _unit_table['kg']

for unit in _unit_table.keys():
    _addPrefixed(unit)

# Fundamental constants

_unit_table['pi'] = umath.pi
_addUnit('c', '299792458.*m/s')      # speed of light
_addUnit('mu0', '4.e-7*pi*N/A**2')   # permeability of vacuum
_addUnit('eps0', '1/mu0/c**2')       # permittivity of vacuum
_addUnit('Grav', '6.67259e-11*m**3/kg/s**2') # gravitational constant
_addUnit('hplanck', '6.6260755e-34*J*s')     # Planck constant
_addUnit('hbar', 'hplanck/(2*pi)')   # Planck constant / 2pi
_addUnit('e', '1.60217733e-19*C')    # elementary charge
_addUnit('me', '9.1093897e-31*kg')   # electron mass
_addUnit('mp', '1.6726231e-27*kg')   # proton mass
_addUnit('Nav', '6.0221367e23/mol')  # Avogadro number
_addUnit('kB', '1.380658e-23*J/K')    # Boltzmann constant

# Time units

_addUnit('min', '60*s')              # minute
_addUnit('h', '60*min')              # hour
_addUnit('d', '24*h')                # day
_addUnit('wk', '7*d')                # week
_addUnit('yr', '365.25*d')           # year

# Length units

_addUnit('inch', '2.54*cm')          # inch
_addUnit('ft', '12*inch')            # foot
_addUnit('yd', '3*ft')               # yard
_addUnit('mi', '5280.*ft')           # (British) mile
_addUnit('nmi', '1852.*m')           # Nautical mile
_addUnit('Ang', '1.e-10*m')          # Angstrom
_addUnit('lyr', 'c*yr')              # light year
_addUnit('Bohr', '4*pi*eps0*hbar**2/me/e**2')  # Bohr radius

# Area units

_addUnit('ha', '10000*m**2')         # hectare
_addUnit('acres', 'mi**2/640')       # acre
_addUnit('b', '1.e-28*m')            # barn

# Volume units

_addUnit('l', 'dm**3')               # liter
_addUnit('dl', '0.1*l')
_addUnit('cl', '0.01*l')
_addUnit('ml', '0.001*l')
_addUnit('tsp', '4.92892159375*ml')  # teaspoon
_addUnit('tbsp', '3*tsp')            # tablespoon
_addUnit('floz', '2*tbsp')           # fluid ounce
_addUnit('cup', '8*floz')            # cup
_addUnit('pt', '16*floz')            # pint
_addUnit('qt', '2*pt')               # quart
_addUnit('galUS', '4*qt')            # US gallon
_addUnit('galUK', '4.54609*l')       # British gallon

# Mass units

_addUnit('amu', '1.6605402e-27*kg')  # atomic mass units
_addUnit('oz', '28.349523125*g')     # ounce
_addUnit('lb', '16*oz')              # pound
_addUnit('ton', '2000*lb')           # ton

# Force units

_addUnit('dyn', '1.e-5*N')           # dyne (cgs unit)

# Energy units

_addUnit('erg', '1.e-7*J')           # erg (cgs unit)
_addUnit('eV', 'e*V')                # electron volt
_addPrefixed('eV')
_addUnit('Hartree', 'me*e**4/16/pi**2/eps0**2/hbar**2')
_addUnit('invcm', 'hplanck*c/cm')    # Wavenumbers/inverse cm
_addUnit('Ken', 'kB*K')               # Kelvin as energy unit
_addUnit('cal', '4.184*J')           # thermochemical calorie
_addUnit('kcal', '1000*cal')         # thermochemical kilocalorie
_addUnit('cali', '4.1868*J')         # international calorie
_addUnit('kcali', '1000*cali')       # international kilocalorie
_addUnit('Btu', '1055.05585262*J')   # British thermal unit

# Power units

_addUnit('hp', '745.7*W')            # horsepower

# Pressure units

_addUnit('bar', '1.e5*Pa')           # bar (cgs unit)
_addUnit('atm', '101325.*Pa')        # standard atmosphere
_addUnit('torr', 'atm/760')          # torr = mm of mercury
_addUnit('psi', '6894.75729317*Pa')  # pounds per square inch

# Angle units

_addUnit('deg', 'pi*rad/180')        # degrees

# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
kelvin = _findUnit ('K')
_addUnit ('degR', '(5./9.)*K')       # degrees Rankine
_addUnit ('degC', PhysicalUnit (None, 1.0, kelvin.powers, 273.15))
_addUnit ('degF', PhysicalUnit (None, 5./9., kelvin.powers, 459.67))
del kelvin


# Some demonstration code. Run with "python -i PhysicalQuantities.py"
# to have this available.

if __name__ == '__main__':

    from umath import *
    l = PhysicalField(10., 'm')
    big_l = PhysicalField(10., 'km')
    print big_l + l
    t = PhysicalField(314159., 's')
    print t.inUnitsOf('d','h','min','s')

    p = PhysicalField # just a shorthand...

    e = p('2.7 Hartree*Nav')
    e.convertToUnit('kcal/mol')
    print e
    print e.inBaseUnits()

    freeze = p('0 degC')
    print freeze.inUnitsOf ('degF')
    
    print PhysicalField("1")
    print PhysicalField(2., " ")
    print PhysicalField(2.)
    
    import Numeric
    a = PhysicalField(Numeric.array(((3.,4.),(5.,6.))),"m")
    print a
    a[0,1] = PhysicalField("6 ft")
    print a
    print a > "13 ft"
    print a > PhysicalField("13 ft")
    print a > PhysicalField(Numeric.array(((3.,13.),(17.,6.))),"ft")
    print a > PhysicalField("1 lb")

