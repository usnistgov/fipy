#!/usr/bin/env python

#-----*-Pyth-*-
####################################################################
# ElPhF - electrochemistry phase field modeling
#
# FILE: "Dimensionalization.py"
#                                   created: 8/6/03 {9:40:12 AM} 
#                               last update: 12/12/03 {5:08:53 PM} 
# Author: Jonathan Guyer
# E-mail: guyer@nist.gov
#   mail: National Institute of Standards and Technology, Metallurgy Division
#         100 Bureau Dr. - Stop 8551, Gaithersburg, MD 20899-8551
#    www: http://www.metallurgy.nist.gov
# 
#========================================================================
#This software was developed at the National Institute of Standards
#and Technology by employees of the Federal Government in the course
#of their official duties.  Pursuant to title 17 Section 105 of the
#United States Code this software is not subject to copyright
#protection and is in the public domain.  ElPhF is an experimental
#system.  NIST assumes no responsibility whatsoever for its use by
#other parties, and makes no guarantees, expressed or implied, about
#its quality, reliability, or any other characteristic.  We would
#appreciate acknowledgement if the software is used.
#
#This software can be redistributed and/or modified freely
#provided that any derivative works bear some notice that they are
#derived from it, and any modified versions bear some notice that
#they have been modified.
#========================================================================
# Description: 
#
# History
#
# modified   by  rev reason
# ---------- --- --- -----------
# 2003-08-06 JEG 1.0 original
####################################################################
#----

from Scientific.Physics.PhysicalQuantities import PhysicalQuantity
from Scientific.Physics.PhysicalQuantities import isPhysicalQuantity

class PhysicalField(PhysicalQuantity):
    def _inMyUnits(self, other):
	if not isPhysicalQuantity(other):
	    if type(other) is type(''):
		other = PhysicalField(other)
	    else:
		raise TypeError, 'Incompatible types'
	return other.inUnitsOf(self.unit)
	
    def __getitem__(self, index): 
	return PhysicalField(self.value[index],self.unit)
	
    def __setitem__(self, index, value):
	if type(value) is type(''):
	    value = PhysicalField(value)
	if isPhysicalQuantity(value):
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
	    

    
def nonDimensionalize(quantity, scaling):
	""" normalize 'quantity' by 'scaling'.
	
	'quantity' can be a PhysicalQuantity, a value-unit string convertable to
	a PhysicalQuantity, or a dimensionless number. A dimensionless number
	is left alone. It is an error for the result to have dimensions.
	"""
	
	quantity = anyQuantity(quantity)
	scaling = anyQuantity(scaling)

	if isPhysicalQuantity(quantity):
		# normalize quantity to scaling
		# error will be thrown if incompatible
		dimensionless = quantity / scaling
	else:
		# Assume quantity is a dimensionless number and return it.
		# Automatically throws an error if it's not a number.
		dimensionless = quantity
		
	if isPhysicalQuantity(dimensionless):
		raise TypeError, `quantity.inBaseUnits().unit` + ' and ' \
				+ `scaling.inBaseUnits().unit` \
				+ ' are incompatible'
		
	return dimensionless
 
def dimensionlessOrInUnits(quantity, units):
	quantity = anyQuantity(quantity)
	if isPhysicalQuantity(quantity):
		quantity.convertToUnit(units)
	return quantity
	
def anyQuantity(quantity):
	""" normalize 'quantity'.
	
	'quantity' can be a PhysicalQuantity, a value-unit string convertable to
	a PhysicalQuantity, or a dimensionless number.
	"""
	
	if type(quantity) == type(''):
		# input is string, so construct a PhysicalQuantity
		quantity = PhysicalQuantity(quantity)
		
	if not isPhysicalQuantity(quantity):
		# Assume quantity is a dimensionless number and return it.
		# Automatically throws an error if it's not a number.
		quantity = float(quantity)
		
	return quantity

if __name__ == '__main__':
    import Numeric
    a = PhysicalField(Numeric.array(((3.,4.),(5.,6.))),"m")
    print a
    a[0,1] = PhysicalQuantity("6 ft")
    print a
    print a > "13 ft"
    print a > PhysicalQuantity("13 ft")
    print a > PhysicalField(Numeric.array(((3.,13.),(17.,6.))),"ft")
    print a > PhysicalQuantity("1 lb")
# 	print 3, anyQuantity(3)
# 	print "3 nm", anyQuantity("3 nm")
# 	print PhysicalQuantity("3 nm"), anyQuantity(PhysicalQuantity("3 nm"))
# 	print "3 ft in m is", nonDimensionalOrInUnits("3 ft", "m")
# 	print "3 ft scaled to 1 m is", nonDimensionalize("3 ft", "1 m")
# 	print [3], anyQuantity([3])	