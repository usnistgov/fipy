#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "array.py"
 #                                    created: 1/10/04 {10:23:17 AM} 
 #                                last update: 6/3/05 {12:01:27 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

"""

The functions provided in ths module replace the `Numeric` module.
The functions work with `Variables`, arrays or numbers. For example,
create a `Variable`.

   >>> from fipy.variables.variable import Variable
   >>> var = Variable(value = 0)

Take the tangent of such a variable. The returned value is itself a
`Variable`.

   >>> v = tan(var)
   >>> v
   array.tan(Variable(value = 0))
   >>> print v
   0.0

Take the tangent of a int.

   >>> tan(0)
   0.0
   
Take the tangent of an array.

   >>> tan(Numeric.array((0,0,0)))
   [ 0., 0., 0.,]
   
This module is building towards a '`Numerix`' module that will be the
only place in the code where `Numeric` is imported.

"""

__docformat__ = 'restructuredtext'

import Numeric
import umath
import MA
import fipy.tools.inline.inline as inline

def _isPhysical(arr):
    """
    Returns `True` if arr is a `Variable` or `PhysicalField`.
    """
    import fipy.variables.variable
    import fipy.tools.dimensions.physicalField

    return isinstance(arr,fipy.variables.variable.Variable) \
	or isinstance(arr,fipy.tools.dimensions.physicalField.PhysicalField)

def take(arr, ids):
    """
    Provides the same functionality as `Numeric.take`.
    """
    
    if _isPhysical(arr):
	return arr.take(ids)    
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.take(arr, ids)
    elif type(arr) is type(MA.array((0))):
	return MA.take(arr, ids)
    else:
	raise TypeError, 'cannot take from object ' + str(arr)
    
def put(arr, ids, values):
    """
    Provides the same functionality as `Numeric.put`.
    """
    if _isPhysical(arr):
	return arr.put(ids, values)
    
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.put(arr, ids, values)
    elif type(arr) is type(MA.array((0))):
	return MA.put(arr, ids, values)
    else:
	raise TypeError, 'cannot put in object ' + str(arr)
    
def reshape(arr, shape):
    """
    Provides the same functionality as `Numeric.reshape`.
    """
    if _isPhysical(arr):
	return arr.reshape(shape)
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.reshape(arr, shape)
    elif type(arr) is type(MA.array((0))):
        return MA.reshape(arr, shape)
    else:
	raise TypeError, 'cannot reshape object ' + str(arr)
	
def sum(arr, index = 0):
    """
    Provides the same functionality as `Numeric.sum`.
    """
    if _isPhysical(arr):
	return arr.sum(index)
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.sum(arr, index)
    elif type(arr) is type(MA.array((0))):
	return MA.sum(arr, index)
    else:        
	raise TypeError, 'cannot sum object ' + str(arr)

def sqrt(arr):
    """
    Provides the same functionality as `Numeric.sqrt`.
    """
    if _isPhysical(arr):
	return arr.sqrt()
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.sqrt(arr)
    else:
	return umath.sqrt(arr)

def exp(arr):
    """
    Provides the same functionality as `Numeric.exp`.
    """
    if _isPhysical(arr):
	return arr.exp()
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.exp(arr)
    else:
	return umath.exp(arr)
	
def tan(arr):
    """
    Provides the same functionality as `Numeric.tan`.
    """
    if _isPhysical(arr):
	return arr.tan()
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.tan(arr)
    else:
	return umath.tan(arr)

def arctan(arr):
    """
    Provides the same functionality as `Numeric.arctan`.
    """
    if _isPhysical(arr):
	return arr.arctan()
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.arctan(arr)
    else:
	return umath.arctan(arr)
		
def arctan2(arr, other):
    """
    Provides the same functionality as `Numeric.arctan2`.
    """
    if _isPhysical(arr):
	return arr.arctan2(other)
    elif _isPhysical(other):
	import fipy.tools.dimensions.physicalField

	return physicalField.PhysicalField(value = arr, unit = "rad").arctan2(other)
    elif type(arr) is type(Numeric.array((0))):
	return Numeric.arctan2(arr,other)
    else:
	return umath.arctan2(arr,other)
	
def crossProd(v1,v2):
    """Return vector cross-product of v1 and v2.
    """
    v1n = Numeric.reshape(v1, (-1, 3))
    v2n = Numeric.reshape(v2, (-1, 3))

    out = Numeric.transpose(Numeric.array((v1n[:,1] * v2n[:,2] - v1n[:,2] * v2n[:,1],
			    v1n[:,2] * v2n[:,0] - v1n[:,0] * v2n[:,2],
			    v1n[:,0] * v2n[:,1] - v1n[:,1] * v2n[:,0])))
    return Numeric.reshape(out, Numeric.shape(v1))

def dot(a1, a2, axis = 1):
    """
    return array of vector dot-products of v1 and v2
    for arrays a1 and a2 of vectors v1 and v2
    
    We can't use Numeric.dot on an array of vectors
    """
    return sum((a1*a2)[:], axis)

def sqrtDot(a1, a2):
    """Return array of square roots of vector dot-products
    for arrays a1 and a2 of vectors v1 and v2
    
    Usually used with v1==v2 to return magnitude of v1.
    """
    ## We can't use Numeric.dot on an array of vectors
##     return Numeric.sqrt(Numeric.sum((a1*a2)[:],1))
##    return fipy.tools.array.sqrt(fipy.tools.array.sum((a1*a2)[:],1))
    return inline._optionalInline(_sqrtDotIn, _sqrtDotPy, a1, a2)

def _sqrtDotPy(a1, a2):
    return sqrt(dot(a1, a2))

##def _sqrtDotIn(a1, a2):
##    ni, nj = Numeric.shape(a1)
##    result = Numeric.zeros((ni,),'d')
##    inline.runInlineLoop1("""
##	int j;
##	result(i) = 0.;
##	for (j = 0; j < nj; j++)
##	{
##	    result(i) += a1(i,j) * a2(i,j);
##	}
##	result(i) = sqrt(result(i));
##    """,result = result, a1 = a1, a2 = a2, ni = ni, nj = nj) 
##    return result

def _sqrtDotIn(a1, a2):
    unit1 = unit2 = 1
    if _isPhysical(a1):
	unit1 = a1.inBaseUnits().getUnit()
        a1 = a1.getNumericValue()
    if _isPhysical(a2):
	unit2 = a2.inBaseUnits().getUnit()
        a2 = a2.getNumericValue()
    ni, nj = Numeric.shape(a1)
    result = Numeric.zeros((ni,),'d')
    inline._runInline("""
        int i;
        for (i = 0; i < ni; i++)
	{
	    int j;
            result(i) = 0.;
            for (j = 0; j < nj; j++)
            {
	        result(i) += a1(i,j) * a2(i,j);
            }
            result(i) = sqrt(result(i));
        }
    """,result = result, a1 = a1, a2 = a2, ni = ni, nj = nj)
    if unit1 != 1 or unit2 != 1:
	from fipy.tools.dimensions.physicalField import PhysicalField
	result = PhysicalField(value = result, unit = (unit1 * unit2)**0.5)
    return result

def allequal(first, second):
    """
    Provides the same functionality as `Numeric.allequal`.
    """
    if _isPhysical(first):
        return first.allequal(second)
## 	return MA.alltrue(first == second)
    elif _isPhysical(second):
        return first.allequal(first)
## 	return MA.alltrue(second == first)
    else:
	return MA.allequal(first, second)
	    
def allclose(first, second, rtol = 1.e-5, atol = 1.e-8):
    """
    Provides the same functionality as `Numeric.allclose`.
    """
    if _isPhysical(first):
	return first.allclose(other = second, atol = atol, rtol = rtol)
    elif _isPhysical(second):
	return second.allclose(other = first, atol = atol, rtol = rtol)
    else:
	return MA.allclose(first, second, atol = atol, rtol = rtol)

def _min(arr):
    arr = Numeric.array(arr)
    return arr[Numeric.argmin(arr)]

def min(arr):
    """
    Find the minimum value in `arr`.
    """
    if _isPhysical(arr):
        return arr.min()
    else:
        return _min(arr)

def _max(arr):
    arr = Numeric.array(arr)
    return arr[Numeric.argmax(arr)]

def max(arr):
    """
    Find the maximum value in `arr`.
    """
    if _isPhysical(arr):
        return arr.max()
    else:
        return _max(arr)

        
# Necessary because LLNL hires stupidheads
def MAtake(array, indices, fill = 0, axis = 0):
    """
    Replaces `MA.take`. `MA.take` does not always work when
    `indices` is a masked array.
    """
    tmp = MA.take(array, MA.filled(indices, fill), axis = axis)
    if indices.mask() is not None and tmp.shape != indices.mask().shape:
        mask = MA.repeat(indices.mask()[...,Numeric.NewAxis],tmp.shape[-1],len(tmp.shape)-1)
        if tmp.mask() is not None:
            mask = Numeric.logical_or(tmp.mask(), mask)
    else:
        mask = indices.mask()
    return MA.array(data = tmp, mask = mask)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__":
    _test() 
