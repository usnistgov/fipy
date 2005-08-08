#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "numerix.py"
 #                                    created: 1/10/04 {10:23:17 AM} 
 #                                last update: 8/8/05 {10:18:15 AM} 
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
   numerix.tan(Variable(value = 0))
   >>> print v
   0.0

Take the tangent of a int.

   >>> tan(0)
   0.0
   
Take the tangent of an array.

   >>> tan(array((0,0,0)))
   [ 0., 0., 0.,]
   
This module is building towards a '`Numerix`' module that will be the
only place in the code where `Numeric` is imported.

"""

__docformat__ = 'restructuredtext'

import Numeric as NUMERIC
from Numeric import *

## import numarray as NUMERIC
## from numarray import *

## from Numeric import pi, array, zeros, ones, ravel, concatenate, where, arange


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
    elif type(arr) is type(array((0))):
	return NUMERIC.take(arr, ids)
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
    
    elif type(arr) is type(array((0))):
	return NUMERIC.put(arr, ids, values)
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
    elif type(arr) is type(array((0))):
	return NUMERIC.reshape(arr, tuple(shape))
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
    elif type(arr) is type(array((0))):
	return NUMERIC.sum(arr, index)
    elif type(arr) is type(MA.array((0))):
	return MA.sum(arr, index)
    else:        
	raise TypeError, 'cannot sum object ' + str(arr)

################################################
#                                              #
#   mathematical and trigonometric functions   #
#                                              #
################################################

def arccos(arr):
    r"""
    Inverse cosine of
    
    .. raw:: latex
    
       $x$, $\cos^{-1} x$
       
    ..
    
        >>> print arccos(0.0)
        1.57079632679
        >>> print arccos(2.0)
        Traceback (most recent call last):
            ...
        OverflowError: math range error
        >>> print arccos(array((0,0.5,1.0)))
        [ 1.57079633, 1.04719755, 0.        ,]
        >>> from fipy.variables.variable import Variable
        >>> arccos(Variable(value = (0,0.5,1.0)))
        numerix.arccos(Variable(value = [ 0. , 0.5, 1. ,]))
        >>> print arccos(Variable(value = (0,0.5,1.0)))
        [ 1.57079633, 1.04719755, 0.        ,] rad
    """
    if _isPhysical(arr):
        return arr.arccos()
    elif type(arr) is type(array((0))):
        return NUMERIC.arccos(arr)
    else:
        return umath.arccos(arr)

def arccosh(arr):
    r"""
    Inverse hyperbolic cosine of
    
    .. raw:: latex
    
       $x$, $\cosh^{-1} x$
       
    ..
    
        >>> print arccosh(1.0)
        0.0
        >>> print arccosh(0.0)
        Traceback (most recent call last):
            ...
        OverflowError: math range error
        >>> print arccosh(array((1,2,3)))
        [ 0.        , 1.3169579 , 1.76274717,]
        >>> from fipy.variables.variable import Variable
        >>> arccosh(Variable(value = (1,2,3)))
        numerix.arccosh(Variable(value = [ 1., 2., 3.,]))
        >>> print arccosh(Variable(value = (1,2,3)))
        [ 0.        , 1.3169579 , 1.76274717,]
    """
    if _isPhysical(arr):
        return arr.arccosh()
    elif type(arr) is type(array((0))):
        return NUMERIC.arccosh(arr)
    else:
        return umath.arccosh(arr)

def arcsin(arr):
    r"""
    Inverse sine of
    
    .. raw:: latex
    
       $x$, $\sin^{-1} x$
       
    ..
    
        >>> print arcsin(1.0)
        1.57079632679
        >>> print arcsin(2.0)
        Traceback (most recent call last):
            ...
        OverflowError: math range error
        >>> print arcsin(array((0,0.5,1.0)))
        [ 0.        , 0.52359878, 1.57079633,]
        >>> from fipy.variables.variable import Variable
        >>> arcsin(Variable(value = (0,0.5,1.0)))
        numerix.arcsin(Variable(value = [ 0. , 0.5, 1. ,]))
        >>> print arcsin(Variable(value = (0,0.5,1.0)))
        [ 0.        , 0.52359878, 1.57079633,] rad
    """
    if _isPhysical(arr):
        return arr.arcsin()
    elif type(arr) is type(array((0))):
        return NUMERIC.arcsin(arr)
    else:
        return umath.arcsin(arr)

def arcsinh(arr):
    r"""
    Inverse hyperbolic sine of
    
    .. raw:: latex
    
       $x$, $\sinh^{-1} x$
       
    ..

        >>> print arcsinh(1.0)
        0.88137358702
        >>> print arcsinh(array((1,2,3)))
        [ 0.88137359, 1.44363548, 1.81844646,]
        >>> from fipy.variables.variable import Variable
        >>> arcsinh(Variable(value = (1,2,3)))
        numerix.arcsinh(Variable(value = [ 1., 2., 3.,]))
        >>> print arcsinh(Variable(value = (1,2,3)))
        [ 0.88137359, 1.44363548, 1.81844646,]
    """
    if _isPhysical(arr):
        return arr.arcsinh()
    elif type(arr) is type(array((0))):
        return NUMERIC.arcsinh(arr)
    else:
        return umath.arcsinh(arr)

def arctan(arr):
    r"""
    Inverse tangent of
    
    .. raw:: latex
    
       $x$, $\tan^{-1} x$
       
    ..
    
        >>> print arctan(1.0)
        0.785398163397
        >>> print arctan(array((0,0.5,1.0)))
        [ 0.        , 0.46364761, 0.78539816,]
        >>> from fipy.variables.variable import Variable
        >>> arctan(Variable(value = (0,0.5,1.0)))
        numerix.arctan(Variable(value = [ 0. , 0.5, 1. ,]))
        >>> print arctan(Variable(value = (0,0.5,1.0)))
        [ 0.        , 0.46364761, 0.78539816,] rad
    """
    if _isPhysical(arr):
        return arr.arctan()
    elif type(arr) is type(array((0))):
        return NUMERIC.arctan(arr)
    else:
        return umath.arctan(arr)
                
def arctan2(arr, other):
    r"""
    Inverse tangent of a ratio
    
    .. raw:: latex
    
       $x/y$, $\tan^{-1} \frac{x}{y}$
       
    ..

        >>> print arctan2(3.0, 3.0)
        0.785398163397
        >>> print arctan2(array((0, 1, 2)), 2)
        [ 0.        , 0.46364761, 0.78539816,]
        >>> from fipy.variables.variable import Variable
        >>> arctan2(Variable(value = (0, 1, 2)), 2)
        (numerix.arctan2(Variable(value = [ 0., 1., 2.,]), 2))
        >>> print arctan2(Variable(value = (0, 1, 2)), 2)
        [ 0.        , 0.46364761, 0.78539816,] rad
    """
    if _isPhysical(arr):
        return arr.arctan2(other)
    elif _isPhysical(other):
        from fipy.tools.dimensions import physicalField

        return physicalField.PhysicalField(value = arr, unit = "rad").arctan2(other)
    elif type(arr) is type(array((0))):
        return NUMERIC.arctan2(arr,other)
    else:
        return umath.arctan2(arr,other)
        
        
def arctanh(arr):
    r"""
    Inverse hyperbolic tangent of
    
    .. raw:: latex
    
       $x$, $\tanh^{-1} x$
       
    ..
    
        >>> print arctanh(0.5)
        0.549306144334
        >>> print arctanh(array((0,0.25,0.5)))
        [ 0.        , 0.25541281, 0.54930614,]
        >>> from fipy.variables.variable import Variable
        >>> arctanh(Variable(value = (0,0.25,0.5)))
        numerix.arctanh(Variable(value = [ 0.  , 0.25, 0.5 ,]))
        >>> print arctanh(Variable(value = (0,0.25,0.5)))
        [ 0.        , 0.25541281, 0.54930614,]
    """
    if _isPhysical(arr):
        return arr.arctanh()
    elif type(arr) is type(array((0))):
        return NUMERIC.arctanh(arr)
    else:
        return umath.arctanh(arr)
        
def cos(arr):
    r"""
    Cosine of
    
    .. raw:: latex
    
       $x$, $\cos x$
       
    ..

        >>> print cos(2*pi/6)
        0.5
        >>> print cos(array((0,2*pi/6,pi/2)))
        [  1.00000000e+00,  5.00000000e-01,  6.12323400e-17,]
        >>> from fipy.variables.variable import Variable
        >>> cos(Variable(value = (0,2*pi/6,pi/2), unit = "rad"))
        numerix.cos(Variable(value = PhysicalField([ 0.        , 1.04719755, 1.57079633,],'rad')))
        >>> print cos(Variable(value = (0,2*pi/6,pi/2), unit = "rad"))
        [  1.00000000e+00,  5.00000000e-01,  6.12323400e-17,]
    """
    if _isPhysical(arr):
        return arr.cos()
    elif type(arr) is type(array((0))):
        return NUMERIC.cos(arr)
    else:
        return umath.cos(arr)

def cosh(arr):
    r"""
    Hyperbolic cosine of
    
    .. raw:: latex
    
       $x$, $\cosh x$
       
    ..

        >>> print cosh(0)
        1.0
        >>> print cosh(array((0,1,2)))
        [ 1.        , 1.54308063, 3.76219569,]
        >>> from fipy.variables.variable import Variable
        >>> cosh(Variable(value = (0,1,2)))
        numerix.cosh(Variable(value = [ 0., 1., 2.,]))
        >>> print cosh(Variable(value = (0,1,2)))
        [ 1.        , 1.54308063, 3.76219569,]
    """
    if _isPhysical(arr):
        return arr.cosh()
    elif type(arr) is type(array((0))):
        return NUMERIC.cosh(arr)
    else:
        return umath.cosh(arr)

def tan(arr):
    r"""
    Tangent of
    
    .. raw:: latex
    
       $x$, $\tan x$
       
    ..

        >>> print tan(pi/3)
        1.73205080757
        >>> print tan(array((0,pi/3,2*pi/3)))
        [ 0.        , 1.73205081,-1.73205081,]
        >>> from fipy.variables.variable import Variable
        >>> tan(Variable(value = (0,pi/3,2*pi/3), unit = "rad"))
        numerix.tan(Variable(value = PhysicalField([ 0.        , 1.04719755, 2.0943951 ,],'rad')))
        >>> print tan(Variable(value = (0,pi/3,2*pi/3), unit = "rad"))
        [ 0.        , 1.73205081,-1.73205081,]
    """
    if _isPhysical(arr):
        return arr.tan()
    elif type(arr) is type(array((0))):
        return NUMERIC.tan(arr)
    else:
        return umath.tan(arr)

def tanh(arr):
    r"""
    Hyperbolic tangent of
    
    .. raw:: latex
    
       $x$, $\tanh x$
       
    ..

        >>> print tanh(1)
        0.761594155956
        >>> print tanh(array((0,1,2)))
        [ 0.        , 0.76159416, 0.96402758,]
        >>> from fipy.variables.variable import Variable
        >>> tanh(Variable(value = (0,1,2)))
        numerix.tanh(Variable(value = [ 0., 1., 2.,]))
        >>> print tanh(Variable(value = (0,1,2)))
        [ 0.        , 0.76159416, 0.96402758,]
    """
    if _isPhysical(arr):
        return arr.tanh()
    elif type(arr) is type(array((0))):
        return NUMERIC.tanh(arr)
    else:
        return umath.tanh(arr)

def log10(arr):
    r"""
    Base-10 logarithm of
    
    .. raw:: latex
    
       $x$, $\log_{10} x$
       
    ..

        >>> print log10(10)
        1.0
        >>> print log10(array((0.1,1,10)))
        [-1., 0., 1.,]
        >>> from fipy.variables.variable import Variable
        >>> log10(Variable(value = (0.1,1,10)))
        numerix.log10(Variable(value = [  0.1,  1. , 10. ,]))
        >>> print log10(Variable(value = (0.1,1,10)))
        [-1., 0., 1.,]
    """
    if _isPhysical(arr):
        return arr.log10()
    elif type(arr) is type(array((0))):
        return NUMERIC.log10(arr)
    else:
        return umath.log10(arr)

def sin(arr):
    r"""
    Sine of
    
    .. raw:: latex
    
       $x$, $\sin x$
       
    ..

        >>> print sin(pi/6)
        0.5
        >>> print sin(array((0,pi/6,pi/2)))
        [ 0. , 0.5, 1. ,]
        >>> from fipy.variables.variable import Variable
        >>> sin(Variable(value = (0,pi/6,pi/2), unit = "rad"))
        numerix.sin(Variable(value = PhysicalField([ 0.        , 0.52359878, 1.57079633,],'rad')))
        >>> print sin(Variable(value = (0,pi/6,pi/2), unit = "rad"))
        [ 0. , 0.5, 1. ,]
    """
    if _isPhysical(arr):
        return arr.sin()
    elif type(arr) is type(array((0))):
        return NUMERIC.sin(arr)
    else:
        return umath.sin(arr)

def sinh(arr):
    r"""
    Hyperbolic sine of
    
    .. raw:: latex
    
       $x$, $\sinh x$
       
    ..

        >>> print sinh(0)
        0.0
        >>> print sinh(array((0,1,2)))
        [ 0.        , 1.17520119, 3.62686041,]
        >>> from fipy.variables.variable import Variable
        >>> sinh(Variable(value = (0,1,2)))
        numerix.sinh(Variable(value = [ 0., 1., 2.,]))
        >>> print sinh(Variable(value = (0,1,2)))
        [ 0.        , 1.17520119, 3.62686041,]
    """
    if _isPhysical(arr):
        return arr.sinh()
    elif type(arr) is type(array((0))):
        return NUMERIC.sinh(arr)
    else:
        return umath.sinh(arr)

def sqrt(arr):
    r"""
    Square root of
    
    .. raw:: latex
    
       $x$, $\sqrt{x}$
       
    ..

        >>> print sqrt(2)
        1.41421356237
        >>> print sqrt(array((1,2,3)))
        [ 1.        , 1.41421356, 1.73205081,]
        >>> from fipy.variables.variable import Variable
        >>> sqrt(Variable(value = (1, 2, 3), unit = "m**2"))
        numerix.sqrt(Variable(value = PhysicalField([ 1., 2., 3.,],'m**2')))
        >>> print sqrt(Variable(value = (1, 2, 3), unit = "m**2"))
        [ 1.        , 1.41421356, 1.73205081,] m

    """
    if _isPhysical(arr):
        return arr.sqrt()
    elif type(arr) is type(array((0))):
        return NUMERIC.sqrt(arr)
    else:
        return umath.sqrt(arr)

def floor(arr):
    r"""
    The largest integer
    
    .. raw:: latex
    
       $\le x$, $\lfloor x \rfloor$
       
    ..

        >>> print floor(2.3)
        2.0
        >>> print floor(array((-1.5,2,2.5)))
        [-2., 2., 2.,]
        >>> from fipy.variables.variable import Variable
        >>> floor(Variable(value = (-1.5,2,2.5), unit = "m**2"))
        numerix.floor(Variable(value = PhysicalField([-1.5, 2. , 2.5,],'m**2')))
        >>> print floor(Variable(value = (-1.5,2,2.5), unit = "m**2"))
        [-2., 2., 2.,] m**2

    """
    if _isPhysical(arr):
        return arr.floor()
    elif type(arr) is type(array((0))):
        return NUMERIC.floor(arr)
    else:
        return umath.floor(arr)

def ceil(arr):
    r"""
    The largest integer
    
    .. raw:: latex
    
       $\ge x$, $\lceil x \rceil$
       
    ..

        >>> print ceil(2.3)
        3.0
        >>> print ceil(array((-1.5,2,2.5)))
        [-1., 2., 3.,]
        >>> from fipy.variables.variable import Variable
        >>> ceil(Variable(value = (-1.5,2,2.5), unit = "m**2"))
        numerix.ceil(Variable(value = PhysicalField([-1.5, 2. , 2.5,],'m**2')))
        >>> print ceil(Variable(value = (-1.5,2,2.5), unit = "m**2"))
        [-1., 2., 3.,] m**2

    """
    if _isPhysical(arr):
        return arr.ceil()
    elif type(arr) is type(array((0))):
        return NUMERIC.ceil(arr)
    else:
        return umath.ceil(arr)

def exp(arr):
    r"""
    Natural exponent of
    
    .. raw:: latex
    
       $x$, $e^x$
       
    ..

    """
    if _isPhysical(arr):
        return arr.exp()
    elif type(arr) is type(array((0))):
        return NUMERIC.exp(arr)
    else:
        return umath.exp(arr)
        

def log(arr):
    r"""
    Natural logarithm of
    
    .. raw:: latex
    
       $x$, $\ln x \equiv \log_e x$
       
    ..

        >>> print log(10)
        2.30258509299
        >>> print log(array((0.1,1,10)))
        [-2.30258509, 0.        , 2.30258509,]
        >>> from fipy.variables.variable import Variable
        >>> log(Variable(value = (0.1,1,10)))
        numerix.log(Variable(value = [  0.1,  1. , 10. ,]))
        >>> print log(Variable(value = (0.1,1,10)))
        [-2.30258509, 0.        , 2.30258509,]
    """
    if _isPhysical(arr):
        return arr.log()
    elif type(arr) is type(array((0))):
        return NUMERIC.log(arr)
    else:
        return umath.log(arr)

def conjugate(arr):
    r"""
    Complex conjugate of
    
    .. raw:: latex
    
       $z = x + i y$, $z^\star = x - i y$
       
    ..

        >>> print conjugate(3 + 4j)
        (3-4j)
        >>> print conjugate(array((3 + 4j, -2j, 10)))
        [  3.-4.j,  0.+2.j, 10.-0.j,]
        >>> from fipy.variables.variable import Variable
        >>> conjugate(Variable(value = (3 + 4j, -2j, 10), unit = "ohm"))
        numerix.conjugate(Variable(value = PhysicalField([  3.+4.j,  0.-2.j, 10.+0.j,],'ohm')))
        >>> print conjugate(Variable(value = (3 + 4j, -2j, 10), unit = "ohm"))
        [  3.-4.j,  0.+2.j, 10.-0.j,] ohm
    """
    if _isPhysical(arr):
        return arr.conjugate()
    elif type(arr) is type(array((0))):
        return NUMERIC.conjugate(arr)
    else:
        return umath.conjugate(arr)

        # conjugate
        
#########################
#                       #
#   Vector operations   #
#                       #
#########################

def crossProd(v1,v2):
    r"""
    Vector cross-product of
    
    .. raw:: latex
    
       $\vec{v}_1$ and $\vec{v}_2$, $\vec{v}_1 \times \vec{v}_2$
       
    ..

    """
    v1n = NUMERIC.reshape(v1, (-1, 3))
    v2n = NUMERIC.reshape(v2, (-1, 3))

    out = NUMERIC.transpose(array((v1n[:,1] * v2n[:,2] - v1n[:,2] * v2n[:,1],
			    v1n[:,2] * v2n[:,0] - v1n[:,0] * v2n[:,2],
			    v1n[:,0] * v2n[:,1] - v1n[:,1] * v2n[:,0])))
    return NUMERIC.reshape(out, NUMERIC.shape(v1))

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
    ni, nj = NUMERIC.shape(a1)
    result = NUMERIC.zeros((ni,),'d')
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

# Necessary because LLNL hires stupidheads
def MAtake(array, indices, fill = 0, axis = 0):
    """
    Replaces `MA.take`. `MA.take` does not always work when
    `indices` is a masked array.
    """
    tmp = MA.take(array, MA.filled(indices, fill), axis = axis)
    if indices.mask() is not None and tmp.shape != indices.mask().shape:
        mask = MA.repeat(indices.mask()[...,NUMERIC.NewAxis],tmp.shape[-1],len(tmp.shape)-1)
        if tmp.mask() is not None:
            mask = NUMERIC.logical_or(tmp.mask(), mask)
    else:
        mask = indices.mask()
    return MA.array(data = tmp, mask = mask)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__":
    _test() 
