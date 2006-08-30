#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "numerix.py"
 #                                    created: 1/10/04 {10:23:17 AM} 
 #                                last update: 7/21/06 {2:36:08 PM} 
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
   >>> print float(v)
   0.0

Take the tangent of a int.

   >>> tan(0)
   0.0
   
Take the tangent of an array.

   >>> tan(array((0,0,0)))
   [ 0., 0., 0.,]
   
Eventually, this module will be the only place in the code where `Numeric` (or
`numarray` (or `scipy_core`)) is explicitly imported.
"""

__docformat__ = 'restructuredtext'

import Numeric as NUMERIC
from Numeric import *

## import numarray as NUMERIC
## from numarray import *

## from Numeric import pi, array, zeros, ones, ravel, concatenate, where, arange


import umath
import MA

def _isPhysical(arr):
    """
    Returns `True` if arr is a `Variable` or `PhysicalField`.
    """
    import fipy.variables.variable
    import fipy.tools.dimensions.physicalField

    return isinstance(arr,fipy.variables.variable.Variable) \
	or isinstance(arr,fipy.tools.dimensions.physicalField.PhysicalField)

def take(arr, ids, axis = 0):
    """
    Selects the elements of `arr` corresponding to `ids`.
    """
    
    if _isPhysical(arr):
	return arr.take(ids, axis = axis)    
    elif type(ids) is type(MA.array((0))):
        return MAtake(arr, ids, axis = axis)
    elif type(arr) is type(array((0))):
	return NUMERIC.take(arr, ids, axis = axis)
    elif type(arr) is type(MA.array((0))):
	return MA.take(arr, ids, axis = axis)
    else:
	raise TypeError, 'cannot take from object ' + str(arr)
    
def put(arr, ids, values):
    """
    The opposite of `take`.  The values of `arr` at the locations specified by
    `ids` are set to the corresponding value of `values`.
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
    Change the shape of `arr` to `shape`, as long as the product of all the
    lenghts of all the axes is constant (the total number of elements does not
    change).
    """
    if _isPhysical(arr):
	return arr.reshape(shape)
    elif type(arr) is type(array((0))):
	return NUMERIC.reshape(arr, tuple(shape))
    elif type(arr) is type(MA.array((0))):
        return MA.reshape(arr, shape)
    else:
	raise TypeError, 'cannot reshape object ' + str(arr)

def getShape(arr):
    """
    Return the shape of `arr`
    
        >>> getShape(1)
        ()
        >>> getShape(1.)
        ()
        >>> from fipy.variables.variable import Variable
        >>> getShape(Variable(1))
        ()
        >>> getShape(Variable(1.))
        ()
        >>> getShape(Variable(1., unit = "m"))
        ()
        >>> getShape(Variable("1 m"))
        ()
    """
    if _isPhysical(arr):
        return arr.getShape()
    elif type(arr) in (type(array(0)), type(MA.array(0))):
        return arr.shape
    elif type(arr) in (type(()), type([])):
        return (len(arr),)
    elif type(arr) in (type(1), type(1.)):
        return ()
    else:
        return array(arr).shape

def sum(arr, index = 0):
    """
    The sum of all the elements of `arr` along the specified axis.
    """
    if _isPhysical(arr):
	return arr.sum(index)
    elif type(arr) is type(array((0))):
	return NUMERIC.sum(arr, index)
    elif type(arr) is type(MA.array((0))):
	return MA.sum(arr, index)
    else:        
	raise TypeError, 'cannot sum object ' + str(arr)

def tostring(arr, max_line_width = None, precision = None, suppress_small = None, separator = ' ', array_output = 0):
    r"""
    Returns a textual representation of a number or field of numbers.  Each
    dimension is indicated by a pair of matching square brackets (`[]`), within
    which each subset of the field is output.  The orientation of the dimensions
    is as follows: the last (rightmost) dimension is always horizontal, so that
    the frequent rank-1 fields use a minimum of screen real-estate.  The
    next-to-last dimesnion is displayed vertically if present and any earlier
    dimension is displayed with additional bracket divisions.
    
    :Parameters:
        - `max\_line\_width`: the maximum number of characters used in a single
          line.  Default is `sys.output_line_width` or 77.
        - `precision`: the number of digits after the decimal point. Default is 
          `sys.float_output_precision` or 8.
        - `suppress_small`: whether small values should be suppressed (and 
          output as `0`). Default is `sys.float_output_suppress_small` or `false`.
        - `separator`: what character string to place between two numbers.
        - `array_output`: Format output for an `eval`. Only used if `arr` is a 
          `Numeric` `array`.
    """
    if _isPhysical(arr):
        return arr.tostring(max_line_width = max_line_width, 
                            precision = precision, 
                            suppress_small = suppress_small, 
                            separator = separator)
    elif type(arr) in [type(array((0))), type(MA.array((0)))]:
        return NUMERIC.array2string(arr, max_line_width = max_line_width, 
                                    precision = precision, 
                                    suppress_small = suppress_small, 
                                    separator = separator, 
                                    array_output = array_output)
    elif type(arr) in [type(0), type(0.)]:
        from ArrayPrinter import _floatFormat, _formatFloat
        format, length = _floatFormat(array([arr]), precision = precision, suppress_small = suppress_small)
        return _formatFloat(arr, format = format)
    else:        
        raise TypeError, 'cannot convert ' + str(arr) + ' to string'
        
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
    
        >>> print tostring(arccos(0.0), precision = 3)
         1.571
         
    If SciPy has been loaded, the next test will return `NaN`, otherwise it will
    generate `OverflowError: math range error`
    
        >>> try: 
        ...     print str(arccos(2.0)) == "nan"
        ... except (OverflowError, ValueError):
        ...     print 1
        1

        >>> print tostring(arccos(array((0,0.5,1.0))), precision = 3)
        [ 1.571  1.047  0.   ]
        >>> from fipy.variables.variable import Variable
        >>> arccos(Variable(value = (0,0.5,1.0)))
        numerix.arccos(Variable(value = [ 0. , 0.5, 1. ,]))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    ..
       
        >>> print tostring(arccos(Variable(value = (0,0.5,1.0))), precision = 3)
        [ 1.571  1.047  0.   ]
        
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

    If SciPy has been loaded, the next test will return `NaN`, otherwise it will
    generate `OverflowError: math range error`
    
        >>> try: 
        ...     print str(arccosh(0.0)) == "nan"
        ... except (OverflowError, ValueError):
        ...     print 1
        1

        >>> print tostring(arccosh(array((1,2,3))), precision = 3)
        [ 0.     1.317  1.763]
        >>> from fipy.variables.variable import Variable
        >>> arccosh(Variable(value = (1,2,3)))
        numerix.arccosh(Variable(value = [ 1., 2., 3.,]))
        >>> print tostring(arccosh(Variable(value = (1,2,3))), precision = 3)
        [ 0.     1.317  1.763]
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
    
        >>> print tostring(arcsin(1.0), precision = 3)
         1.571
         
    If SciPy has been loaded, the next test will return `NaN`, otherwise it will
    generate `OverflowError: math range error`
    
        >>> try: 
        ...     print str(arcsin(2.0)) == "nan"
        ... except (OverflowError, ValueError):
        ...     print 1
        1

        >>> print tostring(arcsin(array((0,0.5,1.0))), precision = 3)
        [ 0.     0.524  1.571]
        >>> from fipy.variables.variable import Variable
        >>> arcsin(Variable(value = (0,0.5,1.0)))
        numerix.arcsin(Variable(value = [ 0. , 0.5, 1. ,]))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    ..
        
        >>> print tostring(arcsin(Variable(value = (0,0.5,1.0))), precision = 3)
        [ 0.     0.524  1.571]
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

        >>> print tostring(arcsinh(1.0), precision = 3)
         0.881
        >>> print tostring(arcsinh(array((1,2,3))), precision = 3)
        [ 0.881  1.444  1.818]
        >>> from fipy.variables.variable import Variable
        >>> arcsinh(Variable(value = (1,2,3)))
        numerix.arcsinh(Variable(value = [ 1., 2., 3.,]))
        >>> print tostring(arcsinh(Variable(value = (1,2,3))), precision = 3)
        [ 0.881  1.444  1.818]
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
    
        >>> print tostring(arctan(1.0), precision = 3)
         0.785
        >>> print tostring(arctan(array((0,0.5,1.0))), precision = 3)
        [ 0.     0.464  0.785]
        >>> from fipy.variables.variable import Variable
        >>> arctan(Variable(value = (0,0.5,1.0)))
        numerix.arctan(Variable(value = [ 0. , 0.5, 1. ,]))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    ..
    
        >>> print tostring(arctan(Variable(value = (0,0.5,1.0))), precision = 3)
        [ 0.     0.464  0.785]
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

        >>> print tostring(arctan2(3.0, 3.0), precision = 3)
         0.785
        >>> print tostring(arctan2(array((0, 1, 2)), 2), precision = 3)
        [ 0.     0.464  0.785]
        >>> from fipy.variables.variable import Variable
        >>> arctan2(Variable(value = (0, 1, 2)), 2)
        (numerix.arctan2(Variable(value = [ 0., 1., 2.,]), 2))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    ..

        >>> print tostring(arctan2(Variable(value = (0, 1, 2)), 2), precision = 3)
        [ 0.     0.464  0.785]
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
    
        >>> print tostring(arctanh(0.5), precision = 3)
         0.549
        >>> print tostring(arctanh(array((0,0.25,0.5))), precision = 3)
        [ 0.     0.255  0.549]
        >>> from fipy.variables.variable import Variable
        >>> arctanh(Variable(value = (0,0.25,0.5)))
        numerix.arctanh(Variable(value = [ 0.  , 0.25, 0.5 ,]))
        >>> print tostring(arctanh(Variable(value = (0,0.25,0.5))), precision = 3)
        [ 0.     0.255  0.549]
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

        >>> print tostring(cos(2*pi/6), precision = 3)
         0.5
        >>> print tostring(cos(array((0,2*pi/6,pi/2))), precision = 3, suppress_small = 1)
        [ 1.   0.5  0. ]
        >>> from fipy.variables.variable import Variable
        >>> cos(Variable(value = (0,2*pi/6,pi/2), unit = "rad"))
        numerix.cos(Variable(value = PhysicalField([ 0.        , 1.04719755, 1.57079633,],'rad')))
        >>> print tostring(cos(Variable(value = (0,2*pi/6,pi/2), unit = "rad")), suppress_small = 1)
        [ 1.   0.5  0. ]
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
        >>> print tostring(cosh(array((0,1,2))), precision = 3)
        [ 1.     1.543  3.762]
        >>> from fipy.variables.variable import Variable
        >>> cosh(Variable(value = (0,1,2)))
        numerix.cosh(Variable(value = [ 0., 1., 2.,]))
        >>> print tostring(cosh(Variable(value = (0,1,2))), precision = 3)
        [ 1.     1.543  3.762]
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

        >>> print tostring(tan(pi/3), precision = 3)
         1.732
        >>> print tostring(tan(array((0,pi/3,2*pi/3))), precision = 3)
        [ 0.     1.732 -1.732]
        >>> from fipy.variables.variable import Variable
        >>> tan(Variable(value = (0,pi/3,2*pi/3), unit = "rad"))
        numerix.tan(Variable(value = PhysicalField([ 0.        , 1.04719755, 2.0943951 ,],'rad')))
        >>> print tostring(tan(Variable(value = (0,pi/3,2*pi/3), unit = "rad")), precision = 3)
        [ 0.     1.732 -1.732]
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

        >>> print tostring(tanh(1), precision = 3)
         0.762
        >>> print tostring(tanh(array((0,1,2))), precision = 3)
        [ 0.     0.762  0.964]
        >>> from fipy.variables.variable import Variable
        >>> tanh(Variable(value = (0,1,2)))
        numerix.tanh(Variable(value = [ 0., 1., 2.,]))
        >>> print tostring(tanh(Variable(value = (0,1,2))), precision = 3)
        [ 0.     0.762  0.964]
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
        >>> print tostring(sinh(array((0,1,2))), precision = 3)
        [ 0.     1.175  3.627]
        >>> from fipy.variables.variable import Variable
        >>> sinh(Variable(value = (0,1,2)))
        numerix.sinh(Variable(value = [ 0., 1., 2.,]))
        >>> print tostring(sinh(Variable(value = (0,1,2))), precision = 3)
        [ 0.     1.175  3.627]
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

        >>> print tostring(sqrt(2), precision = 3)
         1.414
        >>> print tostring(sqrt(array((1,2,3))), precision = 3)
        [ 1.     1.414  1.732]
        >>> from fipy.variables.variable import Variable
        >>> sqrt(Variable(value = (1, 2, 3), unit = "m**2"))
        numerix.sqrt(Variable(value = PhysicalField([ 1., 2., 3.,],'m**2')))
        >>> print tostring(sqrt(Variable(value = (1, 2, 3), unit = "m**2")), precision = 3)
        [ 1.     1.414  1.732] m

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

        >>> print tostring(log(10), precision = 3)
         2.303
        >>> print tostring(log(array((0.1,1,10))), precision = 3)
        [-2.303  0.     2.303]
        >>> from fipy.variables.variable import Variable
        >>> log(Variable(value = (0.1,1,10)))
        numerix.log(Variable(value = [  0.1,  1. , 10. ,]))
        >>> print tostring(log(Variable(value = (0.1,1,10))), precision = 3)
        [-2.303  0.     2.303]
    """
    if _isPhysical(arr):
        return arr.log()
    elif type(arr) is type(array((0))):
        return NUMERIC.log(arr)
    else:
        return umath.log(arr)

pythonmax = max
def max(arr):
    r"""
    max function

    >>> from fipy.tools.dimensions.physicalField import PhysicalField
    >>> print max(PhysicalField(value = (0.1, -0.2, 0.3), unit = 'm'))
    0.3 m
    >>> print max(array((0.1, -0.2, 0.3)))
    0.3
    >>> from fipy.variables.variable import Variable
    >>> print max(Variable(value = (0.1, -0.2, 0.3)))
    0.3
    
    """

    if type(arr) in (type(array(0)), type(MA.array(0)), type(()), type([])) or \
       (_isPhysical(arr) and (not arr.getUnit() is '1') and not arr.getUnit().isDimensionless()):
        return pythonmax(arr)
    else:
        return pythonmax(array(arr))

pythonmin = min
def min(arr):
    r"""
    min function

    >>> from fipy.tools.dimensions.physicalField import PhysicalField
    >>> print min(PhysicalField(value = (0.1, -0.2, 0.3), unit = 'm'))
    -0.2 m
    >>> print min(array((0.1, -0.2, 0.3)))
    -0.2
    >>> from fipy.variables.variable import Variable
    >>> print min(Variable((0.1, -0.2, 0.3)))
    -0.2
    
    """

    if type(arr) in (type(array(0)), type(MA.array(0)), type(()), type([])) or \
       (_isPhysical(arr) and (not arr.getUnit() is '1') and not arr.getUnit().isDimensionless()):
        return pythonmin(arr)
    else:
        return pythonmin(array(arr))

def conjugate(arr):
    r"""
    Complex conjugate of
    
    .. raw:: latex
    
       $z = x + i y$, $z^\star = x - i y$
       
    ..

        >>> print conjugate(3 + 4j)
        (3-4j)
        >>> print allclose(conjugate(array((3 + 4j, -2j, 10))), (3 - 4j, 2j, 10))
        1
        >>> from fipy.variables.variable import Variable
        >>> var = conjugate(Variable(value = (3 + 4j, -2j, 10), unit = "ohm"))
        >>> print var.getUnit()
        <PhysicalUnit ohm>
        >>> print allclose(var.getNumericValue(), (3 - 4j, 2j, 10))
        1
        
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

    Test that Variables are returned as Variables.

       >>> from fipy.meshes.grid2D import Grid2D
       >>> mesh = Grid2D(nx=2, ny=1)
       >>> from fipy.variables.vectorCellVariable import VectorCellVariable
       >>> v1 = VectorCellVariable(mesh=mesh, value=((0,1),(1,2)))
       >>> v2 = array(((0,1),(1,2)))
       >>> dot(v1, v2)._getVariableClass()
       <class 'fipy.variables.cellVariable.CellVariable'>
       >>> dot(v2, v1)._getVariableClass()
       <class 'fipy.variables.cellVariable.CellVariable'>
       >>> print dot(v1, v2)
       [ 1., 5.,]
       >>> dot(v1, v1)._getVariableClass()
       <class 'fipy.variables.cellVariable.CellVariable'>
       >>> print dot(v1, v1)
       [ 1., 5.,]
       >>> type(dot(v2, v2))
       <type 'array'>
       >>> print dot(v2, v2)
       [1,5,]
       
    
    """

    ## have to check MA since MA's have dot() method!!!
    if hasattr(a1, 'dot') and not (type(a1) is type(MA.array(0))):
        return a1.dot(a2)
    elif hasattr(a2, 'dot') and not (type(a2) is type(MA.array(0))):
        return a2.dot(a1)
    else:
        return sum((a1*a2)[:], axis)

def sqrtDot(a1, a2):
    """Return array of square roots of vector dot-products
    for arrays a1 and a2 of vectors v1 and v2
    
    Usually used with v1==v2 to return magnitude of v1.
    """
    from fipy.tools.inline import inline

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
    from fipy.tools.inline import inline
    
    unit1 = unit2 = 1
    if _isPhysical(a1):
	unit1 = a1.inBaseUnits().getUnit()
        a1 = a1.getNumericValue()
    if _isPhysical(a2):
	unit2 = a2.inBaseUnits().getUnit()
        a2 = a2.getNumericValue()
    ni, nj = NUMERIC.shape(a1)
    result1 = NUMERIC.zeros((ni,),'d')

    inline._runInline("""
        int i;
        for (i = 0; i < ni; i++)
	{
	    int j;
            result1(i) = 0.;
            for (j = 0; j < nj; j++)
            {
	        result1(i) += a1(i,j) * a2(i,j);
            }
            result1(i) = sqrt(result1(i));
        }
    """,result1 = result1, a1 = a1, a2 = a2, ni = ni, nj = nj)


    ##result = inline._runInline("""
##        int j;
##        ((double *) result->data)[i] = 0.;
##        for (j = 0; j < NJ; j++)
##        {
##            ((double *) result->data)[i] += a1(i,j) * a2(i,j);
##        }
##        ((double *) result->data)[i] = sqrt(((double *) result->data)[i]);        
##    """, a1 = a1, a2 = a2, ni = ni, NJ = nj)

    
    if unit1 != 1 or unit2 != 1:
	from fipy.tools.dimensions.physicalField import PhysicalField
	result1 = PhysicalField(value = result, unit = (unit1 * unit2)**0.5)
    return result1

def allequal(first, second):
    """
    Returns `true` if every element of `first` is equal to the corresponding
    element of `second`.
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
    r"""
    Tests whether or not `first` and `second` are equal, subect to the given
    relative and absolute tolerances, such that::
        
        | first - second | < atol + rtol * | second |
        
    This means essentially that both elements are small compared to `atol` or
    their difference divided by `second`'s value is small compared to `rtol`.
    """
    if _isPhysical(first):
	return first.allclose(other = second, atol = atol, rtol = rtol)
    elif _isPhysical(second):
	return second.allclose(other = first, atol = atol, rtol = rtol)
    else:
	return MA.allclose(first, second, atol = atol, rtol = rtol)

def MAtake(array, indices, fill = 0, axis = 0):
    """
    Replaces `MA.take`. `MA.take` does not always work when
    `indices` is a masked array.

       
    """

    tmp = MA.take(array, MA.filled(indices, fill), axis = axis)

    if hasattr(indices, 'mask'):
        if indices.mask() is not None and tmp.shape != indices.mask().shape:
            mask = MA.repeat(indices.mask()[...,NUMERIC.NewAxis],tmp.shape[-1],len(tmp.shape)-1)
            if tmp.mask() is not None:
                mask = NUMERIC.logical_or(tmp.mask(), mask)
        else:
            mask = indices.mask()
    else:
        mask = None
        
    return MA.array(data = tmp, mask = mask)

def indices(dimensions, typecode=None):
    """indices(dimensions,typecode=None) returns an array representing a grid
    of indices with row-only, and column-only variation.

       >>> NUMERIC.allclose(NUMERIC.array(indices((4, 6))), NUMERIC.indices((4,6)))
       1
       >>> NUMERIC.allclose(NUMERIC.array(indices((4, 6, 2))), NUMERIC.indices((4, 6, 2)))
       1
       >>> NUMERIC.allclose(NUMERIC.array(indices((1,))), NUMERIC.indices((1,)))
       1
       >>> NUMERIC.allclose(NUMERIC.array(indices((5,))), NUMERIC.indices((5,)))
       1
  
    """

    lst = []

    if len(dimensions) == 1:
        lst.append(NUMERIC.arange(dimensions[0]))
    elif len(dimensions) == 2:
        ## copy() methods are used to force contiguous arrays
        lst = [NUMERIC.swapaxes(NUMERIC.resize(NUMERIC.arange(dimensions[0]), (dimensions[1], dimensions[0])), 0, 1).copy(),
               NUMERIC.resize(NUMERIC.arange(dimensions[1]), dimensions).copy()]
    else:
        tmp = NUMERIC.ones(dimensions, typecode)
        lst = []
        for i in range(len(dimensions)):
            lst.append(NUMERIC.add.accumulate(tmp, i,) - 1)

    ## we don't turn the list back into an array because that is expensive and not required
    return lst

def getTypecode(arr):
    """
    
    Returns the `typecode()` of the array or `Variable`. Also returns a meaningful
    typecode for ints and floats.

        >>> getTypecode(1)
        'l'
        >>> getTypecode(1.)
        'd'
        >>> getTypecode(array(1))
        'l'
        >>> getTypecode(array(1.))
        'd'
        >>> from fipy.variables.variable import Variable
        >>> getTypecode(Variable(1.))
        'd'
        >>> getTypecode(Variable(1))
        'l'
        >>> getTypecode([0])
        Traceback (most recent call last):
              ...
        TypeError: No typecode for object

    """
    
    if hasattr(arr, 'getTypecode'):
        return arr.getTypecode()
    elif type(arr) is type(array(0)):
        return arr.typecode()
    elif type(arr) is type(0):
        return 'l'
    elif type(arr) is type(0.):
        return 'd'
    else:
        raise TypeError, "No typecode for object"

    
if not hasattr(NUMERIC, 'empty'):
    def empty(shape, dtype='d', order='C'):
        """
        `ones()` and `zeros()` are really slow ways to create arrays. NumPy
        provides a routine:
            
            empty((d1,...,dn),dtype=float,order='C') will return a new array of
            shape (d1,...,dn) and given type with all its entries
            uninitialized. This can be faster than zeros.
            
        We approximate this routine when unavailable, but note that `order` is
        ignored when using Numeric.
        """
        from fipy.tools.inline import inline

        return inline._optionalInline(_emptyIn, _emptyPy, shape, dtype)
    
    def _emptyPy(shape, dtype):
        return NUMERIC.zeros(shape, dtype)

    def _emptyIn(shape, dtype):
        from scipy import weave
        
        local_dict = {'shape': shape, 'dtype': dtype}
        
        code = """
PyObject *op;
PyArrayObject *ret;

char type_char='l';
char *type = &type_char;
int nd, dimensions[MAX_DIMS];

if ((nd=PySequence_Length(shape)) == -1) {
    PyErr_Clear();
    if (!(op = PyNumber_Int(shape))) return NULL;
    nd = 1;
    dimensions[0] = PyInt_AsLong(op);
    Py_DECREF(op);
} else {
    if (nd > MAX_DIMS) {
        fprintf(stderr, "Maximum number of dimensions = %d\\n", MAX_DIMS);
        PyErr_SetString(PyExc_ValueError, "Number of dimensions is too large");
        return NULL;
    }
    for(int i=0; i<nd; i++) {
        if( (op=PySequence_GetItem(shape,i))) {
            dimensions[i]=PyInt_AsLong(op);
            Py_DECREF(op);
        }
        if(PyErr_Occurred()) return NULL;
    }
}
if ((ret = (PyArrayObject *)PyArray_FromDims(nd, dimensions, dtype[0])) == NULL) {
    return NULL;
}

return_val = PyArray_Return(ret);

// refcounting bug in weave. See: "weave: Note on ref counts" in
// weave/scxx/notes.txt
while (return_val.refcount() > 1) {
    Py_DECREF((PyObject *) return_val);
}
"""

        return weave.inline(code,
                     local_dict.keys(),
                     local_dict=local_dict,
                     type_converters=weave.converters.blitz,
                     compiler = 'gcc',
                     verbose = 0,
                     support_code = """
#define MAX_DIMS 30
                     """,
                     extra_compile_args =['-O3'])

    
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__":
    _test() 
