#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "numerix.py"
 #
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
   >>> var = Variable(value=0)

Take the tangent of such a variable. The returned value is itself a
`Variable`.

   >>> v = tan(var)
   >>> v
   numerix.tan(Variable(value=array(0)))
   >>> print float(v)
   0.0

Take the tangent of a int.

   >>> tan(0)
   0.0
   
Take the tangent of an array.

   >>> print tan(array((0,0,0)))
   [ 0.  0.  0.]
   
Eventually, this module will be the only place in the code where `Numeric` (or
`numarray` (or `scipy_core`)) is explicitly imported.
"""

__docformat__ = 'restructuredtext'

import numpy as NUMERIX
from numpy.core import umath
from numpy import newaxis as NewAxis
from numpy import *
from numpy import oldnumeric
try:
    from numpy.core import ma as MA
    numpy_version = 'old'
except ImportError:
    # masked arrays have been moved in numpy 1.1
    from numpy import ma as MA
    numpy_version = 'new'

def zeros(a, t='l'):
    return NUMERIX.zeros(a, t)
def ones(a, t='l'):
    return NUMERIX.ones(a, t)

    

def _isPhysical(arr):
    """
    Returns `True` if arr is a `Variable` or `PhysicalField`.
    """
    from fipy.variables.variable import Variable
    from fipy.tools.dimensions.physicalField import PhysicalField

    return isinstance(arr,Variable) or isinstance(arr,PhysicalField)

def getUnit(arr):
    if hasattr(arr, "getUnit") and callable(arr.getUnit):
        return arr.getUnit()
    else:
        from fipy.tools.dimensions import physicalField
        return physicalField._unity
        
def put(arr, ids, values):
    """
    The opposite of `take`.  The values of `arr` at the locations
    specified by `ids` are set to the corresponding value of `values`.

    The following is to test improvments to puts with masked arrays.
    Places in the code were assuming incorrect put behavior.

       >>> maskValue = 999999

       >>> arr = zeros(3, 'l')
       >>> ids = MA.masked_values((2, maskValue), maskValue)
       >>> values = MA.masked_values((4, maskValue), maskValue)
       >>> put(arr, ids, values) ## this should work 
       >>> print arr
       [0 0 4]

       >>> arr = MA.masked_values((maskValue, 5, 10), maskValue)
       >>> ids = MA.masked_values((2, maskValue), maskValue)
       >>> values = MA.masked_values((4, maskValue), maskValue)
       >>> put(arr, ids, values) 
       >>> print arr ## works as expected
       [-- 5 4]
       
       >>> arr = MA.masked_values((maskValue, 5, 10), maskValue)
       >>> ids = MA.masked_values((maskValue, 2), maskValue)
       >>> values = MA.masked_values((4, maskValue), maskValue)
       >>> put(arr, ids, values)
       >>> print arr ## should be [-- 5 --] maybe??
       [-- 5 999999]
    
    """

    if _isPhysical(arr):
        arr.put(ids, values)
    elif MA.isMaskedArray(arr):
        if NUMERIX.sometrue(MA.getmaskarray(ids)):
            if numpy_version == 'old':
                pvalues = MA.array(values, mask=MA.getmaskarray(ids))
            else:
                pvalues = MA.array(values.filled(), mask=MA.getmaskarray(ids))
            MA.put(arr, ids.compressed(), pvalues.compressed())
        else:
            MA.put(arr, ids, values)
    elif MA.isMaskedArray(ids):
        NUMERIX.put(arr, ids.compressed(), MA.array(values, mask=MA.getmaskarray(ids)).compressed())
    else:
        NUMERIX.put(arr, ids, values)
        
def reshape(arr, shape):
    """
    Change the shape of `arr` to `shape`, as long as the product of all the
    lenghts of all the axes is constant (the total number of elements does not
    change).
    """
    if -1 in shape:
        # e.g., NUMERIX.reshape(a, (-1, N)) fails if N == 0
        oldShape = array(getShape(arr))
        oldShape[oldShape == 0] = 1

        if hasattr(shape, 'index'):
            index = shape.index(-1)
        else:
            index = list(shape).index(-1)
            
        left = shape[:index]
        right = shape[index+1:]
        newShape = array(left + right)
        newShape[newShape == 0] = 1
        
        shape = left + (oldShape.prod() / newShape.prod(),) + right
        
    if _isPhysical(arr):
        return arr.reshape(shape)
    elif type(arr) is type(array((0))):
        return NUMERIX.reshape(arr, tuple(shape))
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
        >>> getShape(Variable(1., unit="m"))
        ()
        >>> getShape(Variable("1 m"))
        ()
    """
    if hasattr(arr, "shape"):
        return arr.shape
    elif type(arr) in (type(()), type([])):
        return (len(arr),)
    elif type(arr) in (type(1), type(1.)):
        return ()
    else:
        raise AttributeError, "No attribute 'shape'"

def rank(a):
    """
    Get the rank of sequence a (the number of dimensions, not a matrix rank)
    The rank of a scalar is zero.
    
    .. note::
        
       The rank of a `MeshVariable` is for any single element. E.g., A
       `CellVariable` containing scalars at each cell, and defined on a 9
       element `Grid1D`, has rank 0. If it is defined on a 3x3 `Grid2D`, it is
       still rank 0.
    """
    if hasattr(a, "getRank"):
        return a.getRank()
    else:
        return NUMERIX.rank(a)
        
def sum(arr, axis=0):
    """
    The sum of all the elements of `arr` along the specified axis.
    """
    if _isPhysical(arr):
        return arr.sum(axis)
    elif type(arr) is type(MA.array((0))):
        return MA.sum(arr, axis)
    else:  
        return NUMERIX.sum(arr, axis)
        
def isFloat(arr):
    if isinstance(arr, NUMERIX.ndarray):
        return NUMERIX.issubclass_(arr.dtype.type, float)
    else:
        return NUMERIX.issubclass_(arr.__class__, float)

def isInt(arr):
    if isinstance(arr, NUMERIX.ndarray):
        return NUMERIX.issubclass_(arr.dtype.type, int)
    else:
        return NUMERIX.issubclass_(arr.__class__, int)
    
def tostring(arr, max_line_width=75, precision=8, suppress_small=False, separator=' ', array_output=0):
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


          >>> from fipy import Variable
          >>> print tostring(Variable((1,0,11.2345)), precision=1)
          [  1.    0.   11.2]
          >>> print tostring(array((1,2)), precision=5)
          [1 2]
          >>> print tostring(array((1.12345,2.79)), precision=2)
          [ 1.12  2.79]
          >>> print tostring(1)
          1
          >>> print tostring(array(1))
          1
          >>> print tostring(array([1.23345]), precision=2)
          [ 1.23]
          >>> print tostring(array([1]), precision=2)
          [1]
          >>> print tostring(1.123456, precision=2)
          1.12
          >>> print tostring(array(1.123456), precision=3)
          1.123
          
    
    """

    if _isPhysical(arr):
        return arr.tostring(max_line_width=max_line_width, 
                            precision=precision, 
                            suppress_small=suppress_small, 
                            separator=separator)
    elif isinstance(arr, NUMERIX.ndarray) and arr.shape != ():
        return NUMERIX.array2string(arr,
                                    precision=precision,
                                    max_line_width=max_line_width,
                                    suppress_small=suppress_small,
                                    separator=separator)
    elif isFloat(arr):
        try:
            ## this is for numpy 1.0.4 and above
            ## why has the interface changed again?
            from numpy.core.arrayprint import FloatFormat
            return FloatFormat(NUMERIX.array((arr,)), precision, suppress_small)(arr).strip()       
        except:
            from numpy.core.arrayprint import _floatFormat, _formatFloat
            return _formatFloat(arr, format='%%1.%df' % precision)

    elif isInt(arr):
        from numpy.core.arrayprint import _formatInteger
        return _formatInteger(arr, format='%d')
    else:        
        raise TypeError, 'cannot convert ' + str(arr) + ' to string'
        
################################################
#                                              #
#   mathematical and trigonometric functions   #
#                                              #
################################################

def arccos(arr):
    r"""
    Inverse cosine of :math:`x`, :math:`\cos^{-1} x`
    
    >>> print tostring(arccos(0.0), precision=3)
    1.571
     
    >>> isnan(arccos(2.0))
    True

    >>> print tostring(arccos(array((0,0.5,1.0))), precision=3)
    [ 1.571  1.047  0.   ]
    >>> from fipy.variables.variable import Variable
    >>> arccos(Variable(value=(0,0.5,1.0)))
    numerix.arccos(Variable(value=array([ 0. ,  0.5,  1. ])))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    >>> print tostring(arccos(Variable(value=(0,0.5,1.0))), precision=3)
    [ 1.571  1.047  0.   ]
        
    """
    if _isPhysical(arr):
        return arr.arccos()
    elif type(arr) is type(array((0))):
        return NUMERIX.arccos(arr)
    else:
        return umath.arccos(arr)

def arccosh(arr):
    r"""
    Inverse hyperbolic cosine of :math:`x`, :math:`\cosh^{-1} x`
    
    >>> print arccosh(1.0)
    0.0

    >>> isnan(arccosh(0.0))
    True

    >>> print tostring(arccosh(array((1,2,3))), precision=3)
    [ 0.     1.317  1.763]
    >>> from fipy.variables.variable import Variable
    >>> arccosh(Variable(value=(1,2,3)))
    numerix.arccosh(Variable(value=array([1, 2, 3])))
    >>> print tostring(arccosh(Variable(value=(1,2,3))), precision=3)
    [ 0.     1.317  1.763]
    """
    if _isPhysical(arr):
        return arr.arccosh()
    elif type(arr) is type(array((0))):
        return NUMERIX.arccosh(arr)
    else:
        return umath.arccosh(arr)

def arcsin(arr):
    r"""
    Inverse sine of :math:`x`, :math:`\sin^{-1} x`
    
    >>> print tostring(arcsin(1.0), precision=3)
    1.571
     
    >>> isnan(arcsin(2.0))
    True

    >>> print tostring(arcsin(array((0,0.5,1.0))), precision=3)
    [ 0.     0.524  1.571]
    >>> from fipy.variables.variable import Variable
    >>> arcsin(Variable(value=(0,0.5,1.0)))
    numerix.arcsin(Variable(value=array([ 0. ,  0.5,  1. ])))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    >>> print tostring(arcsin(Variable(value=(0,0.5,1.0))), precision=3)
    [ 0.     0.524  1.571]
    """
    if _isPhysical(arr):
        return arr.arcsin()
    elif type(arr) is type(array((0))):
        return NUMERIX.arcsin(arr)
    else:
        return umath.arcsin(arr)

def arcsinh(arr):
    r"""
    Inverse hyperbolic sine of :math:`x`, :math:`\sinh^{-1} x`

    >>> print tostring(arcsinh(1.0), precision=3)
    0.881
    >>> print tostring(arcsinh(array((1,2,3))), precision=3)
    [ 0.881  1.444  1.818]
    >>> from fipy.variables.variable import Variable
    >>> arcsinh(Variable(value=(1,2,3)))
    numerix.arcsinh(Variable(value=array([1, 2, 3])))
    >>> print tostring(arcsinh(Variable(value=(1,2,3))), precision=3)
    [ 0.881  1.444  1.818]
    """
    if _isPhysical(arr):
        return arr.arcsinh()
    elif type(arr) is type(array((0))):
        return NUMERIX.arcsinh(arr)
    else:
        return umath.arcsinh(arr)

def arctan(arr):
    r"""
    Inverse tangent of :math:`x`, :math:`\tan^{-1} x`

    >>> print tostring(arctan(1.0), precision=3)
    0.785
    >>> print tostring(arctan(array((0,0.5,1.0))), precision=3)
    [ 0.     0.464  0.785]
    >>> from fipy.variables.variable import Variable
    >>> arctan(Variable(value=(0,0.5,1.0)))
    numerix.arctan(Variable(value=array([ 0. ,  0.5,  1. ])))
    
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    >>> print tostring(arctan(Variable(value=(0,0.5,1.0))), precision=3)
    [ 0.     0.464  0.785]
    """
    if _isPhysical(arr):
        return arr.arctan()
    elif type(arr) is type(array((0))):
        return NUMERIX.arctan(arr)
    else:
        return umath.arctan(arr)
                
def arctan2(arr, other):
    r"""
    Inverse tangent of a ratio :math:`x/y`, :math:`\tan^{-1} \frac{x}{y}`

    >>> print tostring(arctan2(3.0, 3.0), precision=3)
    0.785
    >>> print tostring(arctan2(array((0, 1, 2)), 2), precision=3)
    [ 0.     0.464  0.785]
    >>> from fipy.variables.variable import Variable
    >>> arctan2(Variable(value=(0, 1, 2)), 2)
    (numerix.arctan2(Variable(value=array([0, 1, 2])), 2))
        
    .. attention:: 
        
       the next should really return radians, but doesn't
       
    >>> print tostring(arctan2(Variable(value=(0, 1, 2)), 2), precision=3)
    [ 0.     0.464  0.785]
    """
    if _isPhysical(arr):
        return arr.arctan2(other)
    elif _isPhysical(other):
        from fipy.tools.dimensions import physicalField

        return physicalField.PhysicalField(value=arr, unit="rad").arctan2(other)
    elif type(arr) is type(array((0))):
        return NUMERIX.arctan2(arr,other)
    else:
        return umath.arctan2(arr,other)
        
        
def arctanh(arr):
    r"""
    Inverse hyperbolic tangent of :math:`x`, :math:`\tanh^{-1} x`

    >>> print tostring(arctanh(0.5), precision=3)
    0.549
    >>> print tostring(arctanh(array((0,0.25,0.5))), precision=3)
    [ 0.     0.255  0.549]
    >>> from fipy.variables.variable import Variable
    >>> arctanh(Variable(value=(0,0.25,0.5)))
    numerix.arctanh(Variable(value=array([ 0.  ,  0.25,  0.5 ])))
    >>> print tostring(arctanh(Variable(value=(0,0.25,0.5))), precision=3)
    [ 0.     0.255  0.549]
    """
    if _isPhysical(arr):
        return arr.arctanh()
    elif type(arr) is type(array((0))):
        return NUMERIX.arctanh(arr)
    else:
        return umath.arctanh(arr)
        
def cos(arr):
    r"""
    Cosine of :math:`x`, :math:`\cos x`

    >>> print allclose(cos(2*pi/6), 0.5)
    True
    >>> print tostring(cos(array((0,2*pi/6,pi/2))), precision=3, suppress_small=1)
    [ 1.   0.5  0. ]
    >>> from fipy.variables.variable import Variable
    >>> cos(Variable(value=(0,2*pi/6,pi/2), unit="rad"))
    numerix.cos(Variable(value=PhysicalField(array([ 0.        ,  1.04719755,  1.57079633]),'rad')))
    >>> print tostring(cos(Variable(value=(0,2*pi/6,pi/2), unit="rad")), suppress_small=1)
    [ 1.   0.5  0. ]
    """
    if _isPhysical(arr):
        return arr.cos()
    elif type(arr) is type(array((0))):
        return NUMERIX.cos(arr)
    else:
        return umath.cos(arr)

def cosh(arr):
    r"""
    Hyperbolic cosine of :math:`x`, :math:`\cosh x`

    >>> print cosh(0)
    1.0
    >>> print tostring(cosh(array((0,1,2))), precision=3)
    [ 1.     1.543  3.762]
    >>> from fipy.variables.variable import Variable
    >>> cosh(Variable(value=(0,1,2)))
    numerix.cosh(Variable(value=array([0, 1, 2])))
    >>> print tostring(cosh(Variable(value=(0,1,2))), precision=3)
    [ 1.     1.543  3.762]
    """
    if _isPhysical(arr):
        return arr.cosh()
    elif type(arr) is type(array((0))):
        return NUMERIX.cosh(arr)
    else:
        return umath.cosh(arr)

def tan(arr):
    r"""
    Tangent of :math:`x`, :math:`\tan x`

    >>> print tostring(tan(pi/3), precision=3)
    1.732
    >>> print tostring(tan(array((0,pi/3,2*pi/3))), precision=3)
    [ 0.     1.732 -1.732]
    >>> from fipy.variables.variable import Variable
    >>> tan(Variable(value=(0,pi/3,2*pi/3), unit="rad"))
    numerix.tan(Variable(value=PhysicalField(array([ 0.        ,  1.04719755,  2.0943951 ]),'rad')))
    >>> print tostring(tan(Variable(value=(0,pi/3,2*pi/3), unit="rad")), precision=3)
    [ 0.     1.732 -1.732]
    """
    if _isPhysical(arr):
        return arr.tan()
    elif type(arr) is type(array((0))):
        return NUMERIX.tan(arr)
    else:
        return umath.tan(arr)

def tanh(arr):
    r"""
    Hyperbolic tangent of :math:`x`, :math:`\tanh x`

    >>> print tostring(tanh(1), precision=3)
    0.762
    >>> print tostring(tanh(array((0,1,2))), precision=3)
    [ 0.     0.762  0.964]
    >>> from fipy.variables.variable import Variable
    >>> tanh(Variable(value=(0,1,2)))
    numerix.tanh(Variable(value=array([0, 1, 2])))
    >>> print tostring(tanh(Variable(value=(0,1,2))), precision=3)
    [ 0.     0.762  0.964]
    """
    if _isPhysical(arr):
        return arr.tanh()
    elif type(arr) is type(array((0))):
        return NUMERIX.tanh(arr)
    else:
        return umath.tanh(arr)

def log10(arr):
    r"""
    Base-10 logarithm of :math:`x`, :math:`\log_{10} x`

    >>> print log10(10)
    1.0
    >>> print log10(array((0.1,1,10)))
    [-1.  0.  1.]
    >>> from fipy.variables.variable import Variable
    >>> log10(Variable(value=(0.1,1,10)))
    numerix.log10(Variable(value=array([  0.1,   1. ,  10. ])))
    >>> print log10(Variable(value=(0.1,1,10)))
    [-1.  0.  1.]
    """
    if _isPhysical(arr):
        return arr.log10()
    elif type(arr) is type(array((0))):
        return NUMERIX.log10(arr)
    else:
        return umath.log10(arr)

def sin(arr):
    r"""
    Sine of :math:`x`, :math:`\sin x`

    >>> print sin(pi/6)
    0.5
    >>> print sin(array((0,pi/6,pi/2)))
    [ 0.   0.5  1. ]
    >>> from fipy.variables.variable import Variable
    >>> sin(Variable(value=(0,pi/6,pi/2), unit="rad"))
    numerix.sin(Variable(value=PhysicalField(array([ 0.        ,  0.52359878,  1.57079633]),'rad')))
    >>> print sin(Variable(value=(0,pi/6,pi/2), unit="rad"))
    [ 0.   0.5  1. ]
    """
    if _isPhysical(arr):
        return arr.sin()
    elif type(arr) is type(array((0))):
        return NUMERIX.sin(arr)
    else:
        return umath.sin(arr)

def sinh(arr):
    r"""
    Hyperbolic sine of :math:`x`, :math:`\sinh x`

    >>> print sinh(0)
    0.0
    >>> print tostring(sinh(array((0,1,2))), precision=3)
    [ 0.     1.175  3.627]
    >>> from fipy.variables.variable import Variable
    >>> sinh(Variable(value=(0,1,2)))
    numerix.sinh(Variable(value=array([0, 1, 2])))
    >>> print tostring(sinh(Variable(value=(0,1,2))), precision=3)
    [ 0.     1.175  3.627]
    """
    if _isPhysical(arr):
        return arr.sinh()
    elif type(arr) is type(array((0))):
        return NUMERIX.sinh(arr)
    else:
        return umath.sinh(arr)

def sqrt(arr):
    r"""
    Square root of :math:`x`, :math:`\sqrt{x}`

    >>> print tostring(sqrt(2), precision=3)
    1.414
    >>> print tostring(sqrt(array((1,2,3))), precision=3)
    [ 1.     1.414  1.732]
    >>> from fipy.variables.variable import Variable
    >>> sqrt(Variable(value=(1, 2, 3), unit="m**2"))
    numerix.sqrt(Variable(value=PhysicalField(array([1, 2, 3]),'m**2')))
    >>> print tostring(sqrt(Variable(value=(1, 2, 3), unit="m**2")), precision=3)
    [ 1.     1.414  1.732] m
    """
    if _isPhysical(arr):
        return arr.sqrt()
    elif type(arr) is type(array((0))):
        return NUMERIX.sqrt(arr)
    else:
        return umath.sqrt(arr)

def floor(arr):
    r"""
    The largest integer :math:`\le x`, :math:`\lfloor x \rfloor`

    >>> print floor(2.3)
    2.0
    >>> print floor(array((-1.5,2,2.5)))
    [-2.  2.  2.]
    >>> from fipy.variables.variable import Variable
    >>> floor(Variable(value=(-1.5,2,2.5), unit="m**2"))
    numerix.floor(Variable(value=PhysicalField(array([-1.5,  2. ,  2.5]),'m**2')))
    >>> print floor(Variable(value=(-1.5,2,2.5), unit="m**2"))
    [-2.  2.  2.] m**2
    """
    if _isPhysical(arr):
        return arr.floor()
    elif type(arr) is type(array((0))):
        return NUMERIX.floor(arr)
    else:
        return umath.floor(arr)

def ceil(arr):
    r"""
    The largest integer :math:`\ge x`, :math:`\lceil x \rceil`
       
    >>> print ceil(2.3)
    3.0
    >>> print ceil(array((-1.5,2,2.5)))
    [-1.  2.  3.]
    >>> from fipy.variables.variable import Variable
    >>> ceil(Variable(value=(-1.5,2,2.5), unit="m**2"))
    numerix.ceil(Variable(value=PhysicalField(array([-1.5,  2. ,  2.5]),'m**2')))
    >>> print ceil(Variable(value=(-1.5,2,2.5), unit="m**2"))
    [-1.  2.  3.] m**2
    """
    if _isPhysical(arr):
        return arr.ceil()
    elif type(arr) is type(array((0))):
        return NUMERIX.ceil(arr)
    else:
        return umath.ceil(arr)


def sign(arr):
    if _isPhysical(arr):
        return arr.sign()
    elif type(arr) is type(array((0))):
        return NUMERIX.sign(arr)
    else:
        return umath.sign(arr)

def exp(arr):
    r"""
    Natural exponent of :math:`x`, :math:`e^x`
    """
    if _isPhysical(arr):
        return arr.exp()
    elif type(arr) is type(array((0))):
        return NUMERIX.exp(arr)
    else:
        return umath.exp(arr)
        

def log(arr):
    r"""
    Natural logarithm of :math:`x`, :math:`\ln x \equiv \log_e x`

    >>> print tostring(log(10), precision=3)
    2.303
    >>> print tostring(log(array((0.1,1,10))), precision=3)
    [-2.303  0.     2.303]
    >>> from fipy.variables.variable import Variable
    >>> log(Variable(value=(0.1,1,10)))
    numerix.log(Variable(value=array([  0.1,   1. ,  10. ])))
    >>> print tostring(log(Variable(value=(0.1,1,10))), precision=3)
    [-2.303  0.     2.303]
    """
    if _isPhysical(arr):
        return arr.log()
    elif type(arr) is type(array((0))):
        return NUMERIX.log(arr)
    else:
        return umath.log(arr)

def conjugate(arr):
    r"""
    Complex conjugate of :math:`z = x + i y`, :math:`z^\star = x - i y`
       
    >>> print conjugate(3 + 4j) == 3 - 4j
    True
    >>> print allclose(conjugate(array((3 + 4j, -2j, 10))), (3 - 4j, 2j, 10))
    1
    >>> from fipy.variables.variable import Variable
    >>> var = conjugate(Variable(value=(3 + 4j, -2j, 10), unit="ohm"))
    >>> print var.getUnit()
    <PhysicalUnit ohm>
    >>> print allclose(var.getNumericValue(), (3 - 4j, 2j, 10))
    1
    """
    if _isPhysical(arr):
        return arr.conjugate()
    elif type(arr) is type(array((0))):
        return NUMERIX.conjugate(arr)
    else:
        return umath.conjugate(arr)

        # conjugate
        
#########################
#                       #
#   Vector operations   #
#                       #
#########################

def dot(a1, a2, axis=0):
    """
    return array of vector dot-products of v1 and v2
    for arrays a1 and a2 of vectors v1 and v2
    
    We can't use :func:`numpy.dot` on an array of vectors

    Test that Variables are returned as Variables.

    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(nx=2, ny=1)
    >>> from fipy.variables.cellVariable import CellVariable
    >>> v1 = CellVariable(mesh=mesh, value=((0,1),(2,3)), rank=1)
    >>> v2 = CellVariable(mesh=mesh, value=((0,1),(2,3)), rank=1)
    >>> dot(v1, v2)._getVariableClass()
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> dot(v2, v1)._getVariableClass()
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> print rank(dot(v2, v1))
    0
    >>> print dot(v1, v2)
    [ 4 10]
    >>> dot(v1, v1)._getVariableClass()
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> print dot(v1, v1)
    [ 4 10]
    >>> v3 = array(((0,1),(2,3)))
    >>> type(dot(v3, v3))
    <type 'numpy.ndarray'>
    >>> print dot(v3, v3)
    [ 4 10]
    """

    ## have to check MA since MA's have dot() method!!!
    ## and numpy arrays now have a dot() method too!!!!

    def isNumpy(aa):
        return type(aa) in (type(MA.array(0)), type(NUMERIX.array(0)))
 
    if hasattr(a1, 'dot') and not isNumpy(a1):
        return a1.dot(a2)
    elif hasattr(a2, 'rdot') and not isNumpy(a2):
        return a2.rdot(a1)
    elif hasattr(a2, 'dot') and not isNumpy(a2):
        # dot() is not commutative with tensors, but if there's no
        # rdot(), what else can we do? Just throw an error?
        return a2.dot(a1)
    else:
        return sum(a1*a2, axis)

def sqrtDot(a1, a2):
    """Return array of square roots of vector dot-products
    for arrays a1 and a2 of vectors v1 and v2
    
    Usually used with v1==v2 to return magnitude of v1.
    """
    from fipy.tools import inline

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
##      int j;
##      result(i) = 0.;
##      for (j = 0; j < nj; j++)
##      {
##          result(i) += a1(i,j) * a2(i,j);
##      }
##      result(i) = sqrt(result(i));
##    """,result = result, a1 = a1, a2 = a2, ni = ni, nj = nj) 
##    return result

def _sqrtDotIn(a1, a2):
    from fipy.tools import inline
    
    unit1 = unit2 = 1
    if _isPhysical(a1):
        unit1 = a1.inBaseUnits().getUnit()
        a1 = a1.getNumericValue()
    if _isPhysical(a2):
        unit2 = a2.inBaseUnits().getUnit()
        a2 = a2.getNumericValue()
    NJ, ni = NUMERIX.shape(a1)
    result1 = NUMERIX.zeros((ni,),'d')

    inline._runInline("""
        int j;
        result1[i] = 0.;
        for (j = 0; j < NJ; j++)
        {
            // result1[i] += a1[i * NJ + j] * a2[i * NJ + j];
            result1[i] += a1[i + j * ni] * a2[i + j * ni];
        }
        result1[i] = sqrt(result1[i]);        
    """,result1=result1, a1=a1, a2=a2, ni=ni, NJ=NJ)


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
        result1 = PhysicalField(value=result, unit=(unit1 * unit2)**0.5)
    return result1

def allequal(first, second):
    """
    Returns `true` if every element of `first` is equal to the corresponding
    element of `second`.
    """
    if _isPhysical(first):
        return first.allequal(second)
    elif _isPhysical(second):
        return second.allequal(first)
    else:
        return MA.allequal(first, second)
            
def allclose(first, second, rtol=1.e-5, atol=1.e-8):
    r"""
    Tests whether or not `first` and `second` are equal, subect to the given
    relative and absolute tolerances, such that::
        
       |first - second| < atol + rtol * |second|
        
    This means essentially that both elements are small compared to ``atol`` or
    their difference divided by ``second``'s value is small compared to ``rtol``.
    """
    if _isPhysical(first):
        return first.allclose(other=second, atol=atol, rtol=rtol)
    elif _isPhysical(second):
        return second.allclose(other=first, atol=atol, rtol=rtol)
    else:
        return MA.allclose(first, second, atol=atol, rtol=rtol)

def isclose(first, second, rtol=1.e-5, atol=1.e-8):
    r"""
    Returns which elements of `first` and `second` are equal, subect to the given
    relative and absolute tolerances, such that::
        
        |first - second| < atol + rtol * |second|

    This means essentially that both elements are small compared to ``atol`` or
    their difference divided by ``second``'s value is small compared to ``rtol``.
    """
    return abs(first - second) < atol + rtol * abs(second)

def take(a, indices, axis=0, fill_value=None):
    """
    Selects the elements of `a` corresponding to `indices`.
    """
           
    if _isPhysical(a):
        taken = a.take(indices, axis=axis)   
    elif type(indices) is type(MA.array((0))):
        ## Replaces `MA.take`. `MA.take` does not always work when
        ## `indices` is a masked array.
        ##
        taken = MA.take(a, MA.filled(indices, 0), axis=axis)
        
        mask = MA.getmask(indices)
        
        if mask is not MA.nomask:
            mask = MA.getmaskarray(indices)
            if taken.shape != mask.shape:
                mask = MA.repeat(mask[NewAxis, ...], taken.shape[0], axis=0)
                mask = MA.mask_or(MA.getmask(taken), mask)

        
        if mask is not MA.nomask:
            taken = MA.array(data=taken, mask=mask)
        else:
            if MA.getmask(taken) is MA.nomask and numpy_version == 'old':
                # numpy 1.1 returns normal array when masked array is filled
                taken = taken.filled()

    elif type(a) in (type(array((0))), type(()), type([])):
        taken = NUMERIX.take(a, indices, axis=axis)
    elif type(a) is type(MA.array((0))):
        taken = MA.take(a, indices, axis=axis)
    else:
        raise TypeError, 'cannot take from %s object: %s' % (type(a), `a`)
               
    if fill_value is not None and type(taken) is type(MA.array((0))):
        taken = taken.filled(fill_value=fill_value)
        
    return taken

def indices(dimensions, typecode=None):
    """indices(dimensions,typecode=None) returns an array representing a grid
    of indices with row-only, and column-only variation.

       >>> NUMERIX.allclose(NUMERIX.array(indices((4, 6))), NUMERIX.indices((4,6)))
       1
       >>> NUMERIX.allclose(NUMERIX.array(indices((4, 6, 2))), NUMERIX.indices((4, 6, 2)))
       1
       >>> NUMERIX.allclose(NUMERIX.array(indices((1,))), NUMERIX.indices((1,)))
       1
       >>> NUMERIX.allclose(NUMERIX.array(indices((5,))), NUMERIX.indices((5,)))
       1
  
    """

    lst = []

    if len(dimensions) == 1:
        lst.append(NUMERIX.arange(dimensions[0]))
    elif len(dimensions) == 2:
        ## copy() methods are used to force contiguous arrays
        lst = [NUMERIX.swapaxes(NUMERIX.resize(NUMERIX.arange(dimensions[0]), (dimensions[1], dimensions[0])), 0, 1).copy(),
               NUMERIX.resize(NUMERIX.arange(dimensions[1]), dimensions).copy()]
    else:
        tmp = NUMERIX.ones(dimensions, typecode)
        lst = []
        for i in range(len(dimensions)):
            lst.append(NUMERIX.add.accumulate(tmp, i,) - 1)

    ## we don't turn the list back into an array because that is expensive and not required
    return lst

def obj2sctype(rep, default=None):
    if _isPhysical(rep):
        sctype = rep.getsctype(default=default)
    else:
        if type(rep) in (type(()), type([])):
            rep = array(rep)
        sctype = NUMERIX.obj2sctype(rep=rep, default=default)
        
    if sctype is None:
        return obj2sctype(type(rep), default=default)
    else:
        return sctype
        

if not hasattr(NUMERIX, 'empty'):
    print 'defining empty'
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
        from fipy.tools import inline

        return inline._optionalInline(_emptyIn, _emptyPy, shape, dtype)
  
    def _emptyPy(shape, dtype):
        return NUMERIX.zeros(shape, dtype)

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
                     compiler='gcc',
                     verbose=0,
                     support_code="""
#define MAX_DIMS 30
                     """,
                     extra_compile_args =['-O3'])

if not (hasattr(NUMERIX, 'savetxt') and hasattr(NUMERIX, 'loadtxt')):
    # if one is present, but not the other, something is wrong and
    # we shouldn't try to patch.
    
    # The following routines were introduced in NumPy 1.0.3
    # c.f. http://projects.scipy.org/scipy/numpy/changeset/3722
    
    # Adapted from matplotlib 
    
    def _getconv(dtype): 
        typ = dtype.type 
        if issubclass(typ, bool_): 
            return lambda x: bool(int(x)) 
        if issubclass(typ, integer): 
            return int 
        elif issubclass(typ, floating): 
            return float 
        elif issubclass(typ, complex): 
            return complex 
        else: 
            return str 
     
     
    def _string_like(obj): 
        try: obj + '' 
        except (TypeError, ValueError): return 0 
        return 1 
     
    def loadtxt(fname, dtype=float, comments='#', delimiter=None, converters=None, 
                skiprows=0, usecols=None, unpack=False): 
        """ 
        Load ASCII data from fname into an array and return the array. 
     
        The data must be regular, same number of values in every row 
     
        fname can be a filename or a file handle.  Support for gzipped files is 
        automatic, if the filename ends in .gz 
     
        See scipy.loadmat to read and write matfiles. 
     
        Example usage: 
     
          X = loadtxt('test.dat')  # data in two columns 
          t = X[:,0] 
          y = X[:,1] 
     
        Alternatively, you can do the same with "unpack"; see below 
     
          X = loadtxt('test.dat')    # a matrix of data 
          x = loadtxt('test.dat')    # a single column of data 
     
     
        dtype - the data-type of the resulting array.  If this is a 
        record data-type, the the resulting array will be 1-d and each row will 
        be interpreted as an element of the array. The number of columns 
        used must match the number of fields in the data-type in this case.  
         
        comments - the character used to indicate the start of a comment 
        in the file 
     
        delimiter is a string-like character used to seperate values in the 
        file. If delimiter is unspecified or none, any whitespace string is 
        a separator. 
     
        converters, if not None, is a dictionary mapping column number to 
        a function that will convert that column to a float.  Eg, if 
        column 0 is a date string: converters={0:datestr2num} 
     
        skiprows is the number of rows from the top to skip 
     
        usecols, if not None, is a sequence of integer column indexes to 
        extract where 0 is the first column, eg usecols=(1,4,5) to extract 
        just the 2nd, 5th and 6th columns 
     
        unpack, if True, will transpose the matrix allowing you to unpack 
        into named arguments on the left hand side 
     
            t,y = load('test.dat', unpack=True) # for  two column data 
            x,y,z = load('somefile.dat', usecols=(3,5,7), unpack=True) 
     
        """ 
     
        if _string_like(fname): 
            if fname.endswith('.gz'): 
                import gzip 
                fh = gzip.open(fname) 
            else: 
                fh = file(fname) 
        elif hasattr(fname, 'seek'): 
            fh = fname 
        else: 
            raise ValueError('fname must be a string or file handle') 
        X = [] 

        dtype = NUMERIX.dtype(dtype) 
        defconv = _getconv(dtype) 
        converterseq = None     
        if converters is None: 
            converters = {} 
            if dtype.names is not None: 
                converterseq = [_getconv(dtype.fields[name][0]) \
                                for name in dtype.names] 
                 
        for i,line in enumerate(fh): 
            if i<skiprows: continue 
            line = line[:line.find(comments)].strip() 
            if not len(line): continue 
            vals = line.split(delimiter) 
            if converterseq is None: 
               converterseq = [converters.get(j,defconv) \
                               for j in xrange(len(vals))] 
            if usecols is not None: 
                row = [converterseq[j](vals[j]) for j in usecols] 
            else: 
                row = [converterseq[j](val) for j,val in enumerate(vals)] 
            if dtype.names is not None: 
                row = tuple(row) 
            X.append(row) 
     
        X = array(X, dtype) 
        r,c = X.shape 
        if r==1 or c==1: 
            X.shape = max([r,c]), 
        if unpack: return X.T 
        else:  return X 
     
     
    # adjust so that fmt can change across columns if desired.  
     
    def savetxt(fname, X, fmt='%.18e',delimiter=' '): 
        """ 
        Save the data in X to file fname using fmt string to convert the 
        data to strings 
     
        fname can be a filename or a file handle.  If the filename ends in .gz, 
        the file is automatically saved in compressed gzip format.  The load() 
        command understands gzipped files transparently. 
     
        Example usage: 
     
        save('test.out', X)         # X is an array 
        save('test1.out', (x,y,z))  # x,y,z equal sized 1D arrays 
        save('test2.out', x)        # x is 1D 
        save('test3.out', x, fmt='%1.4e')  # use exponential notation 
     
        delimiter is used to separate the fields, eg delimiter ',' for 
        comma-separated values 
        """ 
     
        if _string_like(fname): 
            if fname.endswith('.gz'): 
                import gzip 
                fh = gzip.open(fname,'wb') 
            else: 
                fh = file(fname,'w') 
        elif hasattr(fname, 'seek'): 
            fh = fname 
        else: 
            raise ValueError('fname must be a string or file handle') 
     
     
        X = asarray(X) 
        origShape = None 
        if len(X.shape)==1: 
            origShape = X.shape 
            X.shape = len(X), 1 
        for row in X: 
            fh.write(delimiter.join([fmt%val for val in row]) + '\n') 
     
        if origShape is not None: 
            X.shape = origShape 

    
def L1norm(arr):
    r"""
    :Parameters:
      - `arr`: The `array` to evaluate.
      
    :Returns: 
      :math:`\|\mathtt{arr}\|_1 = \sum_{j=1}^{n} |\mathtt{arr}_j|` is the
      :math:`L^1`-norm of :math:`\mathtt{arr}`.
    """
    return add.reduce(abs(arr))
    
def L2norm(arr):
    r"""
    :Parameters:
      - `arr`: The `array` to evaluate.
      
    :Returns: 
      :math:`\|\mathtt{arr}\|_2 = \sqrt{\sum_{j=1}^{n} |\mathtt{arr}_j|^2}` is
      the :math:`L^2`-norm of :math:`\mathtt{arr}`.
    """
    return sqrt(add.reduce(arr**2))
    
def LINFnorm(arr):
    r"""
    :Parameters:
      - `arr`: The `array` to evaluate.
      
    :Returns: 
      :math:`\|\mathtt{arr}\|_\infty = [\sum_{j=1}^{n}
      |\mathtt{arr}_j|^\infty]^\infty = \over{\max}{j} |\mathtt{arr}_j|` is the
      :math:`L^\infty`-norm of :math:`\mathtt{arr}`.
    """
    return max(abs(arr))

def _compressIndexSubspaces(index, i, broadcastshape = ()):
    """
    Starting at element `i` in a selection tuple `index`, squeeze out all
    sequential elements that can be broadcast to the same shape.
    """
    skip = 0
    while i + skip < len(index):
        element = index[i + skip]
        if element is newaxis or isinstance(element, slice):
            skip -= 1
            break
        else:
            element = array(element, intp)

            # numpy accepts tuples of lists of floats, but not arrays of
            # floats. This test is more liberal than that, but has the
            # benefit of not being insane. Unfortunately, NumPy will throw
            # errors later, on evaluation, that should ideally be caught here.
            #    raise IndexError "arrays used as indices must be of integer (or boolean) type"
            
            broadcastshape = _broadcastShape(broadcastshape, element.shape)
            if broadcastshape is None:
                raise ValueError, "shape mismatch: objects cannot be broadcast to a single shape"
        skip += 1

    return broadcastshape, skip

def _indexShape(index, arrayShape):
    """
    Given a selection object `index` and an `arrayShape` for the the object to
    be sliced into, return the shape of the result. Return `None` if the
    indexing is not legal.
    
    "If the length of the selection tuple is larger than N(=`arrayShape.ndim`)
    an error is raised"
    
        >>> _indexShape(index=(1, 2, 3, 4), 
        ...             arrayShape=(10,20,30))
        Traceback (most recent call last):
            ...
        IndexError: invalid index

    "All sequences and scalars in the selection tuple are converted to intp
    indexing arrays."
    
    "All selection tuple objects must be convertible to intp arrays, or slice
    objects, or the Ellipsis (``...``) object"
    
        >>> _indexShape(index=NUMERIX.index_exp[...,2,"text"], 
        ...             arrayShape=(10,20,30,40,50))
        Traceback (most recent call last):
            ...
        ValueError: setting an array element with a sequence.

    .. note::
        
       The following test should throw::
           
           Traceback (most recent call last):
               ...
           IndexError: arrays used as indices must be of integer (or boolean) type

       but it's not straightforward to achieve that. Moreover, NumPy is not even
       consistent, accepting a `tuple` or `list` of `float`, but not a `float`
       `array`. If it's absolutely crucial to obtain that result, see the
       comment in `_compressIndexSubspaces()`.
       
    ..
    
        >>> ind = zeros((2,3,5), float)
        >>> _indexShape(index=NUMERIX.index_exp[...,ind], 
        ...             arrayShape=(10,20,30,40,50))
        (10, 20, 30, 40, 2, 3, 5)

    "Exactly one Ellipsis object will be expanded, any other Ellipsis objects
    will be treated as full slice (':') objects. The Ellipsis object is replaced
    with as many full slice (':') objects as needed to make the length of the
    selection tuple N."

        >>> _indexShape(index=NUMERIX.index_exp[...,2,...,4], 
        ...             arrayShape=(10,20,30,40,50))
        (10, 20, 40)
    
    "If the selection tuple is smaller than N, then as many ':' objects as
    needed are added to the end of the selection tuple so that the modified 
    selection tuple has length N."

        >>> _indexShape(index=NUMERIX.index_exp[:,2], 
        ...             arrayShape=(10,20,30,40,50))
        (10, 30, 40, 50)

    "The shape of all the integer indexing arrays must be broadcastable to the
    same shape"
    
        >>> ind1 = zeros((2,3,4), intp)
        >>> ind2 = zeros((2,3,5), intp)
        >>> _indexShape(index=NUMERIX.index_exp[:,ind1,ind2], 
        ...             arrayShape=(10,20,30,40,50))
        Traceback (most recent call last):
            ...
        ValueError: shape mismatch: objects cannot be broadcast to a single shape

    "In simple cases (i.e. one indexing array and N - 1 slice objects) it does
    exactly what you would expect (concatenation of repeated application of
    basic slicing)."
    
        >>> ind = zeros((2,3,4), intp)
        >>> _indexShape(index=NUMERIX.index_exp[...,ind,:], 
        ...             arrayShape=(10,20,30))
        (10, 2, 3, 4, 30)
        
    "If the index subspaces are right next to each other, then the broadcasted
    indexing space directly replaces all of the indexed subspaces in X."
       
        >>> ind1 = zeros((2,3,4), intp)
        >>> ind2 = zeros((2,3,4), intp)
        >>> _indexShape(index=NUMERIX.index_exp[:,ind1,ind2], 
        ...             arrayShape=(10,20,30,40,50))
        (10, 2, 3, 4, 40, 50)
        
    "If the indexing subspaces are separated (by slice objects), then the
    broadcasted indexing space is first, followed by the sliced subspace of X."
    
        >>> _indexShape(index=NUMERIX.index_exp[:,ind1,:,ind2,:], 
        ...             arrayShape=(10,20,30,40,50))
        (2, 3, 4, 10, 30, 50)
        
        
    !!! What about Boolean selections ???
    """
    # "when the selection object is not a tuple, it will be referred to as if it
    # had been promoted to a 1-tuple, which will be called the selection tuple"
    if not isinstance(index, tuple):
        if isinstance(index, list):
            index = tuple(index)
        else:
            index = (index,)
    
    Nnewaxes = len([element for element in index if element is newaxis])
    desiredRank = len(arrayShape) + Nnewaxes
    
    ellipsislen = 1 + desiredRank - len(index)
    
    # "Exactly one Ellipsis object will be expanded, any other Ellipsis objects
    # will be treated as full slice (':') objects. The Ellipsis object is replaced
    # with as many full slice (':') objects as needed to make the length of the
    # selection tuple N."
    expanded = ()
    for element in index:
        if element is Ellipsis:
            expanded += (slice(None, None, None),) * ellipsislen
            ellipsislen = 1
        else:
            expanded += (element,)
            
    if len(expanded) > desiredRank:
        # "If the lenth of the selection tuple is larger than N (=X.ndim) an error 
        # is raised."
        if len(arrayShape) == 0:
            raise IndexError, "0-d arrays can't be indexed"
        else:
            raise IndexError, "invalid index"
    else:
        # "If the selection tuple is smaller than N, then as many ':' objects as
        # needed are added to the end of the selection tuple so that the modified 
        # selection tuple has length N."
        expanded += (slice(None, None, None),) * (desiredRank - len(expanded))

    # "The shape of all the integer indexing arrays must be broadcastable to the
    # same shape"
    broadcasted = ()
    arrayIndices = ()
    arrayindex = None
    broadcastshape = None
    i = 0
    j = 0
    while i < len(expanded):
        element = expanded[i]
        if element is newaxis or isinstance(element, slice):
            broadcasted += (element,)
            if isinstance(element, slice):
                arrayIndices += (j,)
        else:
            if broadcastshape is None:
                arrayindex = i
                broadcastshape, skip = _compressIndexSubspaces(index=expanded, i=i)
            else:
                # we're only in this branch if indexing subspaces are separated
                # by slice objects, so indexing subspace is broadcast first
                arrayindex = 0
                broadcastshape, skip = _compressIndexSubspaces(index=expanded, i=i, 
                                                               broadcastshape=broadcastshape)
            i += skip
            j += skip
                
        i += 1
        if element is not newaxis:
            j += 1
        
    indexShape = ()
    j = 0
    for element in broadcasted:
        if element is newaxis:
            indexShape += (1,)
        elif isinstance(element, slice):
            start, stop, stride = element.indices(arrayShape[arrayIndices[j]])
            indexShape += ((stop - start) / stride,)
            j += 1
        else:
            raise IndexError, "invalid index"
                
    if arrayindex is not None:
        indexShape = indexShape[:arrayindex] + broadcastshape + indexShape[arrayindex:]

    return indexShape
    
def _broadcastShapes(shape1, shape2):
    """
    Determine if `shape1` and `shape2` can broadcast to each other, padding if
    necessary, and return their (padded) shapes and the broadcast shape. If the
    shapes cannot broadcast, return a broadcastshape of `None`.
    """
    if len(shape1) > len(shape2):
        shape2 = (1,) * (len(shape1) - len(shape2)) + shape2
    elif len(shape1) < len(shape2):
        shape1 = (1,) * (len(shape2) - len(shape1)) + shape1
    
    if logical_and.reduce([(s == o or s == 1 or o == 1) for s,o in zip(shape1, shape2)]):
        broadcastshape = tuple([max(s,o) for s,o in zip(shape1, shape2)])
    else:
        broadcastshape = None

    return (shape1, shape2, broadcastshape)
    
def _broadcastShape(shape1, shape2):
    """
    Determine if `shape1` and `shape2` can broadcast to each other, padding if
    necessary, and return the broadcast shape. If the shapes cannot broadcast,
    return `None`.
    """
    
    shape1, shape2, broadcastshape = _broadcastShapes(shape1, shape2)
    return broadcastshape
    
if not hasattr(NUMERIX, "in1d"):
    # this handy function was introduced at some point (but it's not in 1.4.1)
    # we define if necessary
    
    def in1d(ar1, ar2, assume_unique=False):
        """
        Test whether each element of a 1D array is also present in a second array.

        Returns a boolean array the same length as `ar1` that is True
        where an element of `ar1` is in `ar2` and False otherwise.

        Parameters
        ----------
        ar1 : array_like, shape (M,)
            Input array.
        ar2 : array_like
            The values against which to test each value of `ar1`.
        assume_unique : bool, optional
            If True, the input arrays are both assumed to be unique, which
            can speed up the calculation.  Default is False.

        Returns
        -------
        mask : ndarray of bools, shape(M,)
            The values `ar1[mask]` are in `ar2`.

        See Also
        --------
        numpy.lib.arraysetops : Module with a number of other functions for
                                performing set operations on arrays.

        Notes
        -----
        `in1d` can be considered as an element-wise function version of the
        python keyword `in`, for 1D sequences. ``in1d(a, b)`` is roughly
        equivalent to ``np.array([item in b for item in a])``.

        .. versionadded:: NOT 1.4.0 !!!!

        Examples
        --------
        >>> test = NUMERIX.array([0, 1, 2, 5, 0])
        >>> states = [0, 2]
        >>> mask = in1d(test, states)
        >>> mask
        array([ True, False,  True, False,  True], dtype=bool)
        >>> test[mask]
        array([0, 2, 0])

        """
        if not assume_unique:

            try:
                ar1, rev_idx = NUMERIX.unique(ar1, return_inverse=True)
            except TypeError:
                ar1, rev_idx = NUMERIX.unique1d(ar1, return_inverse=True)

            ar2 = NUMERIX.unique(ar2)

        ar = NUMERIX.concatenate( (ar1, ar2) )
        # We need this to be a stable sort, so always use 'mergesort'
        # here. The values from the first array should always come before
        # the values from the second array.
        order = ar.argsort(kind='mergesort')
        sar = ar[order]
        equal_adj = (sar[1:] == sar[:-1])
        flag = NUMERIX.concatenate( (equal_adj, [False] ) )
        indx = order.argsort(kind='mergesort')[:len( ar1 )]

        if assume_unique:
            return flag[indx]
        else:
            return flag[indx][rev_idx]

    
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__":
    _test() 
