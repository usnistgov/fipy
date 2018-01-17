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

"""Replacement module for NumPy

.. attention::

   This module should be the only place in the code where :mod:`numpy` is
   explicitly imported and you should always import this module and not
   :mod:`numpy` in your own code. The documentation for :mod:`numpy` remains
   canonical for all functions and classes not explicitly documented here.

The functions provided in ths module replace and augment the `NumPy` module.
The functions work with `Variables`, arrays or numbers. For example,
create a `Variable`.

   >>> from fipy.variables.variable import Variable
   >>> var = Variable(value=0)

Take the tangent of such a variable. The returned value is itself a
`Variable`.

   >>> v = tan(var)
   >>> v
   tan(Variable(value=array(0)))
   >>> print float(v)
   0.0

Take the tangent of a int.

   >>> tan(0)
   0.0

Take the tangent of an array.

   >>> print tan(array((0,0,0)))
   [ 0.  0.  0.]

"""

__docformat__ = 'restructuredtext'

import numpy as NUMERIX

# On all platforms except Win64, int will be either 32 bit
# or 64 bit depending on the platform and the Python;
# however, on Win64 (and python 64-bit), int will always be
# 32-bit. This is a known issue, and StackOverflow has many
# questions and answers.  However, we still need to do
# something about it, so instead of relying on the default
# mapping, we'll be explicit:

import platform
arch=platform.architecture()[0]
if arch == '32bit':
    INT_DTYPE=NUMERIX.int32
elif arch == '64bit':
    INT_DTYPE=NUMERIX.int64
else:
    raise Exception('Cannot set integer dtype because architecture is unknown.')

from numpy.core import umath
from numpy import newaxis as NewAxis
from numpy import *
try:
    from numpy.core import ma as MA
    numpy_version = 'old'
except ImportError:
    # masked arrays have been moved in numpy 1.1
    from numpy import ma as MA
    numpy_version = 'new'

from fipy.tools import inline

# we want NumPy's __all__, with adjustments
import sys
__all__ = list(sys.modules['numpy'].__dict__.setdefault('__all__', []))
__all__.extend(["NUMERIX", "NewAxis", "MA", "numpy_version"])
__all__.extend(sorted(["getUnit", "put", "reshape", "getShape",
                       "rank", "sum", "isFloat", "isInt", "tostring", "dot",
                       "sqrtDot", "nearest", "allequal", "allclose", "all",
                       "isclose", "take", "indices", "empty", "loadtxt",
                       "savetxt", "L1norm", "L2norm", "LINFnorm", "in1d"],
                      key=str.lower))

def _isPhysical(arr):
    """
    Returns `True` if arr is a `Variable` or `PhysicalField`.
    """
    from fipy.variables.variable import Variable
    from fipy.tools.dimensions.physicalField import PhysicalField

    return isinstance(arr,Variable) or isinstance(arr,PhysicalField)

def getUnit(arr):
    if hasattr(arr, "getUnit") and callable(arr.getUnit):
        return arr.unit
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

        shape = left + (oldShape.prod() // newShape.prod(),) + right

    if _isPhysical(arr):
        return arr.reshape(shape)
    elif type(arr) is type(array((0))):
        return NUMERIX.reshape(arr, tuple(shape))
    elif type(arr) is type(MA.array((0))):
        return MA.reshape(arr, shape)
    else:
        return NUMERIX.reshape(array(arr), tuple(shape))
##        raise TypeError, 'cannot reshape object ' + str(arr)

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
    if hasattr(a, "rank"):
        return a.rank
    elif hasattr(a, "getRank"):
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
        if type(arr) in (float, int) or len(arr) == 0 or 0 in arr.shape:
            return NUMERIX.sum(arr, axis)
        else:
            if axis is None:
                axis = 0
            return NUMERIX.tensordot(NUMERIX.ones(arr.shape[axis], 'l'), arr, (0, axis))

def isFloat(arr):
    if isinstance(arr, NUMERIX.ndarray):
        return NUMERIX.issubclass_(arr.dtype.type, NUMERIX.floating)
    else:
        return NUMERIX.issubclass_(arr.__class__, float)

def isInt(arr):
    if isinstance(arr, NUMERIX.ndarray):
        return NUMERIX.issubclass_(arr.dtype.type, NUMERIX.integer)
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
            ## this is for numpy 1.14 and above
            ## why has the interface changed *again*?
            from numpy.core.arrayprint import FloatingFormat
            return FloatingFormat(data=NUMERIX.array((arr,)), precision=precision,
                                  floatmode='maxprec', suppress_small=suppress_small)(arr).strip()
        except ImportError:
            try:
                ## this is for numpy 1.0.4 and above
                ## why has the interface changed again?
                from numpy.core.arrayprint import FloatFormat
                return FloatFormat(NUMERIX.array((arr,)), precision, suppress_small)(arr).strip()
            except ImportError:
                from numpy.core.arrayprint import _floatFormat, _formatFloat
                return _formatFloat(arr, format='%%1.%df' % precision)

    elif isInt(arr):
        try:
            ## this is for numpy 1.7 and above
            ## why has the interface changed again?
            from numpy.core.arrayprint import IntegerFormat
            return IntegerFormat(NUMERIX.array((arr,)))(arr).strip()
        except:
            from numpy.core.arrayprint import _formatInteger
            return _formatInteger(arr, format='%d')
    else:
        raise TypeError, 'cannot convert ' + str(arr) + ' to string'

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

    >>> from fipy.meshes import Grid2D
    >>> mesh = Grid2D(nx=2, ny=1)
    >>> from fipy.variables.cellVariable import CellVariable
    >>> v1 = CellVariable(mesh=mesh, value=((0,1),(2,3)), rank=1)
    >>> v2 = CellVariable(mesh=mesh, value=((0,1),(2,3)), rank=1)
    >>> dot(v1, v2)._variableClass
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> dot(v2, v1)._variableClass
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> print rank(dot(v2, v1))
    0
    >>> print dot(v1, v2)
    [ 4 10]
    >>> dot(v1, v1)._variableClass
    <class 'fipy.variables.cellVariable.CellVariable'>
    >>> print dot(v1, v1)
    [ 4 10]
    >>> v3 = array(((0,1),(2,3)))
    >>> print type(dot(v3, v3)) is type(array(1))
    1
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

if inline.doInline:
    def sqrtDot(a1, a2):
        """Return array of square roots of vector dot-products
        for arrays a1 and a2 of vectors v1 and v2

        Usually used with v1==v2 to return magnitude of v1.
        """
        unit1 = unit2 = 1

        def dimensionlessUnmasked(a):
            unit = 1
            mask = False
            if _isPhysical(a):
                unit = a.inBaseUnits().getUnit()
                a = a.numericValue
            if MA.isMaskedArray(a):
                mask = a.mask
                a = a.filled(fill_value=1)
            if not a.flags['C_CONTIGUOUS']:
                a = a.copy('C')
            return (a, unit, mask)

        a1, unit1, mask1 = dimensionlessUnmasked(a1)
        a2, unit2, mask2 = dimensionlessUnmasked(a2)

        NJ, ni = NUMERIX.shape(a1)
        result1 = NUMERIX.zeros((ni,),'d')

        inline._runInline("""
            int j;
            result1[i] = 0.;
            for (j = 0; j < NJ; j++)
            {
                result1[i] += a1[i + j * ni] * a2[i + j * ni];
            }
            result1[i] = sqrt(result1[i]);
        """,result1=result1, a1=a1, a2=a2, ni=ni, NJ=NJ)

        if NUMERIX.any(mask1) or NUMERIX.any(mask2):
            result1 = MA.array(result1, mask=NUMERIX.logical_or(mask1, mask2))

        if unit1 != 1 or unit2 != 1:
            from fipy.tools.dimensions.physicalField import PhysicalField
            result1 = PhysicalField(value=result, unit=(unit1 * unit2)**0.5)

        return result1
else:
    def sqrtDot(a1, a2):
        """Return array of square roots of vector dot-products
        for arrays a1 and a2 of vectors v1 and v2

        Usually used with v1==v2 to return magnitude of v1.
        """
        ## We can't use Numeric.dot on an array of vectors
        return sqrt(dot(a1, a2))

def nearest(data, points, max_mem=1e8):
    """find the indices of `data` that are closest to `points`

    >>> from fipy import *
    >>> m0 = Grid2D(dx=(.1, 1., 10.), dy=(.1, 1., 10.))
    >>> m1 = Grid2D(nx=2, ny=2, dx=5., dy=5.)
    >>> print nearest(m0.cellCenters.globalValue, m1.cellCenters.globalValue)
    [4 5 7 8]
    >>> print nearest(m0.cellCenters.globalValue, m1.cellCenters.globalValue, max_mem=100)
    [4 5 7 8]
    >>> print nearest(m0.cellCenters.globalValue, m1.cellCenters.globalValue, max_mem=10000)
    [4 5 7 8]
    """
    data = asanyarray(data)
    points = asanyarray(points)

    D = data.shape[0]
    N = data.shape[-1]
    M = points.shape[-1]

    if N == 0:
        return arange(0)

    # given (D, N) data and (D, M) points,
    # break points into (D, C) chunks of points
    # calculatate the full factorial (D, N, C) distances between them
    # and then reduce to the indices of the C closest values of data
    # then assemble chunks C into total M closest indices

    # (D, N) -> (D, N, 1)
    data = data[..., newaxis]

    # there appears to be no benefit to taking chunks that use more
    # than about 100 MiB for D x N x C, and there is a substantial penalty
    # for going much above that (presumably due to swapping, even
    # though this is vastly less than the 4 GiB I had available)
    # see ticket:348

    numChunks = int(round(D * N * data.itemsize * M / max_mem + 0.5))

    nearestIndices = empty((M,), dtype=INT_DTYPE)
    for chunk in array_split(arange(points.shape[-1]), numChunks):
        # last chunk can be empty, but numpy (1.5.0.dev8716, anyway)
        # returns array([], dtype=float64), which can't be used for indexing
        chunk = chunk.astype(int)

        # (D, M) -> (D, C)
        chunkOfPoints = points[..., chunk]
        # (D, C) -> (D, 1, C)
        chunkOfPoints = chunkOfPoints[..., newaxis, :]
        # (D, 1, C) -> (D, N, C)
        chunkOfPoints = NUMERIX.repeat(chunkOfPoints, N, axis=1)

#         print "chunkOfPoints size: ", chunkOfPoints.shape, chunkOfPoints.size, chunkOfPoints.itemsize, chunkOfPoints.size * chunkOfPoints.itemsize

        try:
            tmp = data - chunkOfPoints
        except TypeError:
            tmp = data - PhysicalField(chunkOfPoints)

        # (D, N, C) -> (N, C)
        tmp = dot(tmp, tmp, axis=0)

        # (N, C) -> C
        nearestIndices[chunk] = argmin(tmp, axis=0)

    return nearestIndices

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

def all(a, axis=None, out=None):
    r"""
    Test whether all array elements along a given axis evaluate to True.

    Parameters
    ----------
    a : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which an logical AND is performed.
        The default (`axis` = `None`) is to perform a logical AND
        over a flattened input array. `axis` may be negative, in which
        case it counts from the last to the first axis.
    out : ndarray, optional
        Alternative output array in which to place the result.
        It must have the same shape as the expected output and
        the type is preserved.

    """
    if _isPhysical(a):
        return a.all(axis=axis)
    else:
        return MA.all(a=a, axis=axis, out=out)

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

if not hasattr(NUMERIX, 'empty'):
    print 'defining empty'
    if inline.doInline:
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
            import weave

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
    else:
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
            return NUMERIX.zeros(shape, dtype)

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
        ...             arrayShape=(10,20,30,40,50))            #doctest: +IGNORE_EXCEPTION_DETAIL
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
            indexShape += ((stop - start) // stride,)
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

    Broadcasting zero lenght arrays must also be accounted for.

    >>> _broadcastShapes((1,), (0,))[2]
    (0,)
    >>> _broadcastShapes((2, 0,), (1,))[2]
    (2, 0)

    """

    if len(shape1) > len(shape2):
        shape2 = (1,) * (len(shape1) - len(shape2)) + shape2
    elif len(shape1) < len(shape2):
        shape1 = (1,) * (len(shape2) - len(shape1)) + shape1

    def maxzero(s, o):
        if s == 0 or o == 0:
            return 0
        else:
            return max(s,o)

    if logical_and.reduce([(s == o or s == 1 or o == 1) for s,o in zip(shape1, shape2)]):
        broadcastshape = tuple([maxzero(s,o) for s,o in zip(shape1, shape2)])
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
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
