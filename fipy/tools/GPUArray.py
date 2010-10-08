#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "numerix.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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

import pycuda.gpuarray as gpuarr
import pycuda.driver as cuda
import pycuda.cumath as cumath
import pycuda.autoinit
import numpy
from functools import wraps
from numpy import * # this module sits atop numpy

# Decorators
# ----------


def incomingGPUArrToNumpyArr(method):
    """Decorator which converts an incoming GPUArray to a numpy array. To be
    applied to `GPUArray` methods."""

    def wrapped(self, *args, **kwargs):
        assert isinstance(self, GPUArray), \
                "Must be operating on a GPUArray method."
        return method(self.get(), *args, **kwargs)
    return wraps(method)(wrapped) # preserve method's signature

def ensureArgsOnGPU(method):
    """Decorator which ensures that every argument is a GPUArray. If an argument
    isn't, it will be replaced with an equivalent GPUArray."""

    def wrapped(self, *args, **kwargs):
        newArgs = [a if isinstance(a, GPUArray) else array(a) for a in args]
        return method(self, *newArgs, **kwargs)
    return wraps(method)(wrapped)

def roundTripFromGPU(modifySelf=False):
    """Decorator to apply to any function for which GPUArray doesn't have an
    analogue. 
    Assumes first argument of function is a `GPUArray`.

    If the incoming array is a GPUArray, extract its value from the graphics
    card, make the modification using Numpy functions, and then communicate it
    back out to the card. If `modifySelf` is True, modify `self` instead of
    returning a modified copy of `self`."""

    def withDecArgs(method):
        def wrapped(self, *args, **kwargs):
            assert isinstance(self, GPUArray), \
                    "Must be operating on a GPUArray method."
            result = method(self.get(), *args, **kwargs)
            if not result.flags.contiguous:
                result = result.copy()
            if modifySelf:
                self.set(result)
                return None
            else:
                return array(result)
        return wraps(method)(wrapped)
    return withDecArgs
 
# GPUArray Class
# --------------

class GPUArray(gpuarr.GPUArray):
    """FiPy's very own rip-off of PyCuda's GPUArray class."""

    def __init__(self, *args, **kwargs):

        super(GPUArray, self).__init__(*args, **kwargs)

        self.itemsize = self.dtype.itemsize


    @incomingGPUArrToNumpyArr
    def getNdim(self):
        return self.ndim

    ndim = property(getNdim)

    @roundTripFromGPU()
    def reshape(self, *shape):
        return self.reshape(tuple(shape))

    @incomingGPUArrToNumpyArr
    def getRank(self):
        return numpy.rank(self)

    @incomingGPUArrToNumpyArr
    def max(self):
        return self.max()

    @incomingGPUArrToNumpyArr
    def min(self):
        return self.min()

    @incomingGPUArrToNumpyArr
    def sum(self):
        return self.sum()

    @roundTripFromGPU()
    def put(self, ids, values):
        return numpy.put(self, ids, values)

    def arccos(self):
        return cumath.acos(self)

    @roundTripFromGPU()
    def arccosh(self):
        return numpy.arccosh(self)

    def arcsin(self):
        return cumath.asin(self)

    @roundTripFromGPU()
    def arcsinh(self):
        return numpy.archsinh(self)

    def arctan(self):
        return cumath.atan(self)

    @roundTripFromGPU()
    def arctan2(self, other):
        return numpy.arctan2(self, other)

    @roundTripFromGPU()
    def arctanh(self):
        return numpy.arctanh(self)

    def cos(self):
        return cumath.cos(self)

    def cosh(self):
        return cumath.cosh(self)

    def tan(self):
        return cumath.tan(self)

    def tanh(self):
        return cumath.tanh(self)

    def log10(self):
        return cumath.log10(self)

    def sin(self):
        return cumath.sin(self)

    def sinh(self):
        return cumath.sinh(self)

    def sqrt(self):
        return cumath.sqrt(self)

    def floor(self):
        return cumath.floor(self)

    def ceil(self):
        return cumath.ceil(self)

    @roundTripFromGPU()
    def sign(self):
        return numpy.sign(self)

    def exp(self):
        return cumath.exp(self)

    def log(self):
        return cumath.log(self)

    def conjugate(self):
        return self.conj()

    def take(self, indices):
        return array(gpuarr.take(self, indices))

    @roundTripFromGPU()
    def swapaxes(self, a, b):
        """
        >>> gpua = array([[1, 2, 3], [4, 5, 6]])
        >>> gpua
        array([[1, 2, 3],
               [4, 5, 6]])
        >>> gpua.swapaxes(0,1)
        array([[1, 4],
               [2, 5],
               [3, 6]])
        """
        #tempArr = self.get()
        ##print "  self.swap", self.swapaxes(*args)
        #return array(tempArr.swapaxes(*args))
        return numpy.swapaxes(self, a, b)
    
     # TODO: is the FiPy conventional ``dot'' correct?
    @roundTripFromGPU()
    def dot(self, *others):
        """Compute the scalar product with variable arity."""
        c = numpy.sum(self*others[0], axis=0)
        for o in others[1:]:
            c = numpy.sum(c*o, axis=0)

        return c

    @ensureArgsOnGPU
    def sqrtDot(self, *others):
        """Compute the square root of the scalar product of a variable number of
        arrays."""
        dotted = self.dot(*others)
        return cumath.sqrt(dotted)

    def __getitem__(self, idx):
        """If `idx` is an int, return a scalar a la numpy; else, hand-off
        to parent's method."""
        if type(idx) is int:
            return self.get()[idx]
        elif type(idx) in (numpy.ndarray, tuple):
            return self._numpyIdx(idx)
        elif type(idx) is GPUArray:
            return self._numpyIdx(idx.get())
        else:
            return super(GPUArray, self).__getitem__(idx)

    @roundTripFromGPU()
    def _numpyIdx(self, idx):
        """Used when `self` is indexed with a `numpy.ndarray`."""
        return self[idx]

    @roundTripFromGPU(modifySelf=True)
    def __setitem__(self, idx, value):
        if isinstance(value, GPUArray):
            value = value.get()
        self[idx] = value
        return self

    def __mul__(self, other):
        """pycuda.gpuarray.GPUArray's can't handle (for some reason) direct
        multiplication by large numbers. This can be circumvented by fetching
        the numpy.ndarray equivalent, doing the multiplication, and then
        throwing the result to a GPUArray."""
        try:
            return super(GPUArray, self).__mul__(other)
        except OverflowError:
            temp = self.get()
            return array(temp * other)
 
# NumPy emulation functions
# -------------------------

def array(arr, *args, **kwargs):
    """Return a GPUArray. Argument `arr` can either be a list or a
    `numpy.ndarray`."""
    if type(arr) is GPUArray:
        arr = arr.get()
    else:
        arr = numpy.array(arr, *args, **kwargs)

    # assert isinstance(arr, numpy.ndarray), \
        # "Must be given a list, int, or `numpy.ndarray`; was given `%s`." \
            # % type(arr)

    newArr = GPUArray(arr.shape, dtype=arr.dtype)
    newArr.set(arr)

    return newArr

def zeros(shape, dtype=None):
    """Return a `GPUArray` full of zeros."""
    newArr = GPUArray(shape, dtype=dtype)
    newArr.fill(0)

    return newArr

def ones(shape, dtype=None):
    """Return a `GPUArray` full of ones."""
    newArr = GPUArray(shape, dtype=dtype)
    newArr.fill(1)

    return newArr

def empty(shape, dtype=None):
    return array(numpy.empty(shape, dtype=dtype))

def arange(*args, **kwargs):
    return array(numpy.arange(*args, **kwargs))

def issubclass_(arrType, type):
    return issubclass(arrType, type)

