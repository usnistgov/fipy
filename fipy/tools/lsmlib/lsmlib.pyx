import cython
import numpy as np
from numpy.core import intc	
cimport numpy as np
from numpy.compat import asbytes
from numpy import int, double, int32
from lsmlib cimport computeDistanceFunction2d
from libc.stdlib cimport malloc, free

def computeDistanceFunction2d_(np.ndarray[double, ndim=1] phi, nx=1, ny=1, dx=1., dy=1., order=2):
    cdef np.ndarray[double, ndim=1] distance_function = np.zeros((nx * ny,))
    cdef np.ndarray[double, ndim=1] mask = np.zeros((nx * ny,))
    cdef int spatial_derivative_order = order				
    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy))
    
    error = computeDistanceFunction2d(
	<double *> distance_function.data,
    	<double *> phi.data,	
    	<double *> mask.data,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    return distance_function

def computeExtensionFields2d_(np.ndarray[double, ndim=1] phi, np.ndarray[double, ndim=2] extensionFields, nx=1, ny=1, dx=1., dy=1.):

    cdef int num_ext_fields = extensionFields.shape[0]

    cdef np.ndarray[double, ndim=1] distance_function = np.zeros((nx * ny,))
    cdef np.ndarray[double, ndim=1] mask = np.zeros((nx * ny,))
    cdef int spatial_derivative_order = 2							
    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy))
    
    cdef double **ext_fields = <double **> malloc(num_ext_fields * sizeof(double*))
    cdef double **source_fields = <double **> malloc(num_ext_fields * sizeof(double*))
    cdef np.ndarray[double, ndim=2] extReturnFields = np.zeros((num_ext_fields, nx * ny))

    for i in range(num_ext_fields):
        ext_fields[i] = &extReturnFields[i,0]
        source_fields[i] = &extensionFields[i,0]
        
    error = computeExtensionFields2d(
	<double *> distance_function.data,
	<double **> ext_fields,
    	<double *> phi.data,	
    	<double *> mask.data,
	<double **> source_fields,
	num_ext_fields,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    free(ext_fields)
    free(source_fields)

    return distance_function, extReturnFields

def solveEikonalEquation2d_(np.ndarray[double, ndim=1] phi,
                            np.ndarray[double, ndim=1] speed, nx=1, ny=1, dx=1., dy=1.):

    cdef np.ndarray[double, ndim=1] mask = np.zeros((nx * ny,))
    cdef int spatial_derivative_order = 2							
    cdef np.ndarray[int, ndim=1] grid_dims = np.array((nx, ny), dtype=int32)
    cdef np.ndarray[double, ndim=1] _dx = np.array((dx, dy))

    minphi = min(phi)
    phi[phi > 0] = minphi - 1.
    phi[:] = phi - minphi

    error = solveEikonalEquation2d(
	<double *> phi.data,
    	<double *> speed.data,	
    	<double *> mask.data,
    	spatial_derivative_order,
    	<int *> grid_dims.data,
    	<double *> _dx.data)

    phi[:] = phi + minphi

    return phi
