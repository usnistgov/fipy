import numpy as np
cimport numpy as np

cdef extern int computeDistanceFunction2d(
     double *distance_function,
     double *_phi,
     double *mask,
     int spatial_derivative_order,
     int *grid_dims,
     double *_dx)

cdef extern int computeExtensionFields2d(
     double *distance_function,
     double **ext_fields,
     double *_phi,
     double *mask,
     double **source_fields,
     int num_ext_fields,
     int spatial_derivative_order,
     int *grid_dims,
     double *dX)

cdef extern int solveEikonalEquation2d(
     double *phi,
     double *speed,
     double *mask,
     int spatial_derivative_order,
     int *grid_dims,
     double *dX)
