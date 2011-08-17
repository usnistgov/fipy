#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "clapack.pxd"
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


import numpy as np
cimport numpy as np

cdef extern void dgeev_(char* jobvl, char* jobvr, int* n, double a[],
                        int* lda, double wr[], double wi[], double vl[], int* ldvl,
                        double vr[], int* ldvr, double work[], int* lwork, int* info,
                        int len_jobvl, int len_jobvr)

cdef extern void dgesv_(int* n, int* nrhs, double a[], int* lda, int ipiv[],
                        double b[], int* ldb, int* info)

cdef extern void zgesv_(int* n, int* nrhs, complex a[], int* lda, int ipiv[],
                        complex b[], int* ldb, int* info)

