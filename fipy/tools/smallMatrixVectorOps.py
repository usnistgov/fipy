#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "smallMatrixVectorOps.py"
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
>>> import numpy as np

>>> np.random.seed(1)
>>> N = 5
>>> def test(M):
...     A = -np.random.random((N, M, M))
...     for i in range(M):
...         A[i, i] += 1 + M
...     B = np.random.random((N, M, M))
...     I = np.identity(A.shape[-1])
...
...     e, R = sortedeig(A)
...     eslow, Rslow = sortedeig(A, slow=True)
...     E = e[:,:,np.newaxis] * I
...     return np.all([np.allclose(invmul(A, A), I[np.newaxis]),
...                    np.allclose(e, eslow),
...                    np.allclose(mul(A, R), e[:,np.newaxis] * R),
...                    np.allclose(A, mulinv(mul(R, E), R))])

>>> b = []
>>> for M in (2, 3, 4, 5):
...     b += [test(M)]
>>> np.all(b)
True

"""

import numpy as np
import pyximport
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

pyximport.install(build_dir=os.path.join(os.path.expanduser('~'), '.pyxbld'),
                  setup_args = {'options' : {'build_ext' : {'libraries' : 'lapack'}}})
                             
from smallMatrixVectorOpsExt import solve, fasteigvec

def fastsum(arr, axis=0):
    if type(arr) in (float, int) or len(arr) == 0 or 0 in arr.shape:
        return np.sum(arr, axis)
    else:
        return np.tensordot(np.ones(arr.shape[axis], 'l'), arr, (0, axis))

def mul(A, B):
    """
    Matrix multiply N, MxM matrices, A.shape = B.shape = (N, M, M).
    """
    return fastsum(A.swapaxes(1,2)[:, :, :, np.newaxis] * B[:, :, np.newaxis], 1)

def slowinv(A): 
    """ 
    Inverts N, MxM matrices, A.shape = (M, M, N). 
    """ 
    return np.array(map(np.linalg.inv, A)) 

def invmul(A, B):
    """
    Calculates the product (A^-1.B) of N, M (rows) x P (cols) matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (N, M, N), B.shape = (N, M, P).
    """
    N, M, P = B.shape
    N, M, M = A.shape

    return solve(A, B).transpose(0, 2, 1)

def mulinv(A, B):
    """
    Calculates the product (A.B^-1) of N, MxM matrices (A) and the inverse of N, MxM matrices (B).
    A.shape = (N, M, M), B.shape = (N, M, M).
    """
    return invmul(B.transpose(0, 2, 1), A.transpose(0, 2, 1)).transpose(0, 2, 1)

def invmatvec(A, v):
    """
    Calculates the product (A^-1.v) of N, MxM matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (N, M, M), v.shape = (N, M).
    """
    return invmul(A, v[:,:,np.newaxis]).reshape(v.shape)

def matvec(A, v):
    """
    Multiply N, MxM matrices by N, M length vectors, A.shape = (M, M, N), v.shape = (M, N).
    """
    return fastsum(A * v, 1)

def dot(v0, v1):
    """
    Dot product of two N, M-length vectors, v0.shape = (N, M), v1.shape = (N, M).
    """
    return fastsum(v0 * v1, 1)

def sloweigvec(A):
    """
    Calculate the eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    tmp = zip(*map(np.linalg.eig, A))
    return np.array(tmp[0]), np.array(tmp[1])

def eigvec(A, slow=False):
    N, M, M = A.shape
    if slow:
        return sloweigvec(A)
    else:
        return fasteigvec(A)

def sortedeig(A, slow=False):
    """
    Caclulates the sorted eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    N = A.shape[0]
    eigenvalues, R = eigvec(A, slow=slow)
    order = eigenvalues.argsort(1)
    Nlist = np.arange(N).reshape((N, 1))
    return eigenvalues[Nlist, order], R[Nlist, :, order].transpose(0, 2, 1)

def _test(): 
    import doctest
    return doctest.testmod()
 
if __name__ == "__main__":
    _test()
