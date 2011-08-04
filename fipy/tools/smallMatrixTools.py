#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "smallMatrixTools.py"
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
 #  
 # ###################################################################
 ##

r"""
Vectorized routines for small matrices.

>>> import numpy
>>> numpy.random.seed(1)
>>> A = numpy.random.random((7, 2, 2))
>>> A[0,0,0] = 0
>>> A[0,1,1] = 0
>>> A[1,1,0] = 0
>>> A[1,0,1] = 0
>>> A[2,0,0] = 0
>>> A[3,1,1] = 0
>>> A[4,0,1] = 0
>>> A[5,1,0] = 0
>>> B = numpy.random.random((7, 2, 2))
>>> I = numpy.identity(A.shape[-1])
>>> print numpy.allclose(mul(slowinv(A), B), invmul(A, B))
True
>>> print numpy.allclose(mul(A, slowinv(B)), mulinv(A, B))
True
>>> print numpy.allclose(mulinv(A, A), I[numpy.newaxis])
True
>>> print numpy.allclose(invmul(A, A), I[numpy.newaxis])
True
>>> e, R = sortedeig(A)
>>> eslow, Rslow = sortedeig(A, slow=True)
>>> print numpy.allclose(e, eslow)
True
>>> print numpy.allclose(mul(A, R), e[:,numpy.newaxis] * R)
True
>>> E = e[:,:,numpy.newaxis] * I
>>> print numpy.allclose(A, mulinv(mul(R, E), R))
True

>>> A = numpy.random.random((5, 3, 3))
>>> B = numpy.random.random((5, 3, 3))
>>> I = numpy.identity(A.shape[-1])
>>> print numpy.allclose(mulinv(A, A), I[numpy.newaxis])
True
>>> print numpy.allclose(invmul(A, A), I[numpy.newaxis])
True
>>> e, R = sortedeig(A)
>>> eslow, Rslow = sortedeig(A, slow=True)
>>> print numpy.allclose(e, eslow)
True
>>> print numpy.allclose(R / R[:,0][:,numpy.newaxis], Rslow / Rslow[:,0][:,numpy.newaxis])
True
>>> E = e[:,:,numpy.newaxis] * I
>>> print numpy.allclose(A, mulinv(mul(R, E), R))
True

>>> A = numpy.random.random((5, 4, 4))
>>> B = numpy.random.random((5, 4, 4))
>>> I = numpy.identity(A.shape[-1])
>>> print numpy.allclose(mul(slowinv(A), B), invmul(A, B))
True
>>> print numpy.allclose(mul(A, slowinv(B)), mulinv(A, B))
True
>>> print numpy.allclose(mulinv(A, A), I[numpy.newaxis])
True
>>> print numpy.allclose(invmul(A, A), I[numpy.newaxis])
True
>>> e, R = sortedeig(A)
>>> eslow, Rslow = sortedeig(A, slow=True)
>>> print numpy.allclose(e, eslow)
True
>>> print numpy.allclose(R / R[:,0][:,numpy.newaxis], Rslow / Rslow[:,0][:,numpy.newaxis])
True
>>> E = e[:,:,numpy.newaxis] * I
>>> print numpy.allclose(A, mulinv(mul(R, E), R))
True

"""

__docformat__ = 'restructuredtext'

import scipy
import scipy.sparse
import scipy.sparse.linalg
import numpy

def slowinv(A): 
    """ 
    Inverts N, MxM matrices, A.shape = (M, M, N). 
    """ 
    return numpy.array(map(numpy.linalg.inv, A)) 

def fastsum(arr, axis=0):
    if type(arr) in (float, int) or len(arr) == 0 or 0 in arr.shape:
        return numpy.sum(arr, axis)
    else:
        return numpy.tensordot(numpy.ones(arr.shape[axis], 'l'), arr, (0, axis))

def mul(A, B):
    """
    Matrix multiply N, MxM matrices, A.shape = B.shape = (N, M, M).
    """
    return fastsum(A.swapaxes(1,2)[:, :, :, numpy.newaxis] * B[:, :, numpy.newaxis], 1)

def invmul(A, B):
    """
    Calculates the product (A^-1.B) of N, M (rows) x P (cols) matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (N, M, N), B.shape = (N, M, P).
    """
    N, M, P = B.shape
    N, M, M = A.shape
    import scipy
    Amat = scipy.sparse.bsr_matrix((A, numpy.arange(N), numpy.arange(N + 1)), shape=(M * N, M * N))
    ans = numpy.zeros((N, M, P), 'd')

    import scipy
    Amat = scipy.sparse.bsr_matrix((A, numpy.arange(N), numpy.arange(N + 1)), shape=(M * N, M * N))
    LU = scipy.sparse.linalg.splu(Amat.asformat("csc"))
    for iCol in range(P):
        ans[:,:,iCol] = LU.solve(B[:,:,iCol].ravel().real).reshape((N, M))

    return ans

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
    return invmul(A, v[:,:,numpy.newaxis]).reshape(v.shape)
    
def sloweigvec(A):
    """
    Calculate the eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    tmp = zip(*map(numpy.linalg.eig, A))
    return numpy.array(tmp[0]), numpy.array(tmp[1])

def matvec(A, v):
    """
    Multiply N, MxM matrices by N, M length vectors, A.shape = (M, M, N), v.shape = (M, N).
    """
    return fastsum(A * v, 1)

def polar(x, y):
    ## From http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html#Sec2.2
    return scipy.hypot(x, y), scipy.arctan2(y, x) 

def quadratic(a, b, c):
    ## From http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html#Sec2.2
    t = b / a / 2.0
    r = t**2 - c / a
    y1 = scipy.sqrt(r)
    y2 = -y1
    return y1 - t, y2 - t 

def cbrt(x):
    return x**(1./ 3.)

def cubic(a, b, c, d):
    ## From http://www.physics.rutgers.edu/~masud/computing/WPark_recipes_in_python.html#Sec2.2
    A, B, C = b / float(a), c / float(a), d / float(a)
    t = A / 3.0
    p, q = B - 3 * t**2, C - B * t + 2 * t**3
    u, v = quadratic(1., q, -(p/3.0)**3)

    r, w = polar(u.real, u.imag)
    y1 = numpy.where(scipy.iscomplex(u),
                     2 * cbrt(r) * scipy.cos(w / 3.0),
                     cbrt(u) + cbrt(v))
                     
    y2, y3 = quadratic(1., y1, p + y1**2)
    return y1 - t, y2 - t, y3 - t 

def determinant(A):
    N, M, M = A.shape
    if M == 1:
        det = A[:,0,0]
    else:
        det = 0
        for iCol in range(M):
            value = A[:,0,iCol] * determinant(numpy.concatenate((A[:,1:,:iCol], A[:,1:,iCol + 1:]), axis=2))
            if iCol % 2 == 0:
                det += value
            else:
                det -= value
    return det

def trace(A):
    N, M, M = A.shape
    tr = 0
    for iCol in range(M):
        tr += A[:,iCol, iCol]
    return tr

## def vecN(A, eigs):
##     """
##     Calculates the eigenvectors for N, MxM matrices (A) using the eigenvalues (eigs).
##     A.shape = (N, M, M), eigs.shape = (N, M).
##     """
##     N, M, M = A.shape
##     vectors = numpy.zeros((N, M, M), 'd')
##     for iCol in range(M):
##         Amod = A - eigs[:,iCol,numpy.newaxis,numpy.newaxis] * numpy.identity(M)
##         Amat = scipy.sparse.bsr_matrix((Amod, numpy.arange(N), numpy.arange(N + 1)), shape=(M * N, M * N))
##         X = numpy.ones((N, M), 'd') * 0.1
##         X[:,iCol] = 1.
##         vectors[:,:,iCol] = scipy.sparse.linalg.gmres(Amat.asformat("csc"), numpy.zeros(M * N), X.ravel())[0].reshape((N, M))
##         vectors[:,:,iCol] = vectors[:,:,iCol] / numpy.sqrt(numpy.sum(vectors[:,:,iCol]**2, 1))[...,numpy.newaxis]
##     return vectors

## def vecN(A, eigs):
##     """
##     Calculates the eigenvectors for N, MxM matrices (A) using the eigenvalues (eigs).
##     A.shape = (N, M, M), eigs.shape = (N, M).
##     """
##     N, M, M = A.shape
##     vectors = numpy.zeros((N, M, M), 'd')
##     for iCol in range(M):
##         Amod = A - eigs[:,iCol,numpy.newaxis,numpy.newaxis] * numpy.identity(M)
##         Amat = scipy.sparse.bsr_matrix((Amod, numpy.arange(N), numpy.arange(N + 1)), shape=(M * N, M * N))
##         X = numpy.zeros((N, M), 'd')
##         X[:,iCol] = 1.
##         print Amat
##         LU = scipy.sparse.linalg.splu(Amat.asformat("csc"))
##         error = matvec(L, X)
##         vectors[:,:,iCol] -= LU.solve(error.ravel()).reshape((N, M))
## ##        vectors[:,:,iCol] = scipy.sparse.linalg.gmres(Amat.asformat("csc"), numpy.zeros(M * N), X.ravel())[0].reshape((N, M))
##         vectors[:,:,iCol] = vectors[:,:,iCol] / numpy.sqrt(numpy.sum(vectors[:,:,iCol]**2, 1))[...,numpy.newaxis]
##     return vectors

def eigvec2(A):
    """
    Calculates the eigenvalues and eigenvectors of N, 2x2 matrices, A.shape = (N, 2, 2).
    """
    N, M, M = A.shape

    eigs = numpy.zeros((N, M), dtype=complex)
    eigs[:,0], eigs[:,1] = quadratic(1., -trace(A), determinant(A))

    vectors = numpy.zeros((N, M, M), dtype=complex)
    a, b, c, d, = A[:,0,0][:,numpy.newaxis], A[:,0,1][:,numpy.newaxis], A[:,1,0][:,numpy.newaxis], A[:,1,1][:,numpy.newaxis]

    tmp = numpy.zeros((N, M, M), dtype=complex)
    tmp[:,0,0] = 1
    tmp[:,1,1] = 1
    tmp[:,0,1] = 0
    tmp[:,1,0] = 0

    bb = numpy.concatenate((b[...,numpy.newaxis], b[...,numpy.newaxis]), axis=2)
    cc = numpy.concatenate((c[...,numpy.newaxis], c[...,numpy.newaxis]), axis=2)
    
    vectors = numpy.where(bb == 0,
                          numpy.where(cc == 0,
                                      numpy.where((a - eigs[:,0][:,numpy.newaxis])[:,numpy.newaxis,:] == 0,
                                                  tmp,
                                                  tmp[:,::-1]),
                                      numpy.concatenate(((eigs - d).reshape(N, 1, M), cc), axis=1)),
                          numpy.concatenate((bb, (eigs - a).reshape(N, 1, M)), axis=1))
                         
    vectors[:,:,:] = vectors / numpy.sqrt(fastsum(vectors**2, 1))[:,numpy.newaxis, :]

    return eigs, vectors

def eigvec3(A):
    """
    Calculates the eigenvalues and eigenvectors of N, 3x3 matrices, A.shape = (N, 3, 3).
    """
    N, M, M = A.shape

    eigs = numpy.zeros((N, M), dtype=complex)
    trA = trace(A)
    eigs[:,0], eigs[:,1], eigs[:,2] = cubic(-1., trA, (trace(mul(A, A)) - trA**2) / 2, determinant(A))

    vectors = numpy.zeros((N, M, M), dtype=complex)
    vectors[:, 0, :] = 1.
    a, b, c = A[:,0,0], A[:,0,1], A[:,0,2]
    d, e, f = A[:,0,0], A[:,0,1], A[:,0,2]
    vectors[:, 1, :] = 1
    d0 = (c * e - b * f)[:,numpy.newaxis]
    d1 = (a * f - c * d)[:,numpy.newaxis]
    d2 = (a * e - b * d)[:,numpy.newaxis]
    
    vectors[:, 1, :] = (eigs * (c - f)[:,numpy.newaxis] + d1) / d0
    vectors[:, 2, :] = (eigs * (b - e)[:,numpy.newaxis] + d2) / -d0
    vectors[:,:,:] = vectors / numpy.sqrt(fastsum(vectors**2, 1))[:,numpy.newaxis, :]

    return eigs, vectors

def eigvec(A, slow=False):
    N, M, M = A.shape
    if M == 2 and not slow:
        return eigvec2(A)
    elif M == 3 and not slow:
        return eigvec3(A)
    else:
        return sloweigvec(A)

def sortedeig(A, slow=False):
    """
    Caclulates the sorted eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    N = A.shape[0]
    eigenvalues, R = eigvec(A, slow=slow)
    order = eigenvalues.argsort(1)
    Nlist = numpy.arange(N).reshape((N, 1))
    return eigenvalues[Nlist, order], R[Nlist, :, order].transpose(0, 2, 1)

def dot(v0, v1):
    """
    Dot product of two N, M-length vectors, v0.shape = (N, M), v1.shape = (N, M).
    """
    return fastsum(v0 * v1, 1)


def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
