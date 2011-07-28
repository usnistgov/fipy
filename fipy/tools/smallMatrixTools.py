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

import scipy
import scipy.sparse
import scipy.sparse.linalg
from fipy.tools import numerix

def mul(A, B):
    """
    Matrix multiply N, MxM matrices, A.shape = B.shape = (M, M, N).
    """
    return numerix.sum(A.swapaxes(0,1)[:, :, numerix.newaxis] * B[:, numerix.newaxis], 0)
    
def invmul(A, B):
    """
    Calculates the product (A^-1.B) of N, M (rows) x P (cols) matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (M, M, N), B.shape = (M, M, M).
    """
    tmp = A.transpose(2, 0, 1)
    M, P, N = B.shape
    N, M, M = tmp.shape
    
    import scipy
    Amat = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
    ans = numerix.zeros((N, M, P), 'd')
    Bvec = B.transpose(2, 0, 1)

    LU = scipy.sparse.linalg.splu(Amat.asformat("csc"))    
    for iCol in range(P):
        import scipy.sparse.linalg
        ans[:,:,iCol] = LU.solve(Bvec[:,:,iCol].ravel()).reshape((N, M))

    return ans.transpose(1, 2, 0)

def mulinv(A, B):
    """
    Calculates the product (A.B^-1) of N, MxM matrices (A) and the inverse of N, MxM matrices (B).
    A.shape = (M, M, N), B.shape = (M, M, N).
    """
    return invmul(B.transpose(1, 0, 2), A.transpose(1, 0, 2)).transpose(1, 0, 2)

def invmatvec(A, v):
    """
    Calculates the product (A^-1.v) of N, MxM matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (M, M, N), v.shape = (M, N).
    """
    return invmul(A, v[:,numerix.newaxis,:]).reshape(v.shape)
    
def eigvecN(A):
    """
    Calculate the eigenvalues and eigenvectors of N, MxM matrices, A.shape = (M, M, N).
    """
    tmp = zip(*map(numerix.linalg.eig, A.transpose(2, 0, 1)))
    return numerix.array(tmp[0]).swapaxes(0,1), numerix.array(tmp[1]).transpose(1,2,0)

def vecN(A, eigs):
    """
    Calculates the eigenvectors for N, MxM matrices (A) using the eigenvalues (eigs).
    A.shape = (M, M, N), eigs.shape = (M, N).
    """
    M, M, N = A.shape
    vectors = numerix.zeros((N, M, M), 'd')
    for iCol in range(M):
        Amod = A - eigs[iCol] * numerix.identity(M)[..., numerix.newaxis]
        tmp = Amod.transpose(2, 0, 1)
        Amat = scipy.sparse.bsr_matrix((tmp, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
        X = numerix.zeros((N, M), 'd')
        X[:,iCol] = 1.
        vectors[:,:,iCol] = scipy.sparse.linalg.gmres(Amat.asformat("csc"), numerix.zeros(M * N), X.ravel())[0].reshape((N, M))
        vectors[:,:,iCol] = vectors[:,:,iCol] / numerix.sqrt(numerix.sum(vectors[:,:,iCol]**2, 1))[...,numerix.newaxis]
    return vectors.transpose(1, 2, 0)

def eigvec2(A):
    """
    Calculates the eigenvalues and eigenvectors of N, 2x2 matrices, A.shape = (2, 2, N).
    """
    tr = A[0,0] + A[1,1]
    det = A[0,0] * A[1,1] - A[1,0] * A[0,1]
    DeltaRt = numerix.sqrt(tr**2 - 4 * det)
    eigs = numerix.array(((tr + DeltaRt) / 2., (tr - DeltaRt) / 2.))
    return eigs, vecN(A, eigs)

def eigvec(A):
    M, M, N = A.shape
    if M == 2:
        return eigvec2(A)
    else:
        return eigvecN(A)

def sortedeig(A):
    """
    Caclulates the sorted eigenvalues and eigenvectors of N, MxM matrices, A.shape = (M, M, N).
    """
    N = A.shape[-1]
    
    eigenvalues, R = eigvec(A)

    order = eigenvalues.argsort(0).swapaxes(0, 1)
    Nlist = [[i] for i in xrange(N)]
    return (eigenvalues[order, Nlist].swapaxes(0, 1),
            R[:, order, Nlist].swapaxes(1, 2))

def matvec(A, v):
    """
    Multiply N, MxM matrices by N, M length vectors, A.shape = (M, M, N), v.shape = (M, N).
    """
    return numerix.sum(A * v, 1)

def dot(v0, v1):
    """
    Dot product of two N, M-length vectors, v0.shape = (M, N), v1.shape = (M, N).
    """
    return numerix.sum(v0 * v1, 0)


