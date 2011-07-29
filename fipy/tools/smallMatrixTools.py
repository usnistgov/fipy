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

def mulNew(A, B):
    """
    Matrix multiply N, MxM matrices, A.shape = B.shape = (N, M, M).
    """
    return numerix.sum(A.swapaxes(1,2)[:, :, :, numerix.newaxis] * B[:, :, numerix.newaxis], 1)
    
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

def invmulNew(A, B):
    """
    Calculates the product (A^-1.B) of N, M (rows) x P (cols) matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (N, M, N), B.shape = (N, M, M).
    """
##    tmp = A.transpose(2, 0, 1)
    N, M, P = B.shape
    N, M, M = A.shape
    
    import scipy
    Amat = scipy.sparse.bsr_matrix((A, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
    ans = numerix.zeros((N, M, P), 'd')
##    Bvec = B.transpose(2, 0, 1)

    LU = scipy.sparse.linalg.splu(Amat.asformat("csc"))
    for iCol in range(P):
        import scipy.sparse.linalg
        ans[:,:,iCol] = LU.solve(B[:,:,iCol].ravel()).reshape((N, M))

    return ans

def mulinv(A, B):
    """
    Calculates the product (A.B^-1) of N, MxM matrices (A) and the inverse of N, MxM matrices (B).
    A.shape = (M, M, N), B.shape = (M, M, N).
    """
    return invmul(B.transpose(1, 0, 2), A.transpose(1, 0, 2)).transpose(1, 0, 2)

def mulinvNew(A, B):
    """
    Calculates the product (A.B^-1) of N, MxM matrices (A) and the inverse of N, MxM matrices (B).
    A.shape = (N, M, M), B.shape = (N, M, M).
    """
    return invmulNew(B.transpose(0, 2, 1), A.transpose(0, 2, 1)).transpose(0, 2, 1)

def invmatvec(A, v):
    """
    Calculates the product (A^-1.v) of N, MxM matrices (B) and the inverse of N, MxM matrices (A).
    A.shape = (M, M, N), v.shape = (M, N).
    """
    return invmul(A, v[:,numerix.newaxis,:]).reshape(v.shape)
    
def eigvecN(A):
    """
    Calculate the eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    tmp = zip(*map(numerix.linalg.eig, A))
    return numerix.array(tmp[0]), numerix.array(tmp[1])

def matvec(A, v):
    """
    Multiply N, MxM matrices by N, M length vectors, A.shape = (M, M, N), v.shape = (M, N).
    """
    return numerix.sum(A * v, 1)

def vecN(A, eigs):
    """
    Calculates the eigenvectors for N, MxM matrices (A) using the eigenvalues (eigs).
    A.shape = (N, M, M), eigs.shape = (N, M).
    """
    N, M, M = A.shape
    vectors = numerix.zeros((N, M, M), 'd')
    for iCol in range(M):
        Amod = A - eigs[:,iCol,numerix.newaxis,numerix.newaxis] * numerix.identity(M)
        Amat = scipy.sparse.bsr_matrix((Amod, numerix.arange(N), numerix.arange(N + 1)), shape=(M * N, M * N))
        X = numerix.zeros((N, M), 'd')
        X[:,iCol] = 1.
        vectors[:,:,iCol] = scipy.sparse.linalg.gmres(Amat.asformat("csc"), numerix.zeros(M * N), X.ravel())[0].reshape((N, M))
        vectors[:,:,iCol] = vectors[:,:,iCol] / numerix.sqrt(numerix.sum(vectors[:,:,iCol]**2, 1))[...,numerix.newaxis]
    return vectors

def eigvec2(A):
    """
    Calculates the eigenvalues and eigenvectors of N, 2x2 matrices, A.shape = (N, 2, 2).
    """
    tr = A[:,0,0] + A[:,1,1]
    det = A[:,0,0] * A[:,1,1] - A[:,1,0] * A[:,0,1]
    DeltaRt = numerix.sqrt(tr**2 - 4 * det)
    eigs = numerix.array(((tr + DeltaRt) / 2., (tr - DeltaRt) / 2.)).transpose(1, 0)
    return eigs, vecN(A, eigs)

def eigvec(A):
    N, M, M = A.shape
    if M == 2:
        return eigvec2(A)
    else:
        return eigvecN(A)

def sortedeig(A):
    """
    Caclulates the sorted eigenvalues and eigenvectors of N, MxM matrices, A.shape = (N, M, M).
    """
    N = A.shape[0]
    eigenvalues, R = eigvec(A)
    order = eigenvalues.argsort(1)
    Nlist = [[i] for i in xrange(N)]
    return eigenvalues[Nlist, order], R[Nlist, :, order].transpose(0, 2, 1)


def dot(v0, v1):
    """
    Dot product of two N, M-length vectors, v0.shape = (M, N), v1.shape = (M, N).
    """
    return numerix.sum(v0 * v1, 0)


