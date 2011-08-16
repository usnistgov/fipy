import numpy as np
import pyximport
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension('smt', ['smt.pyx'],
                         libraries=['lapack'])]

pyximport.install(setup_args = {'cmdclass' : {'build_ext' : build_ext},
                                'ext_modules' : ext_modules})
                             
from smt import solve, fasteigvec

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
