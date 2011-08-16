#cython: profile=True
import cython
import numpy as np
from numpy.core import intc	
cimport numpy as np
from numpy.compat import asbytes

@cython.boundscheck(False)
cdef np.ndarray[complex, ndim=3] zsolve(np.ndarray[complex, ndim=3] A, np.ndarray[complex, ndim=3] B):
    cdef int N = A.shape[0]
    cdef int M = A.shape[1]
    cdef unsigned int i = 0
    cdef unsigned int j = 0
    cdef unsigned int k = 0
    cdef np.ndarray[complex, ndim=2] Bi = np.zeros((M, M), dtype=complex)
    cdef np.ndarray[complex, ndim=2] Ai = np.zeros((M, M), dtype=complex)
    cdef int NRHS = B.shape[1]
    cdef int info = 0     
    cdef np.ndarray ipiv = np.zeros(M, intc)

    for i from 0 <= i < N:

        for j from 0 <= j < M:
            for k from 0 <= k < M:
                Ai[j, k] = A[i, j, k]
                Bi[j, k] = B[i, j, k]

        zgesv_(&M, &NRHS, <complex *> Ai.data, &M, <int *> ipiv.data, <complex *> Bi.data, &M, &info)

        for j from 0 <= j < M:
            for k from 0 <= k < M:
                B[i, j, k] = Bi[j, k] 

    return B

@cython.boundscheck(False)
cdef np.ndarray[double, ndim=3] dsolve(np.ndarray[double, ndim=3] A, np.ndarray[double, ndim=3] B):
    cdef int N = A.shape[0]
    cdef int M = A.shape[1]
    cdef unsigned int i = 0
    cdef unsigned int j = 0
    cdef unsigned int k = 0
    cdef np.ndarray[double, ndim=2] Bi = np.zeros((M, M), 'd')
    cdef np.ndarray[double, ndim=2] Ai = np.zeros((M, M), 'd')
    cdef int NRHS = B.shape[1]
    cdef int info = 0     
    cdef np.ndarray ipiv = np.zeros(M, intc)

    for i from 0 <= i < N:
        
        for j from 0 <= j < M:
            for k from 0 <= k < M:
                Ai[j, k] = A[i, j, k]
                Bi[j, k] = B[i, j, k]

        dgesv_(&M, &NRHS, <double *> Ai.data, &M, <int *> ipiv.data, <double *> Bi.data, &M, &info)

        for j from 0 <= j < M:
            for k from 0 <= k < M:
                B[i, j, k] = Bi[j, k] 

    return B

def solve(A, B):
    A = A.transpose(0, 2, 1).copy()
    B = B.transpose(0, 2, 1).copy()

    if np.iscomplexobj(A) or np.iscomplexobj(B):
        A.dtype = complex
        B.dtype = complex
        return zsolve(A, B)
    else:
        return dsolve(A, B)

@cython.boundscheck(False)
def fasteigvec(A):
    cdef int N
    cdef int M
    N, M, M = A.shape
    cdef np.ndarray[double, ndim=3] Anew = A.transpose(0, 2, 1).copy()

    cdef np.ndarray[double, ndim=2] Ai = np.zeros((M, M))
    cdef np.ndarray[double, ndim=1] eigsi = np.zeros((M,))
    cdef np.ndarray[double, ndim=1] Ieigsi = np.zeros((M,))
    cdef np.ndarray[double, ndim=2] vri = np.zeros((M, M))

    cdef np.ndarray[double, ndim=2] eigs = np.zeros((N, M))
    cdef np.ndarray[double, ndim=2] Ieigs = np.zeros((N, M))
    cdef np.ndarray[double, ndim=3] vr = np.zeros((N, M, M))
    cdef np.ndarray[double, ndim=1] dummy = np.zeros((1,), 'd')
    cdef np.ndarray[complex, ndim=3] zvecs

    cdef np.ndarray[int, ndim=1] indi
    cdef np.ndarray[int, ndim=1] indj

    cdef unsigned int i = 0
    cdef unsigned int j = 0
    cdef unsigned int k = 0
    cdef unsigned int iInd = 0
    cdef unsigned int jInd = 0
    cdef unsigned int jjInd = 0

    cdef int lwork = 1
    cdef np.ndarray[double, ndim=1] work = np.zeros((lwork,), 'd')
    cdef int vl = 1 
    cdef int info = 0  
    cdef int len_jobvl = 1
    cdef int len_jobvr = 1
    cdef int minusone = -1
    cdef char _N = 'N'
    cdef char _V = 'V'

    i = 0
    for j from 0 <= j < M:
        for k from 0 <= k < M:
            Ai[j, k] = Anew[i, j, k]

    ##Ai = Anew[0]
    ##vri = vr[0]
    ##eigsi = eigs[0]
    ##Ieigsi = Ieigs[0]
    dgeev_(&_N, &_V, &M, <double *> Ai.data,
           &M, <double *> eigsi.data, <double *> Ieigsi.data, <double *> dummy.data, &vl,
           <double *> vri.data, &M, <double *> work.data, &minusone, &info,
           len_jobvl, len_jobvr)
    
    lwork = int(abs(work[0]))
    work = np.zeros((lwork,), 'd')

    for i from 0 <= i < N:
        for j from 0 <= j < M:
            for k from 0 <= k < M:
                Ai[j, k] = Anew[i, j, k]
                
        dgeev_(&_N, &_V, &M, <double *> Ai.data,
               &M, <double *> eigsi.data, <double *> Ieigsi.data, <double *> dummy.data, &vl,
               <double *> vri.data, &M, <double *> work.data, &lwork, &info,
               len_jobvl, len_jobvr)

        for j from 0 <= j < M:
            eigs[i, j] = eigsi[j]
            Ieigs[i, j] = Ieigsi[j] 
            for k from 0 <= k < M:
                vr[i, j, k] = vri[j, k]

    if (Ieigs.flatten() != 0).any():
        zeigs = eigs + 1j * Ieigs
        zvecs = np.array(vr, dtype=complex)
        indi, indj = np.nonzero(Ieigs != 0)
        i = 0
        while i < (len(indi) - 1):
            iInd = indi[i]
            jInd = indj[i]
            jjInd = indj[i + 1]
            if iInd == indi[i + 1]:
                for j from 0 <= j < M:
                    zvecs[iInd, jInd, j] = vr[iInd, jInd, j] + 1j * vr[iInd, jjInd, j]
                    zvecs[iInd, jjInd, j] = vr[iInd, jInd, j] - 1j * vr[iInd, jjInd, j]
                i += 2
            else:
                i += 1
        return zeigs, zvecs.transpose(0, 2, 1)
    else:
        return eigs, vr.transpose(0, 2, 1)

