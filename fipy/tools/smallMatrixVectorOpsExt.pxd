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

