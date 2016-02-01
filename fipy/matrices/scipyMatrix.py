#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "scipyMatrix.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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

__docformat__ = 'restructuredtext'

__all__ = []

import scipy.sparse as sp
from fipy.tools import numerix

from fipy.matrices.sparseMatrix import _SparseMatrix

class _ScipyMatrix(_SparseMatrix):

    """class wrapper for a scipy sparse matrix.

    `_ScipyMatrix` is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, matrix):
        """Creates a `_ScipyMatrix`.

        :Parameters:
          - `matrix`: The starting `spmatrix`
        """
        self.matrix = matrix

    def getCoupledClass(self):
        return _CoupledScipyMeshMatrix

    def copy(self):
        return _ScipyMatrix(matrix=self.matrix.copy())

    def __getitem__(self, index):
        m = self.matrix[index]
        if type(m) is type(0) or type(m) is type(0.):
            return m
        else:
            return _ScipyMatrix(matrix=m)

    def __iadd__(self, other):
        return self._iadd(other)

    def _iadd(self, other, sign=1):
        if hasattr(other, "matrix"):
            self.matrix = self.matrix + (sign * other.matrix)
        elif type(other) in [float, int]:
            fillVec = numerix.repeat(other, self.matrix.nnz)

            self.matrix = self.matrix \
                          + sp.csr_matrix((fillVec, self.matrix.nonzero()),
                                          self.matrix.shape)
        else:
            self.matrix = self.matrix + (sign * other)

        return self

    def __add__(self, other):
        """
        Add two sparse matrices

            >>> L = _ScipyMatrixFromShape(size=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L + _ScipyIdentityMatrix(size=3)
             1.000000  10.000000   3.000000  
                ---     4.141593      ---    
             2.500000      ---     1.000000  

            >>> print L + 0
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    

            >>> print L + 3
            Traceback (most recent call last):
            ...
            AttributeError: 'int' object has no attribute 'matrix'
        """

        if other == 0:
            return self
        else:
            return _ScipyMatrix(matrix=self.matrix + other.matrix)

    __radd__ = __add__

    def __sub__(self, other):

        if other == 0:
            return self
        else:
            L = self.matrix.copy()
            L -= other.matrix
            return _ScipyMatrix(matrix=L)

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
            return self._iadd(other, sign=-1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix

        >>> L1 = _ScipyMatrixFromShape(size=3)
        >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
        >>> L2 = _ScipyIdentityMatrix(size=3)
        >>> L2.put([4.38,12357.2,1.1], [2,1,0], [1,0,2])

        >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
        ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
        ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

        >>> numerix.allclose((L1 * L2).numpyArray, tmp)
        1

        or a sparse matrix by a vector

        >>> tmp = numerix.array((29., 6.28318531, 2.5))
        >>> numerix.allclose(L1 * numerix.array((1,2,3),'d'), tmp)
        1

        or a vector by a sparse matrix

        >>> tmp = numerix.array((7.5, 16.28318531,  3.))
        >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp) ## The multiplication is broken. Numpy is calling __rmul__ for every element instead of with  the whole array.
        1
        """
        N = self.matrix.shape[0]

        if isinstance(other, _ScipyMatrix):
            return _ScipyMatrix(matrix=(self.matrix * other.matrix))
        else:
            shape = numerix.shape(other)
            if shape == ():
                return _ScipyMatrix(matrix=(self.matrix * other))
            elif shape == (N,):
                return self.matrix * other
            else:
                raise TypeError

    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            y = self.matrix.transpose() * other.copy()
            return y
        else:
            return self * other

    def asformat(self, *args, **kwargs):
        return self.matrix.asformat(*args, **kwargs)

    @property
    def _shape(self):
        return self.matrix.shape

    @property
    def _range(self):
        return range(self._shape[1]), range(self._shape[0])

    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)

            >>> L = _ScipyMatrixFromShape(size=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        assert(len(id1) == len(id2) == len(vector))

        # done in such a way to vectorize everything
        tempVec = numerix.array(vector) - self.matrix[id1, id2].flat
        tempMat = sp.csr_matrix((tempVec, (id1, id2)), self.matrix.shape)

        self.matrix = self.matrix + tempMat

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix

            >>> L = _ScipyMatrixFromShape(size=3)
            >>> L.putDiagonal([3.,10.,numerix.pi])
            >>> print L
             3.000000      ---        ---    
                ---    10.000000      ---    
                ---        ---     3.141593  
            >>> L.putDiagonal([10.,3.])
            >>> print L
            10.000000      ---        ---    
                ---     3.000000      ---    
                ---        ---     3.141593  
        """
        if type(vector) in [int, float]:
            vector = numerix.repeat(vector, self._shape[0])

        self.matrix.setdiag(vector)

    def take(self, id1, id2):
        return self.matrix[id1, id2]

    def takeDiagonal(self):
        return self.matrix.diagonal()

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)

            >>> L = _ScipyMatrixFromShape(size=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
            >>> print L
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        assert(len(id1) == len(id2) == len(vector))

        temp = sp.csr_matrix((vector, (id1, id2)), self.matrix.shape)

        self.matrix = self.matrix + temp

    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            vector = numerix.repeat(vector, self._shape[0])

        ids = numerix.arange(len(vector))
        self.addAt(vector, ids, ids)

    @property
    def numpyArray(self):
        return self.matrix.toarray()

    def matvec(self, x):
        """
        This method is required for scipy solvers.
        """
        return self * x

    def __getitem__(self, indices):
        return self.matrix[indices]

class _ScipyMatrixFromShape(_ScipyMatrix):

    def __init__(self, size, bandwidth=0, sizeHint=None, matrix=None, storeZeros=True):
        """Instantiates and wraps a scipy sparse matrix

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `storeZeros`: Instructs scipy to store zero values if possible.

        """
        if matrix is None:
            matrix = sp.csr_matrix((size, size))

        _ScipyMatrix.__init__(self, matrix=matrix)

class _ScipyMeshMatrix(_ScipyMatrixFromShape):

    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, numberOfVariables=1, numberOfEquations=1, storeZeros=True):

        """Creates a `_ScipyMatrixFromShape` associated with a `Mesh`.

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `numberOfVariables`: The columns of the matrix is determined by numberOfVariables * self.mesh.numberOfCells.
          - `numberOfEquations`: The rows of the matrix is determined by numberOfEquations * self.mesh.numberOfCells.
          - `storeZeros`: Instructs scipy to store zero values if possible.


        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        size = self.numberOfVariables * self.mesh.numberOfCells
        assert numberOfEquations == self.numberOfVariables
        _ScipyMatrixFromShape.__init__(self, size=size, matrix=matrix)

    def __mul__(self, other):
        if isinstance(other, _ScipyMeshMatrix):
            return _ScipyMeshMatrix(mesh=self.mesh,
                                    matrix=(self.matrix * other.matrix))
        else:
            return _ScipyMatrixFromShape.__mul__(self, other)

    def asTrilinosMeshMatrix(self):
        """Transforms a scipy matrix into a trilinos matrix and maintains the
        trilinos matrix as an attribute.

        :Returns:
          The trilinos matrix.

        """
        """
        A = self.matrix.copy()
        values, irow, jcol = A.find()

        if not hasattr(self, 'trilinosMatrix'):
            if A.shape[0] == 0:
                bandwidth = 0
            else:
                bandwidth = int(numerix.ceil(float(len(values)) / float(A.shape[0])))
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrixKeepStencil
            self.trilinosMatrix = _TrilinosMeshMatrixKeepStencil(mesh=self.mesh, bandwidth=bandwidth, numberOfVariables=self.numberOfVariables)

        self.trilinosMatrix.addAt(values, irow, jcol)
        self.trilinosMatrix.finalize()

        return self.trilinosMatrix
        """

        raise NotImplementedError

    def flush(self):
        """
        Deletes the copy of the scipy matrix held and calls `self.trilinosMatrix.flush()` if necessary.

        if hasattr(self, 'trilinosMatrix'):
            if hasattr(self.matrix, 'storeZeros'):
                self.trilinosMatrix.flush(cacheStencil=self.matrix.storeZeros)
            else:
                self.trilinosMatrix.flush(cacheStencil=False)
        """

        if (not hasattr(self, 'cache')) or (self.cache is False):
            del self.matrix

    def _test(self):
        """
        Tests

        >>> m = _ScipyMatrixFromShape(size=3, storeZeros=True)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> nonZeroIdx = m.matrix.nonzero()
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(nonZeroIdx, [(0, 1), (1, 0), (2, 2)])
        True
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(m.matrix[nonZeroIdx].toarray(), [1., 2., 0.])
        True
        >>> m = _ScipyMatrixFromShape(size=3, storeZeros=False)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> nonZeroIdx = m.matrix.nonzero()
        >>> print numerix.allequal(nonZeroIdx, [(0, 1), (1, 0)])
        True
        >>> print numerix.allequal(numerix.array(m.matrix[nonZeroIdx]), numerix.array([1.0, 2.0]))
        True

        """
        pass

class _ScipyIdentityMatrix(_ScipyMatrixFromShape):
    """
    Represents a sparse identity matrix for scipy.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal

            >>> print _ScipyIdentityMatrix(size=3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _ScipyMatrixFromShape.__init__(self, size=size, bandwidth = 1)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)

class _ScipyIdentityMeshMatrix(_ScipyIdentityMatrix):
    def __init__(self, mesh):
        """
        Create a sparse matrix associated with a `Mesh` with '1' in the diagonal

            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print _ScipyIdentityMeshMatrix(mesh=mesh)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _ScipyIdentityMatrix.__init__(self, size=mesh.numberOfCells)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
