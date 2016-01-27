#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "pysparseMatrix.py"
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

__docformat__ = 'restructuredtext'

__all__ = []

from pysparse import spmatrix
from fipy.tools import numerix

from fipy.matrices.sparseMatrix import _SparseMatrix

class _PysparseMatrix(_SparseMatrix):

    def __init__(self, matrix):
        """Creates a wrapper for a pysparse matrix

        Allows basic python operations __add__, __sub__ etc.
        Facilitate matrix populating in an easy way.

        :Parameters:
          - `matrix`: The starting `spmatrix`
        """
        self.matrix = matrix

    def getCoupledClass(self):
        return _CoupledPysparseMeshMatrix

    def copy(self):
        return _PysparseMatrix(matrix=self.matrix.copy())

    def __getitem__(self, index):
        m = self.matrix[index]
        if type(m) is type(0) or type(m) is type(0.):
            return m
        else:
            return _PysparseMatrix(matrix=m)

    def __iadd__(self, other):
        self._iadd(self.matrix, other)
        return self

    def _iadd(self, L, other, sign = 1):
        if other != 0:
            L.shift(sign, other.matrix)

    def __add__(self, other):
        """
        Add two sparse matrices

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L + _PysparseIdentityMatrix(size=3)
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
            L = self.matrix.copy()
            L.shift(1, other.matrix)
            return _PysparseMatrix(matrix=L)

    __radd__ = __add__

    def __sub__(self, other):

        if other == 0:
            return self
        else:
            L = self.matrix.copy()
            L.shift(-1, other.matrix)
            return _PysparseMatrix(matrix=L)

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
            return self._iadd(self.matrix, other, -1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix

            >>> L1 = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PysparseIdentityMatrix(size=3)
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
        N = self.matrix.shape[1]

        if isinstance(other, _PysparseMatrix):
            return _PysparseMatrix(matrix=spmatrix.matrixmultiply(self.matrix, other.matrix))
        else:
            shape = numerix.shape(other)
            if shape == ():
                L = spmatrix.ll_mat(N, N, N)
                L.put(other * numerix.ones(N, 'l'))
                return _PysparseMatrix(matrix=spmatrix.matrixmultiply(self.matrix, L))
            elif shape == (N,):
                y = numerix.empty((self.matrix.shape[0],))
                self.matrix.matvec(other, y)
                return y
            else:
                raise TypeError

    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            y = other.copy()
            self.matrix.matvec_transp(other, y)
            return y
        else:
            return self * other

    @property
    def _shape(self):
        return self.matrix.shape

    @property
    def _range(self):
        return range(self._shape[1]), range(self._shape[0])

    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        self.matrix.put(vector, id1, id2)

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
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
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._shape[0])
            tmp = numerix.zeros((self._shape[0],), 'd')
            tmp[:] = vector
            self.put(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.put(vector, ids, ids)

    def take(self, id1, id2):
        vector = numerix.zeros(len(id1), 'd')
        self.matrix.take(vector, id1, id2)
        return vector

    def takeDiagonal(self):
        ids = numerix.arange(self._shape[0])
        return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
            >>> print L
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        self.matrix.update_add_at(vector, id1, id2)

    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._shape[0])
            tmp = numerix.zeros((self._shape[0],), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)

    @property
    def numpyArray(self):
        shape = self._shape
        indices = numerix.indices(shape)
        numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
        return numerix.reshape(numMatrix, shape)

    def matvec(self, x):
        """
        This method is required for scipy solvers.
        """
        return self * x

    def exportMmf(self, filename):
        """
        Exports the matrix to a Matrix Market file of the given filename.
        """
        self.matrix.export_mtx(filename)

class _PysparseMatrixFromShape(_PysparseMatrix):

    def __init__(self, rows, cols, bandwidth=0, sizeHint=None, matrix=None, storeZeros=True):
        """Instantiates and wraps a pysparse `ll_mat` matrix

        :Parameters:
          - `rows`: The number of matrix rows
          - `cols`: The number of matrix columns
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `storeZeros`: Instructs pysparse to store zero values if possible.

        """
        sizeHint = sizeHint or max(rows, cols) * bandwidth
        if matrix is None:
            tmpMatrix = spmatrix.ll_mat(1, 1, 1)
            if hasattr(tmpMatrix, 'storeZeros'):
                matrix = spmatrix.ll_mat(rows, cols, sizeHint, storeZeros)
            else:
                matrix = spmatrix.ll_mat(rows, cols, sizeHint)

        _PysparseMatrix.__init__(self, matrix=matrix)

class _PysparseMeshMatrix(_PysparseMatrixFromShape):
    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, numberOfVariables=1, numberOfEquations=1, storeZeros=True):

        """Creates a `_PysparseMatrixFromShape` associated with a `Mesh`. Allows for different number of equations and/or variables

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `numberOfVariables`: The columns of the matrix is determined by numberOfVariables * self.mesh.numberOfCells.
          - `numberOfEquations`: The rows of the matrix is determined by numberOfEquations * self.mesh.numberOfCells.
          - `storeZeros`: Instructs pysparse to store zero values if possible.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations
        rows = numberOfEquations * self.mesh.numberOfCells
        cols = numberOfVariables * self.mesh.numberOfCells
        _PysparseMatrixFromShape.__init__(self, rows=rows, cols=cols, bandwidth=bandwidth, sizeHint=sizeHint, matrix=matrix, storeZeros=storeZeros)

    def __mul__(self, other):
        if isinstance(other, _PysparseMeshMatrix):
            return _PysparseMeshMatrix(mesh=self.mesh,
                                       matrix=spmatrix.matrixmultiply(self.matrix, other.matrix))
        else:
            return _PysparseMatrixFromShape.__mul__(self, other)

    def asTrilinosMeshMatrix(self):
        """Transforms a pysparse matrix into a trilinos matrix and maintains the
        trilinos matrix as an attribute.

        :Returns:
          The trilinos matrix.

        """
        A = self.matrix.copy()
        values, irow, jcol = A.find()

        if not hasattr(self, 'trilinosMatrix'):
            if A.shape[0] == 0:
                bandwidth = 0
            else:
                bandwidth = int(numerix.ceil(float(len(values)) / float(A.shape[0])))
            bandwidth = 1
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrixKeepStencil
            self.trilinosMatrix = _TrilinosMeshMatrixKeepStencil(mesh=self.mesh, bandwidth=bandwidth,
                                                                 numberOfVariables=self.numberOfVariables,
                                                                 numberOfEquations=self.numberOfEquations)

        self.trilinosMatrix.addAt(values, irow, jcol)
        self.trilinosMatrix.finalize()

        return self.trilinosMatrix

    @property
    def numpyArray(self):
        from fipy.tools import parallelComm
        if parallelComm.Nproc == 1:
            return super(_PysparseMeshMatrix, self).numpyArray
        else:
            return self.asTrilinosMeshMatrix().numpyArray


    def flush(self):
        """
        Deletes the copy of the pysparse matrix held and calls `self.trilinosMatrix.flush()` if necessary.
        """

        if hasattr(self, 'trilinosMatrix'):
            if hasattr(self.matrix, 'storeZeros'):
                self.trilinosMatrix.flush(cacheStencil=self.matrix.storeZeros)
            else:
                self.trilinosMatrix.flush(cacheStencil=False)

        if (not hasattr(self, 'cache')) or (self.cache is False):
            del self.matrix

    def _test(self):
        """
        Tests

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3, storeZeros=True)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(m.matrix.keys(), [(0, 1), (1, 0), (2, 2)])
        True
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(m.matrix.values(), [1., 2., 0.])
        True
        >>> m = _PysparseMatrixFromShape(rows=3, cols=3, storeZeros=False)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> print numerix.allequal(m.matrix.keys(), [(0, 1), (1, 0)])
        True
        >>> print numerix.allequal(m.matrix.values(), numerix.array([1.0, 2.0]))
        True

        """
        pass

class _PysparseIdentityMatrix(_PysparseMatrixFromShape):
    """
    Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size):
        """Create a sparse matrix with '1' in the diagonal

            >>> print _PysparseIdentityMatrix(size=3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PysparseMatrixFromShape.__init__(self, rows=size, cols=size, bandwidth = 1)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)

class _PysparseIdentityMeshMatrix(_PysparseIdentityMatrix):
    def __init__(self, mesh):
        """Create a sparse matrix associated with a `Mesh` with '1' in the diagonal

            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print _PysparseIdentityMeshMatrix(mesh=mesh)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PysparseIdentityMatrix.__init__(self, size=mesh.numberOfCells)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
