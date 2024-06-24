from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

__all__ = []

import scipy.sparse as sp
from fipy.tools import numerix

from fipy.matrices.sparseMatrix import _SparseMatrix

class _ScipyMatrix(_SparseMatrix):

    """class wrapper for a scipy sparse matrix.

    `_ScipyMatrix` is always `NxN`.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, matrix):
        """Creates a `_ScipyMatrix`.

        Parameters
        ----------
        matrix : ~scipy.sparse.csr_matrix
            The internal SciPy matrix
        """
        self.matrix = matrix

        super(_ScipyMatrix, self).__init__()

    def copy(self):
        return _ScipyMatrix(matrix=self.matrix.copy())

    def __getitem__(self, index):
        m = self.matrix[index]
        if isinstance(m, (type(0), type(0.))):
            return m
        else:
            return _ScipyMatrix(matrix=m)

    def __iadd__(self, other):
        return self._iadd(other)

    def _iadd(self, other, sign=1):
        if hasattr(other, "matrix"):
            self.matrix = self.matrix + (sign * other.matrix)
        elif isinstance(other, (float, int)):
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

            >>> L = _ScipyMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> print(L + _ScipyIdentityMatrix(size=3))
             1.000000  10.000000   3.000000  
                ---     4.141593      ---    
             2.500000      ---     1.000000  

            >>> print(L + 0)
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    

            >>> print(L + 3)
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

        >>> L1 = _ScipyMatrixFromShape(rows=3, cols=3)
        >>> L1.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
        >>> L2 = _ScipyIdentityMatrix(size=3)
        >>> L2.put([4.38, 12357.2, 1.1], [2, 1, 0], [1, 0, 2])

        >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
        ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
        ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

        >>> numerix.allclose((L1 * L2).numpyArray, tmp)
        1

        or a sparse matrix by a vector

        >>> tmp = numerix.array((29., 6.28318531, 2.5))
        >>> numerix.allclose(L1 * numerix.array((1, 2, 3), 'd'), tmp)
        1

        or a vector by a sparse matrix

        >>> tmp = numerix.array((7.5, 16.28318531,  3.))
        >>> numerix.allclose(numerix.array((1, 2, 3), 'd') * L1, tmp)
        1

        (The multiplication is broken.  Numpy is calling __rmul__ for every
        element instead of with the whole array.)
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
        if isinstance(numerix.ones(1, 'l'), type(other)):
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
        return list(range(self._shape[1])), list(range(self._shape[0]))

    def put(self, vector, id1, id2, overlapping=False):
        """Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)

        Parameters
        ----------
        vector : array_like
            The values to insert.
        id1 : array_like
            The row indices.
        id2 : array_like
            The column indices.
        overlapping : bool
            Whether to insert ghosted values or not (Ignored. Default False).

        Examples
        --------

            >>> L = _ScipyMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> print(L)
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        assert len(id1) == len(id2) == len(vector)

        # done in such a way to vectorize everything
        tempVec = numerix.array(vector) - self.matrix[id1, id2].flat
        tempMat = sp.csr_matrix((tempVec, (id1, id2)), self.matrix.shape)

        self.matrix = self.matrix + tempMat

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix

            >>> L = _ScipyMatrixFromShape(rows=3, cols=3)
            >>> L.putDiagonal([3., 10., numerix.pi])
            >>> print(L)
             3.000000      ---        ---    
                ---    10.000000      ---    
                ---        ---     3.141593  
            >>> L.putDiagonal([10., 3.])
            >>> print(L)
            10.000000      ---        ---    
                ---     3.000000      ---    
                ---        ---     3.141593  
        """
        if isinstance(vector, (int, float)):
            vector = numerix.repeat(vector, self._shape[0])

        self.matrix.setdiag(vector)

    def take(self, id1, id2):
        return self.matrix[id1, id2]

    def takeDiagonal(self):
        return self.matrix.diagonal()

    def addAt(self, vector, id1, id2):
        """Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)

        Parameters
        ----------
        vector : array_like
            The values to insert.
        id1 : array_like
            The row indices.
        id2 : array_like
            The column indices.
        overlapping : bool
            Whether to add ghosted values or not (Ignored. Default False).

        Examples
        --------

            >>> L = _ScipyMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> L.addAt([1.73, 2.2, 8.4, 3.9, 1.23], [1, 2, 0, 0, 1], [2, 2, 0, 0, 2])
            >>> print(L)
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        assert len(id1) == len(id2) == len(vector)

        temp = sp.csr_matrix((vector, (id1, id2)), self.matrix.shape)

        self.matrix = self.matrix + temp

    def addAtDiagonal(self, vector):
        if isinstance(vector, (int, float)):
            vector = numerix.repeat(vector, self._shape[0])

        ids = numerix.arange(len(vector))
        self.addAt(vector, ids, ids)

    @property
    def numpyArray(self):
        return self.matrix.toarray()

    def matvec(self, x):
        """This method is required for scipy solvers.
        """
        return self * x

    def exportMmf(self, filename):
        """Exports the matrix to a Matrix Market file of the given `filename`.
        """
        from scipy.io import mmio

        mmio.mmwrite(filename, self.matrix)

    @property
    def CSR(self):
        """The Compact Sparse Row description of the matrix

        Returns
        -------
        ptrs : array_like of int
            Locations in `cols` and `data` vectors that start a row,
            terminated with len(data) + 1
        cols : array_like of int
            Sequence of non-sparse column indices.
        data : array_like of float
            Sequence of non-sparse values.

        Examples
        --------

        >>> L = _ScipyMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
        >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
        >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
        >>> ptrs, cols, data = L.CSR
        >>> print(numerix.asarray(ptrs))
        [0 3 5 7]
        >>> print(numerix.asarray(cols))
        [0 1 2 1 2 0 2]
        >>> print(numerix.asarray(data))
        [ 12.3         10.           3.           3.14159265   2.96
           2.5          2.2       ]
        """
        return self.matrix.indptr, self.matrix.indices, self.matrix.data

    @property
    def LIL(self):
        """The List of Lists description of the local matrix

        Returns
        -------
        rows : list of sequence of int
            List of non-sparse column indices on each row.
        data : list of sequence of float
            List of non-sparse values on each row.

        Examples
        --------

        >>> L = _ScipyMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
        >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
        >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
        >>> rows, data = L.LIL
        >>> from scipy.stats.mstats import argstoarray # doctest: +SCIPY
        >>> print(argstoarray(*rows)) # doctest: +SCIPY
        [[0.0 1.0 2.0]
         [1.0 2.0 --]
         [0.0 2.0 --]]
        >>> print(argstoarray(*data)) # doctest: +SCIPY
        [[12.3 10.0 3.0]
         [3.141592653589793 2.96 --]
         [2.5 2.2 --]]
        """
        lil = self.matrix.tolil()

        return lil.rows, lil.data

    @property
    def T(self):
        """Transpose matrix

        Returns
        -------
        ~fipy.matrices.scipyMatrix._ScipyMatrix

        Examples
        --------

        >>> import fipy as fp

        >>> mesh = fp.Grid1D(nx=10)
        >>> ids = fp.CellVariable(mesh=mesh, value=mesh._globalOverlappingCellIDs)

        >>> mat = _ScipyColMeshMatrix(mesh=mesh, rows=1)
        >>> mat.put(vector=ids.value,
        ...         id1=[fp.parallelComm.procID] * mesh.numberOfCells,
        ...         id2=mesh._localOverlappingCellIDs,
        ...         overlapping=True)

        >>> print(mat.T.numpyArray)
        [[ 0.]
         [ 1.]
         [ 2.]
         [ 3.]
         [ 4.]
         [ 5.]
         [ 6.]
         [ 7.]
         [ 8.]
         [ 9.]]
        """
        return _ScipyMatrix(matrix=self.matrix.transpose(copy=True))

class _ScipyMatrixFromShape(_ScipyMatrix):

    def __init__(self, rows, cols,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Instantiates and wraps a scipy sparse matrix

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        nonZerosPerRow : int or array_like of int
            *ignored*
        exactNonZeros : bool
            *ignored*
        matrix : ~scipy.sparse.csr_matrix
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        if matrix is None:
            matrix = sp.csr_matrix((rows, cols))

        super(_ScipyMatrixFromShape, self).__init__(matrix=matrix)

class _ScipyBaseMeshMatrix(_ScipyMatrixFromShape):
    def __init__(self, mesh, rows, cols,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_ScipyMatrixFromShape` associated with a `Mesh`.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        rows : int
            The number of local matrix rows.
        cols : int
            The number of local matrix columns.
        nonZerosPerRow : int or array_like of int
            *ignored*
        exactNonZeros : bool
            *ignored*
        matrix : ~scipy.sparse.csr_matrix
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        self.mesh = mesh

        super(_ScipyBaseMeshMatrix, self).__init__(rows=rows,
                                                   cols=cols,
                                                   nonZerosPerRow=nonZerosPerRow,
                                                   exactNonZeros=exactNonZeros,
                                                   matrix=matrix,
                                                   storeZeros=storeZeros)

    def _getGhostedValues(self, var):
        """Obtain current ghost values from across processes

        Nothing to do for serial matrix.

        Returns
        -------
        ndarray
            Ghosted values
        """
        return var.value

class _ScipyRowMeshMatrix(_ScipyBaseMeshMatrix):
    def __init__(self, mesh, cols, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_ScipyBaseMeshMatrix` with rows associated with equations.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        cols : int
            The number of matrix columns.
        numberOfEquations : int
            The rows of the matrix are determined by
            `numberOfEquations * mesh.numberOfCells`.
        nonZerosPerRow : int or array_like of int
            *ignored*
        exactNonZeros : bool
            *ignored*
        matrix : ~scipy.sparse.csr_matrix
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        self.numberOfEquations = numberOfEquations

        super(_ScipyRowMeshMatrix, self).__init__(mesh=mesh,
                                                  rows=numberOfEquations * mesh.numberOfCells,
                                                  cols=cols,
                                                  nonZerosPerRow=nonZerosPerRow,
                                                  exactNonZeros=exactNonZeros,
                                                  matrix=matrix,
                                                  storeZeros=storeZeros)

class _ScipyColMeshMatrix(_ScipyBaseMeshMatrix):
    def __init__(self, mesh, rows, numberOfVariables=1,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_ScipyBaseMeshMatrix` with columns associated with solution variables.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        rows : int
            The number of matrix rows.
        numberOfVariables : int
            The columns of the matrix are determined by
            `numberOfVariables * mesh.globalNumberOfCells`.
        nonZerosPerRow : int or array_like of int
            *ignored*
        exactNonZeros : bool
            *ignored*
        matrix : ~scipy.sparse.csr_matrix
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        self.numberOfVariables = numberOfVariables

        super(_ScipyColMeshMatrix, self).__init__(mesh=mesh,
                                                  rows=rows,
                                                  cols=numberOfVariables * mesh.numberOfCells,
                                                  nonZerosPerRow=nonZerosPerRow,
                                                  exactNonZeros=exactNonZeros,
                                                  matrix=matrix,
                                                  storeZeros=storeZeros)

class _ScipyMeshMatrix(_ScipyRowMeshMatrix):
    def __init__(self, mesh, numberOfVariables=1, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False, matrix=None, storeZeros=True):
        """Creates a `_ScipyBaseMeshMatrix` associated with equations and variables.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        numberOfVariables : int
            The columns of the matrix are determined by
            `numberOfVariables * mesh.numberOfCells`.
        numberOfEquations : int
            The rows of the matrix are determined by
            `numberOfEquations * mesh.numberOfCells`.
        nonZerosPerRow : int or array_like of int
            *ignored*
        exactNonZeros : bool
            *ignored*
        matrix : ~scipy.sparse.csr_matrix
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        self.numberOfVariables = numberOfVariables

        super(_ScipyMeshMatrix, self).__init__(mesh=mesh,
                                               cols=numberOfVariables * mesh.numberOfCells,
                                               nonZerosPerRow=nonZerosPerRow,
                                               exactNonZeros=exactNonZeros,
                                               matrix=matrix,
                                               numberOfEquations=numberOfEquations,
                                               storeZeros=storeZeros)

    def __mul__(self, other):
        if isinstance(other, _ScipyMeshMatrix):
            return _ScipyMeshMatrix(mesh=self.mesh,
                                    matrix=(self.matrix * other.matrix))
        else:
            return _ScipyMatrixFromShape.__mul__(self, other)

    def asTrilinosMeshMatrix(self):
        """Transforms a scipy matrix into a trilinos matrix and maintains the
        trilinos matrix as an attribute.

        Returns
        -------
        ~fipy.matrices.trilinosMatrix._TrilinosMatrix
        """
        raise NotImplementedError

    def flush(self):
        """Deletes the scipy matrix and calls `self.trilinosMatrix.flush()` if necessary.
        """
        if not getattr(self, 'cache', False):
            del self.matrix

    def _test(self):
        """
        Tests

        >>> m = _ScipyMatrixFromShape(rows=3, cols=3, storeZeros=True)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> nonZeroIdx = m.matrix.nonzero()
        >>> print(not hasattr(m.matrix, 'storeZeros')
        ...       or numerix.allequal(nonZeroIdx, [(0, 1), (1, 0), (2, 2)]))
        True
        >>> print(not hasattr(m.matrix, 'storeZeros')
        ...       or numerix.allequal(m.matrix[nonZeroIdx].toarray(),
        ...                           [1., 2., 0.]))
        True
        >>> m = _ScipyMatrixFromShape(rows=3, cols=3, storeZeros=False)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> nonZeroIdx = m.matrix.nonzero()
        >>> print(numerix.allequal(nonZeroIdx, [(0, 1), (1, 0)]))
        True
        >>> print(numerix.allequal(numerix.array(m.matrix[nonZeroIdx]), numerix.array([1.0, 2.0])))
        True

        """
        pass

class _ScipyIdentityMatrix(_ScipyMatrixFromShape):
    """
    Represents a sparse identity matrix for scipy.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with `1` in the diagonal

            >>> print(_ScipyIdentityMatrix(size=3))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _ScipyMatrixFromShape.__init__(self, rows=size, cols=size, nonZerosPerRow=1)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)

class _ScipyIdentityMeshMatrix(_ScipyIdentityMatrix):
    def __init__(self, mesh):
        """
        Create a sparse matrix associated with a `Mesh` with `1` in the diagonal

            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print(_ScipyIdentityMeshMatrix(mesh=mesh))
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
