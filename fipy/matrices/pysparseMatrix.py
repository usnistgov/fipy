from __future__ import unicode_literals
from builtins import range
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

        Parameters
        ----------
        matrix : ~pysparse.spmatrix.ll_mat
            The internal Pysparse matrix
        """
        self.matrix = matrix

        super(_PysparseMatrix, self).__init__()

    def copy(self):
        return _PysparseMatrix(matrix=self.matrix.copy())

    def __getitem__(self, index):
        m = self.matrix[index]
        if not isinstance(m, (type(0), type(0.))):
            m = _PysparseMatrix(matrix=m)
        return m

    @staticmethod
    def _iadd(L, other, sign=1):
        if other != 0:
            L.shift(sign, other.matrix)

    def __iadd__(self, other):
        self._iadd(self.matrix, other)
        return self

    def __add__(self, other):
        """
        Add two sparse matrices

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> print(L + _PysparseIdentityMatrix(size=3))
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
            >>> L1.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> L2 = _PysparseIdentityMatrix(size=3)
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
        if isinstance(numerix.ones(1, 'l'), type(other)):
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
        return list(range(self._shape[1])), list(range(self._shape[0]))

    def put(self, vector, id1, id2):
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

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> print(L)
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        self.matrix.put(vector, id1, id2)

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
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

            >>> L = _PysparseMatrixFromShape(rows=3, cols=3)
            >>> L.put([3., 10., numerix.pi, 2.5], [0, 0, 1, 2], [2, 1, 1, 0])
            >>> L.addAt([1.73, 2.2, 8.4, 3.9, 1.23], [1, 2, 0, 0, 1], [2, 2, 0, 0, 2])
            >>> print(L)
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        self.matrix.update_add_at(vector,
                                  numerix.asarray(id1, dtype='int32'),
                                  numerix.asarray(id2, dtype='int32'))

    def addAtDiagonal(self, vector):
        if isinstance(vector, (int, float)):
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
        Exports the matrix to a Matrix Market file of the given `filename`.
        """
        self.matrix.export_mtx(filename)

    @property
    def CSR(self):
        """The Compact Sparse Row description of the local matrix

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

        >>> L = _PysparseMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
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
        rows, lildata = self.LIL

        ptrs = [0] + [len(row) for row in rows if row]
        ptrs = list(numerix.cumsum(ptrs))
        cols = [col for row in rows for col in row]
        data = [datum for row in lildata for datum in row]

        return ptrs, cols, data

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

        >>> L = _PysparseMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
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
        nrows, _ = self.matrix.shape
        rows = [[] for i in range(nrows)]
        data = [[] for i in range(nrows)]
        for (row, col), datum in self.matrix.items():
            rows[row].append(col)
            data[row].append(datum)

        return rows, data

    @property
    def T(self):
        """Transpose matrix

        Returns
        -------
        ~fipy.matrices.pysparseMatrix._PysparseMatrix

        Examples
        --------

        >>> import fipy as fp

        >>> mesh = fp.Grid1D(nx=10)
        >>> ids = fp.CellVariable(mesh=mesh, value=mesh._globalOverlappingCellIDs)

        >>> mat = _PysparseColMeshMatrix(mesh=mesh, rows=1)
        >>> mat.put(vector=ids.value,
        ...         id1=[fp.parallelComm.procID] * mesh.numberOfCells,
        ...         id2=mesh._localOverlappingCellIDs,
        ...         overlapping=True) # doctest: +SERIAL

        >>> print(mat.T.numpyArray) # doctest: +SERIAL
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
        val, irow, jcol = self.matrix.find()
        rows, cols = self.matrix.shape
        if hasattr(self.matrix, 'storeZeros'):
            A_T = spmatrix.ll_mat(cols, rows, self.matrix.nnz, self.matrix.storeZeros)
        else:
            A_T = spmatrix.ll_mat(cols, rows, self.matrix.nnz)

        A_T.put(val, jcol, irow)

        return _PysparseMatrix(matrix=A_T)

class _PysparseMatrixFromShape(_PysparseMatrix):

    def __init__(self, rows, cols,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Instantiates and wraps a Pysparse `ll_mat` matrix

        Parameters
        ----------
        rows : int
            The number of matrix rows
        cols : int
            The number of matrix columns
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            *ignored*
        matrix : ~pysparse.spmatrix.ll_mat
            Pre-assembled Pysparse matrix to use for storage
        storeZeros : bool
            Instructs pysparse to store zero values if possible.
        """
        if matrix is None:
            tmpMatrix = spmatrix.ll_mat(1, 1, 1)
            try:
                assert len(nonZerosPerRow) == rows
                sizeHint = numerix.sum(numerix.asarray(nonZerosPerRow))
            except TypeError:
                sizeHint = nonZerosPerRow * rows
            if hasattr(tmpMatrix, 'storeZeros'):
                matrix = spmatrix.ll_mat(rows, cols, sizeHint, storeZeros)
            else:
                matrix = spmatrix.ll_mat(rows, cols, sizeHint)

        super(_PysparseMatrixFromShape, self).__init__(matrix=matrix)

class _PysparseBaseMeshMatrix(_PysparseMatrixFromShape):
    def __init__(self, mesh, rows, cols,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_PysparseMatrixFromShape` associated with a `Mesh`.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        rows : int
            The number of local matrix rows.
        cols : int
            The number of local matrix columns.
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~pysparse.spmatrix.ll_mat
            Pre-assembled SciPy matrix to use for storage.
        storeZeros : bool
            Instructs scipy to store zero values if possible.
        """
        self.mesh = mesh

        super(_PysparseBaseMeshMatrix, self).__init__(rows=rows,
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

    def put(self, vector, id1, id2, overlapping=False):
        """Insert local overlapping values and coordinates into global

        Parameters
        ----------
        vector : array_like
            The overlapping values to insert.
        id1 : array_like
            The local overlapping row indices.
        id2 : array_like
            The local overlapping column indices.
        overlapping : bool
            Whether to insert ghosted values or not (Ignored)
        """
        super(_PysparseBaseMeshMatrix, self).put(vector=vector, id1=id1, id2=id2)

    def addAt(self, vector, id1, id2, overlapping=False):
        """Accumulate local overlapping values and coordinates into global

        Parameters
        ----------
        vector : array_like
            The overlapping values to insert.
        id1 : array_like
            The local overlapping row indices.
        id2 : array_like
            The local overlapping column indices.
        overlapping : bool
            Whether to add ghosted values or not (Ignored)
        """
        super(_PysparseBaseMeshMatrix, self).addAt(vector=vector, id1=id1, id2=id2)

class _PysparseRowMeshMatrix(_PysparseBaseMeshMatrix):
    def __init__(self, mesh, cols, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_PysparseBaseMeshMatrix` with rows associated with equations.

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
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~pysparse.spmatrix.ll_mat
            Pre-assembled Pysparse matrix to use for storage.
        storeZeros : bool
            Instructs Pysparse to store zero values if possible.
        """
        self.numberOfEquations = numberOfEquations

        super(_PysparseRowMeshMatrix, self).__init__(mesh=mesh,
                                                     rows=numberOfEquations * mesh.numberOfCells,
                                                     cols=cols,
                                                     nonZerosPerRow=nonZerosPerRow,
                                                     exactNonZeros=exactNonZeros,
                                                     matrix=matrix,
                                                     storeZeros=storeZeros)

class _PysparseColMeshMatrix(_PysparseBaseMeshMatrix):
    def __init__(self, mesh, rows, numberOfVariables=1,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_PysparseBaseMeshMatrix` with columns associated with solution variables.

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
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~pysparse.spmatrix.ll_mat
            Pre-assembled Pysparse matrix to use for storage.
        storeZeros : bool
            Instructs Pysparse to store zero values if possible.
        """
        self.numberOfVariables = numberOfVariables

        super(_PysparseColMeshMatrix, self).__init__(mesh=mesh,
                                                     rows=rows,
                                                     cols=numberOfVariables * mesh.numberOfCells,
                                                     nonZerosPerRow=nonZerosPerRow,
                                                     exactNonZeros=exactNonZeros,
                                                     matrix=matrix,
                                                     storeZeros=storeZeros)

class _PysparseMeshMatrix(_PysparseRowMeshMatrix):
    def __init__(self, mesh, numberOfVariables=1, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, storeZeros=True):
        """Creates a `_PysparseBaseMeshMatrix` with associated with equations and variables.

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
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~pysparse.spmatrix.ll_mat
            Pre-assembled Pysparse matrix to use for storage
        storeZeros : bool
            Instructs Pysparse to store zero values if possible.
        """
        self.numberOfVariables = numberOfVariables

        super(_PysparseMeshMatrix, self).__init__(mesh=mesh,
                                                  cols=numberOfVariables * mesh.numberOfCells,
                                                  numberOfEquations=numberOfEquations,
                                                  nonZerosPerRow=nonZerosPerRow,
                                                  exactNonZeros=exactNonZeros,
                                                  matrix=matrix,
                                                  storeZeros=storeZeros)

    def __mul__(self, other):
        if isinstance(other, _PysparseMeshMatrix):
            return _PysparseMeshMatrix(mesh=self.mesh,
                                       matrix=spmatrix.matrixmultiply(self.matrix, other.matrix))
        else:
            return _PysparseMatrixFromShape.__mul__(self, other)

    def asTrilinosMeshMatrix(self):
        """Transforms a Pysparse matrix into a Trilinos matrix

        Maintains the Trilinos matrix as an attribute.

        Returns
        -------
        ~fipy.matrices.trilinosMatrix._TrilinosMatrix
        """
        A = self.matrix.copy()
        values, irow, jcol = A.find()

        if not hasattr(self, 'trilinosMatrix'):
            if A.shape[0] == 0:
                nonZerosPerRow = 0
            else:
                nonZerosPerRow = int(numerix.ceil(float(len(values)) / float(A.shape[0])))
            nonZerosPerRow = 1
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrixKeepStencil
            tmp = _TrilinosMeshMatrixKeepStencil(mesh=self.mesh, nonZerosPerRow=nonZerosPerRow,
                                                 numberOfVariables=self.numberOfVariables,
                                                 numberOfEquations=self.numberOfEquations)
            self.trilinosMatrix = tmp

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
        """Deletes the pysparse matrix and calls `self.trilinosMatrix.flush()` if necessary.
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
        >>> print(not hasattr(m.matrix, 'storeZeros')
        ...       or numerix.allequal(list(m.matrix.keys()),
        ...                           [(0, 1), (1, 0), (2, 2)]))
        True
        >>> print(not hasattr(m.matrix, 'storeZeros')
        ...       or numerix.allequal(list(m.matrix.values()), [1., 2., 0.]))
        True
        >>> m = _PysparseMatrixFromShape(rows=3, cols=3, storeZeros=False)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> print(numerix.allequal(list(m.matrix.keys()), [(0, 1), (1, 0)]))
        True
        >>> print(numerix.allequal(list(m.matrix.values()), numerix.array([1.0, 2.0])))
        True

        Storing more than preallocated is not an error when `exactNonZeros` is set

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1, exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0]) # doctest: +IGNORE_EXCEPTION_DETAIL

        This is also true if multiple values are accumulated into the
        same matrix entry.

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1, exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,0], [2,1,1,1]) # doctest: +IGNORE_EXCEPTION_DETAIL

        Preallocation can be specified row-by-row

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3,
        ...                              nonZerosPerRow=[2, 1, 1])
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])

        Preallocating on the wrong rows is not an error

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3,
        ...                              nonZerosPerRow=[1, 2, 1])
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])

        even when `exactNonZeros` is specified.

        >>> m = _PysparseMatrixFromShape(rows=3, cols=3,
        ...                              nonZerosPerRow=[1, 2, 1],
        ...                              exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
        """
        pass

class _PysparseIdentityMatrix(_PysparseMatrixFromShape):
    """
    Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size):
        """Create a sparse matrix with `1` in the diagonal

            >>> print(_PysparseIdentityMatrix(size=3))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PysparseMatrixFromShape.__init__(self, rows=size, cols=size, nonZerosPerRow=1)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)

class _PysparseIdentityMeshMatrix(_PysparseIdentityMatrix):
    def __init__(self, mesh):
        """Create a sparse matrix associated with a `Mesh` with `1` in the diagonal

            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print(_PysparseIdentityMeshMatrix(mesh=mesh))
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
