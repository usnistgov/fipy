from __future__ import division
from builtins import zip
from builtins import range
__docformat__ = 'restructuredtext'

__all__ = []

from collections.abc import Iterable
from petsc4py import PETSc

from fipy.tools import numerix
from fipy.matrices.sparseMatrix import (_SparseMatrix, _RowMesh2Matrix,
                                        _ColMesh2Matrix, _RowColMesh2Matrix)

class _PETScMatrix(_SparseMatrix):

    def __init__(self, matrix):
        """Creates a wrapper for a PETSc matrix

        Allows basic python operations __add__, __sub__ etc.
        Facilitate matrix populating in an easy way.

        :Parameters:
          - `matrix`: The starting `PETSc.Mat`
        """
        self.matrix = matrix

    def __del__(self):
        self.matrix.destroy()

    def copy(self):
        return _PETScMatrix(matrix=self.matrix.copy())

    def __getitem__(self, index):
        self.matrix.assemble()
        m = self.matrix[index]
        if numerix.shape(m) != ():
            m = _PETScMatrix(matrix=m)
        return m

    def __str__(self):
        self.matrix.assemble()

        return _SparseMatrix.__str__(self)

    def __iadd__(self, other):
        self.matrix.assemble()
        if other != 0:
            other.matrix.assemble()
            self.matrix = self.matrix + other.matrix
        return self

    def __add__(self, other):
        """
        Add two sparse matrices

            >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print(L + _PETScIdentityMatrix(size=3))
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
            self.matrix.assemble()
            other.matrix.assemble()
            return _PETScMatrix(matrix=self.matrix + other.matrix)

    __radd__ = __add__

    def __sub__(self, other):

        if other == 0:
            return self
        else:
            self.matrix.assemble()
            other.matrix.assemble()
            return _PETScMatrix(matrix=self.matrix - other.matrix)

    def __rsub__(self, other):
        return -self + other

    def __isub__(self, other):
        if other != 0:
            self.matrix.assemble()
            other.matrix.assemble()
            self.matrix = self.matrix - other.matrix
        return self

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix

            >>> L1 = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=2)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PETScIdentityMatrix(size=3, nonZerosPerRow=3)
            >>> L2.put([4.38], [2], [1])
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
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp)
            1

        (The multiplication is broken.  Numpy is calling __rmul__ for every
        element instead of with the whole array.)
        """
        N = self._shape[1]

        self.matrix.assemble()

        if isinstance(other, _PETScMatrix):
            other.matrix.assemble()
            copy = self.copy()
            copy.matrix = self.matrix.matMult(other.matrix)
            return copy
        elif isinstance(other, PETSc.Vec):
            y = other.duplicate()
            self.matrix.mult(other, y)
            return y
        else:
            shape = numerix.shape(other)
            if shape == ():
                result = self.copy()
                result.matrix = self.matrix * other
                return result
            elif shape == (N,):
                y = PETSc.Vec().createWithArray(other, comm=self.matrix.comm)
                result = y.duplicate()
                self.matrix.mult(y, result)
                y.destroy()
                return result
            else:
                raise TypeError

    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            N = self._shape[1]
            x = PETSc.Vec().createMPI(N, comm=self.matrix.comm)
            y = x.duplicate()
            x[:] = other
            self.matrix.assemble()
            x.assemble()
            self.matrix.multTranspose(x, y)
            arr = numerix.asarray(y)
            x.destroy()
            y.destroy()
            return arr
        else:
            return self * other

    @property
    def _shape(self):
        return self.matrix.sizes[0][0], self.matrix.sizes[1][1]

    @property
    def _range(self):
        return list(range(self._shape[1])), list(range(self._shape[0]))

    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)

            >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=2)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print(L)
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        self.matrix.assemble(self.matrix.AssemblyType.FLUSH)
        self.matrix.setValuesCSR(*self._ijv2csr(id2, id1, vector))

    def _ijv2csr(self, i, j, v):
        """Convert arrays of matrix indices and values into CSR format

        see: http://netlib.org/linalg/html_templates/node91.html#SECTION00931100000000000000

        .. note::
           petsc4py only understands CSR formatted matrices (setValuesCSR and
           setValuesIJV both inexplicably call the same underlying routine).

        Parameters
        ----------
        i : array_like
            column indices
        j : array_like
            row indices
        v : array_like
            non-zero values

        Returns
        -------
        ptrs : array_like
            locations in the val vector that start a row,
            terminated with len(val) + 1
        cols : array_like
            column indices
        data : array_like
            non-zero values
        """
        i = numerix.asarray(i)
        j = numerix.asarray(j)
        v = numerix.asarray(v)
        start_row, end_row = self.matrix.getOwnershipRange()

        ix = numerix.lexsort([i, j])
        ix = ix[(j[ix] >= start_row) & (j[ix] < end_row)]
        cols = i[ix]
        row_ptr = numerix.searchsorted(j[ix],
                                       numerix.arange(start_row, end_row + 1))
        vals = v[ix]

        # note: PETSc (at least via pip) only seems to handle 32 bit addressing
        return row_ptr.astype('int32'), cols.astype('int32'), vals

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix

            >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1)
            >>> L.putDiagonal([3.,10.,numerix.pi])
            >>> print(L)
             3.000000      ---        ---    
                ---    10.000000      ---    
                ---        ---     3.141593  
            >>> L.putDiagonal([10.,3.])
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
        # FIXME: This is a terrible implementation
        self.matrix.assemble()
        vector = [self.matrix[i, j] for i, j in zip(id1, id2)]
        vector = numerix.array(vector, 'd')
        return vector

    def takeDiagonal(self):
        self.matrix.assemble()

        return self.matrix.getDiagonal().array

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)

            >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
            >>> print(L)
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        self.matrix.assemble(self.matrix.AssemblyType.FLUSH)
        self.matrix.setValuesCSR(*self._ijv2csr(id2, id1, vector),
                                 addv=True)

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

        >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
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
        self.matrix.assemble()

        return self.matrix.getValuesCSR()

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

        >>> L = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=3)
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
        ptrs, cols, csrdata = self.CSR

        rows = [cols[start:stop] for start, stop in zip(ptrs[:-1], ptrs[1:])]
        data = [csrdata[start:stop] for start, stop in zip(ptrs[:-1], ptrs[1:])]

        return rows, data

    @property
    def _scipy_csr(self):
        """Return the PETSc-ordered CSR matrix
        """
        from scipy import sparse
        mpi4pycomm = self.matrix.comm.tompi4py()

        indptr, indices, data = self.CSR

        # getValuesCSR() returns entries local to node
        # with node-relative indptr.
        # sparse.csr_matrix() requires all elements for construction
        # and global indptr

        offset = numerix.cumsum([0] + mpi4pycomm.allgather(len(data)))
        offset = mpi4pycomm.scatter(offset[:-1])

        indices = numerix.concatenate(mpi4pycomm.allgather(indices))
        data = numerix.concatenate(mpi4pycomm.allgather(data))

        # strip local end markers and append global end marker
        indptr = mpi4pycomm.allgather(indptr[:-1] + offset) + [[len(data)]]
        indptr = numerix.concatenate(indptr)

        (rows, globalRows), (cols, globalCols) = self.matrix.getSizes()

        return sparse.csr_matrix((data, indices, indptr),
                                 shape=(globalRows, globalCols))

    @property
    def _scipy_coo(self):
        """Return the application-ordered COO matrix
        """
        from scipy import sparse

        coo = self._scipy_csr.tocoo()
        return sparse.coo_matrix((coo.data,
                                  (self._matrix2mesh(coo.row),
                                   self._matrix2mesh(coo.col))),
                                 shape=coo.shape)

    @property
    def numpyArray(self):
        # self.matrix.getDenseArray() raises
        # [0] MatDenseGetArray() line 1782 in .../src/mat/impls/dense/seq/dense.c
        # [0] No support for this operation for this object type
        # [0] Cannot locate function MatDenseGetArray_C in object
        return self._scipy_coo.toarray()

    def matvec(self, x):
        """This method is required for scipy solvers.
        """
        return self * x

    def exportMmf(self, filename):
        """Exports the matrix to a Matrix Market file of the given `filename`.
        """
        from scipy.io import mmio

        mmio.mmwrite(filename, self._scipy_coo)

    @property
    def T(self):
        """Transpose matrix

        Returns
        -------
        ~fipy.matrices.petscMatrix._PETScMatrix

        Examples
        --------

        >>> import fipy as fp

        >>> mesh = fp.Grid1D(nx=10)
        >>> ids = fp.CellVariable(mesh=mesh, value=mesh._globalOverlappingCellIDs)

        >>> mat = _PETScColMeshMatrix(mesh=mesh, rows=1)
        >>> mat.put(vector=ids.value,
        ...         id1=[fp.parallelComm.procID] * mesh.numberOfCells,
        ...         id2=mesh._localOverlappingCellIDs,
        ...         overlapping=True)

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
        >>> print(mat.T.numpyArray) # doctest: +PARALLEL_2
        [[ 0.  0.]
         [ 1.  0.]
         [ 2.  0.]
         [ 3.  3.]
         [ 4.  4.]
         [ 5.  5.]
         [ 6.  6.]
         [ 0.  7.]
         [ 0.  8.]
         [ 0.  9.]]
        """
        self.matrix.assemble()
        return _PETScMatrix(matrix=self.matrix.transpose())

class _PETScMatrixFromShape(_PETScMatrix):

    def __init__(self, rows, cols,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None, comm=PETSc.COMM_SELF):
        """Instantiates and wraps a PETSc `Mat` matrix

        Parameters
        ----------
        rows : int
            The number of local matrix rows.
        cols : int
            The number of local matrix columns
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~petsc4py.PETSc.Mat
            Pre-assembled PETSc matrix to use for storage.
        comm : ~PETSc.Comm
            The MPI communicator to use.
        """
        if matrix is None:
            matrix = PETSc.Mat()
            matrix.create(comm)
            # rows are owned per process
            # cols are owned by everyone
            matrix.setSizes([[rows, None], [cols, None]])
            matrix.setType('aij') # sparse
            matrix.setUp()
            if isinstance(nonZerosPerRow, Iterable) or (nonZerosPerRow > 0):
                if isinstance(nonZerosPerRow, Iterable):
                    nonZerosPerRow = numerix.asarray(nonZerosPerRow, dtype=PETSc.IntType)
                matrix.setPreallocationNNZ(nonZerosPerRow)
                if not exactNonZeros:
                    matrix.setOption(matrix.Option.NEW_NONZERO_ALLOCATION_ERR, False)

        super(_PETScMatrixFromShape, self).__init__(matrix=matrix)

class _PETScBaseMeshMatrix(_PETScMatrixFromShape):
    def __init__(self, mesh, rows, cols, m2m,
                 nonZerosPerRow=0, exactNonZeros=False,
                 matrix=None):
        """Creates a `_PETScMatrixFromShape` associated with a `Mesh`.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        rows : int
            The number of local matrix rows.
        cols : int
            The number of local matrix columns.
        m2m : ~fipy.matrices.sparseMatrix._Mesh2Matrix
            Object to convert between mesh coordinates and matrix coordinates.
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~petsc4py.PETSc.Mat
            Pre-assembled PETSc matrix to use for storage.
        """
        self.mesh = mesh
        self._m2m = m2m

        super(_PETScBaseMeshMatrix, self).__init__(rows=rows,
                                                   cols=cols,
                                                   nonZerosPerRow=nonZerosPerRow,
                                                   exactNonZeros=exactNonZeros,
                                                   matrix=matrix,
                                                   comm=mesh.communicator.petsc4py_comm)

    def copy(self):
        tmp = super(_PETScBaseMeshMatrix, self).copy()
        copy = self.__class__(mesh=self.mesh) # FIXME: ??? , nonZerosPerRow=self.nonZerosPerRow)
        copy.matrix = tmp.matrix
        return copy

    def __del__(self):
        if hasattr(self, "_ao_"):
            self._ao_.destroy()
        super(_PETScBaseMeshMatrix, self).__del__()

    @property
    def _ao(self):
        """Application Ordering to relate FiPy matrix rows to PETSc matrix rows

        FiPy naturally blocks matrix rows, one set of Equations (or Variables) at a time.
        PETSc requires that all rows pertaining to a particular MPI node be contiguous.
        This PETSc `AO` (Application Ordering) object converts between them.

        Only needed for FiPy to PETSc. We can efficiently slice from PETSc to
        FiPy, but PETSc requires us to know the row IDs.
        """
        if not hasattr(self, "_ao_"):
            comm = self.mesh.communicator

            from mpi4py import MPI

            fipyIDs = self._m2m.globalNonOverlappingRowIDs
            N = len(fipyIDs)

            count = numerix.zeros((comm.Nproc,), dtype=int)
            count[comm.procID] = N
            comm.mpi4py_comm.Allreduce(sendbuf=MPI.IN_PLACE, recvbuf=count, op=MPI.MAX)

            petscIDs = numerix.arange(N) + numerix.sum(count[:comm.procID])

            self._ao_ = PETSc.AO().createBasic(petsc=petscIDs.astype('int32'),
                                               app=fipyIDs.astype('int32'),
                                               comm=comm.petsc4py_comm)

        return self._ao_

    def _matrix2mesh(self, ids):
        """Convert matrix row indices to mesh cell indices
        """
        return self._ao.petsc2app(ids)

    def _mesh2matrix(self, ids):
        """Convert mesh cell indices to matrix row indices
        """
        return self._ao.app2petsc(ids.astype('int32'))

    def _fipy2petscGhost(self, var):
        """Convert a FiPy Variable to a PETSc `GhostVec`

        Moves the ghosts to the end, as necessary.
        `var` may be coupled/vector and so moving the ghosts is a bit subtle.

        Given a (2x4) vector variable `vij`

        ```
        v00  v01 (v02)        processor 0
        v10  v11 (v12)

            (v01) v02 v03     processor 1
            (v11) v12 v13
        ```

        where i is the vector index and j is the global index.
        Elements in () are ghosted

        We end up with the `GhostVec`

        ```
        v00 v01 v10 v11 (v02) (v12)   [4, 6]  processor 0
        v02 v03 v12 v13 (v01) (v11)   [1, 3]  processor 1
        ```

        where the [a, b] are the global ghost indices
        """
        corporeal = numerix.asarray(var[..., self._m2m.bodies]).ravel()
        incorporeal = numerix.asarray(var[..., ~self._m2m.bodies]).ravel()
        array = numerix.concatenate([corporeal, incorporeal])

        comm = self.mesh.communicator.petsc4py_comm
        vec = PETSc.Vec().createGhostWithArray(ghosts=self._m2m.ghosts.astype('int32'),
                                               array=array,
                                               comm=comm)

        return vec

    def _petsc2fipyGhost(self, vec):
        """Convert a PETSc `GhostVec` to a FiPy Variable (form)

        Moves the ghosts from the end, as necessary.
        The return Variable may be coupled/vector and so moving the ghosts
        is a bit subtle.

        Given an 8-element `GhostVec` `vj`

        ```
        v0 v1 v2 v3 (v4) (v6)   [4, 6]  processor 0
        v4 v5 v6 v7 (v1) (v3)   [1, 3]  processor 1
        ```

        where j is the global index and the `[a, b]` are the global ghost
        indices. Elements in () are ghosted

        We end up with the (2x4) FiPy Variable

        ```
        v0  v1 (v4)        processor 0
        v2  v3 (v6)

           (v1) v4 v4      processor 1
           (v3) v6 v7
        ```
        """
        N = len(self.mesh._globalOverlappingCellIDs)
        M = self._m2m.numberOfEquations
        var = numerix.empty((M, N), dtype=vec.array.dtype)
        bodies = numerix.asarray(vec)
        if M > 1:
            bodies = numerix.reshape(bodies, (M, -1))
        var[..., self._m2m.bodies] = bodies
        vec.ghostUpdate()
        with vec.localForm() as lf:
            if len(self._m2m.ghosts) > 0:
                ids = numerix.arange(-len(self._m2m.ghosts), 0)
                ghosts = numerix.reshape(numerix.array(lf)[ids], (M, -1))
                var[..., ~self._m2m.bodies] = ghosts

        return var.flatten()

    def _getGhostedValues(self, var):
        """Obtain current ghost values from across processes

        Returns
        -------
        ndarray
            Ghosted values
        """
        vec = self._fipy2petscGhost(var)
        arr = self._petsc2fipyGhost(vec)
        vec.destroy()
        return arr

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
            Whether to insert ghosted values or not (default False)
        """
        vector, id1, id2 = self._m2m.globalVectorAndIDs(vector, id1, id2, overlapping)
        super(_PETScBaseMeshMatrix, self).put(vector=vector, id1=id1, id2=id2)

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
            Whether to add ghosted values or not (default False)
        """
        vector, id1, id2 = self._m2m.globalVectorAndIDs(vector, id1, id2, overlapping)
        super(_PETScBaseMeshMatrix, self).addAt(vector=vector, id1=id1, id2=id2)

class _PETScRowMeshMatrix(_PETScBaseMeshMatrix):
    def __init__(self, mesh, cols, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False, matrix=None, m2m=None):
        """Creates a `_PETScMatrixFromShape` with rows associated with equations.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        cols : int
            The number of matrix columns
        numberOfEquations : int
            The local rows of the matrix are determined by
            `numberOfEquations * mesh._localNonOverlappingCellIDs`.
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~petsc4py.PETSc.Mat
            Pre-assembled PETSc matrix to use for storage.
        m2m : ~fipy.matrices.sparseMatrix._RowMesh2Matrix
            Object to convert between mesh coordinates and matrix coordinates.
        """
        if m2m is None:
            m2m = _RowMesh2Matrix(mesh=mesh, matrix=self,
                                  numberOfEquations=numberOfEquations)

        super(_PETScRowMeshMatrix, self).__init__(mesh=mesh,
                                                  rows=numberOfEquations * len(mesh._localNonOverlappingCellIDs),
                                                  cols=cols,
                                                  m2m=m2m,
                                                  nonZerosPerRow=nonZerosPerRow,
                                                  exactNonZeros=exactNonZeros,
                                                  matrix=matrix)

class _PETScColMeshMatrix(_PETScBaseMeshMatrix):
    def __init__(self, mesh, rows, numberOfVariables=1,
                 nonZerosPerRow=0, exactNonZeros=False, matrix=None):
        """Creates a `_PETScMatrixFromShape` with columns associated with solution variables.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        rows : int
            The number of matrix rows.
        numberOfVariables : int
            The local columns of the matrix are determined by
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
        matrix : ~petsc4py.PETSc.Mat
            Pre-assembled PETSc matrix to use for storage.
        """
        m2m = _ColMesh2Matrix(mesh=mesh, matrix=self,
                              numberOfVariables=numberOfVariables)

        super(_PETScColMeshMatrix, self).__init__(mesh=mesh,
                                                  rows=rows,
                                                  cols=numberOfVariables * len(mesh._localNonOverlappingCellIDs),
                                                  m2m=m2m,
                                                  nonZerosPerRow=nonZerosPerRow,
                                                  exactNonZeros=exactNonZeros,
                                                  matrix=matrix)

class _PETScMeshMatrix(_PETScRowMeshMatrix):
    def __init__(self, mesh, numberOfVariables=1, numberOfEquations=1,
                 nonZerosPerRow=0, exactNonZeros=False, matrix=None):
        """Creates a `_PETScBaseMeshMatrix` associated with equations and variables.

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to assemble the matrix for.
        numberOfVariables : int
            The local columns of the matrix are determined by
            `numberOfVariables * len(mesh._localNonOverlappingCellIDs)`.
        numberOfEquations : int
            The local rows of the matrix are determined by
            `numberOfEquations * len(mesh._localNonOverlappingCellIDs)`.
        nonZerosPerRow : int or array_like of int
            The approximate number of sparse entries per row.  Either a
            typical number, or an iterable of values for each row
            (default: 0).
        exactNonZeros : bool
            Whether `nonZerosPerRow` is exact or approximate.
            Performance is improved preallocation is exact, but errors
            can result if additional allocations are necessary.
            (default: False).
        matrix : ~petsc4py.PETSc.Mat
            Pre-assembled PETSc matrix to use for storage.
        """
        m2m = _RowColMesh2Matrix(mesh=mesh, matrix=self,
                                 numberOfVariables=numberOfVariables,
                                 numberOfEquations=numberOfEquations)

        super(_PETScMeshMatrix, self).__init__(mesh=mesh,
                                               cols=numberOfVariables * len(mesh._localNonOverlappingCellIDs),
                                               numberOfEquations=numberOfEquations,
                                               nonZerosPerRow=nonZerosPerRow,
                                               exactNonZeros=exactNonZeros,
                                               matrix=matrix,
                                               m2m=m2m)

    def __mul__(self, other):
        """Multiply a sparse matrix by another sparse matrix

            >>> L1 = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=2)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PETScIdentityMatrix(size=3, nonZerosPerRow=3)
            >>> L2.put([4.38], [2], [1])
            >>> L2.put([4.38,12357.2,1.1], [2,1,0], [1,0,2])

            >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
            ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
            ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

            >>> numerix.allclose((L1 * L2).numpyArray, tmp)
            1

        or a sparse matrix by a vector

            >>> tmp = numerix.array((29., 6.28318531, 2.5))
            >>> vec = L1 * numerix.array((1,2,3),'d')
            >>> numerix.allclose(vec, tmp)
            1
            >>> _ = vec.destroy()

        or a vector by a sparse matrix

            >>> tmp = numerix.array((7.5, 16.28318531,  3.))
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp)
            1

        (The multiplication is broken.  Numpy is calling __rmul__ for every
        element instead of with the whole array.)
        """
        N = self._shape[1]

        self.matrix.assemble()

        if isinstance(other, (_PETScMatrix, PETSc.Vec)):
            return _PETScMatrixFromShape.__mul__(self, other=other)
        else:
            shape = numerix.shape(other)

            if shape == ():
                result = self.copy()
                result.matrix = self.matrix * other
                return result
            else:
                x = other[self._m2m.localNonOverlappingColIDs]
                x = PETSc.Vec().createWithArray(x, comm=self.matrix.comm)

                y = PETSc.Vec().createGhost(ghosts=self._m2m.ghosts.astype('int32'),
                                            size=(len(self._m2m.localNonOverlappingColIDs), None),
                                            comm=self.matrix.comm)
                self.matrix.mult(x, y)
                arr = self._petsc2fipyGhost(vec=y)
                x.destroy()
                y.destroy()
                return arr

    def takeDiagonal(self):
        self.matrix.assemble()

        y = PETSc.Vec().createGhost(ghosts=self._m2m.ghosts.astype('int32'),
                                    size=(len(self._m2m.localNonOverlappingColIDs), None),
                                    comm=self.matrix.comm)
        self.matrix.getDiagonal(result=y)
        arr = self._petsc2fipyGhost(vec=y)
        y.destroy()
        return arr

    def flush(self):
        """Deletes the copy of the PETSc matrix held.
        """

        if not getattr(self, 'cache', False):
            del self.matrix

    def _test(self):
        """Tests

        >>> m = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> m.matrix.assemble()

        # FIXME: are these names even right? is this a good test?

        >>> col, row, val = m.matrix.getValuesCSR()
        >>> print(numerix.allequal(col, [0, 1, 2, 3]))
        True
        >>> print(numerix.allequal(row, [1, 0, 2]))
        True
        >>> print(numerix.allclose(val, [1., 2., 0.]))
        True

        Storing more than preallocated is an error when `exactNonZeros` is set

        >>> m = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1, exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0]) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        petsc4py.PETSc.Error: error code 63

        This is also true if multiple values are accumulated into the
        same matrix entry.

        >>> m = _PETScMatrixFromShape(rows=3, cols=3, nonZerosPerRow=1, exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,0], [2,1,1,1]) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        petsc4py.PETSc.Error: error code 63

        Preallocation can be specified row-by-row

        >>> m = _PETScMatrixFromShape(rows=3, cols=3,
        ...                           nonZerosPerRow=[2, 1, 1])
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])

        Preallocating on the wrong rows is not an error

        >>> m = _PETScMatrixFromShape(rows=3, cols=3,
        ...                           nonZerosPerRow=[1, 2, 1])
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])

        but it is when `exactNonZeros` is specified.

        >>> m = _PETScMatrixFromShape(rows=3, cols=3,
        ...                           nonZerosPerRow=[1, 2, 1],
        ...                           exactNonZeros=True)
        >>> m.addAt([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0]) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        petsc4py.PETSc.Error: error code 63
        """
        pass

class _PETScIdentityMatrix(_PETScMatrixFromShape):
    """Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size, nonZerosPerRow=1, comm=PETSc.COMM_SELF):
        """Create a sparse matrix with `1` in the diagonal

            >>> print(_PETScIdentityMatrix(size=3))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScMatrixFromShape.__init__(self, rows=size, cols=size, nonZerosPerRow=nonZerosPerRow, comm=comm)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)

class _PETScIdentityMeshMatrix(_PETScIdentityMatrix):
    def __init__(self, mesh, nonZerosPerRow=1):
        """Create a sparse matrix associated with a `Mesh` with `1` in the diagonal

            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print(_PETScIdentityMeshMatrix(mesh=mesh))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScIdentityMatrix.__init__(self, size=mesh.numberOfCells, nonZerosPerRow=nonZerosPerRow,
                                      comm=mesh.communicator.petsc4py_comm)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
