from __future__ import division
from builtins import zip
from builtins import range
__docformat__ = 'restructuredtext'

__all__ = []

from petsc4py import PETSc

from fipy.tools import numerix
from fipy.matrices.sparseMatrix import _SparseMatrix

class _PETScMatrix(_SparseMatrix):
    
    def __init__(self, matrix):
        """Creates a wrapper for a PETSc matrix

        Allows basic python operations __add__, __sub__ etc.
        Facilitate matrix populating in an easy way.
        
        :Parameters:
          - `matrix`: The starting `PETSc.Mat` 
        """
        self.matrix = matrix
        
    def getCoupledClass(self):
        return _CoupledPETScMeshMatrix
    
    def copy(self):
        return _PETScMatrix(matrix=self.matrix.copy())
        
    def __getitem__(self, index):
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()
        m = self.matrix[index]
        if numerix.shape(m) == ():
            return m
        else:
            return _PETScMatrix(matrix=m)
    
    def __str__(self):
        self.matrix.assemble()

        return _SparseMatrix.__str__(self)

    def __iadd__(self, other):
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()
        if other != 0:
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
            self.matrix = self.matrix + other.matrix
        return self

    def __add__(self, other):
        """
        Add two sparse matrices
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=3)
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
            self.matrix.assemblyBegin()
            self.matrix.assemblyEnd()
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
            return _PETScMatrix(matrix=self.matrix + other.matrix)
        
    __radd__ = __add__
    
    def __sub__(self, other):

        if other == 0:
            return self
        else:
            self.matrix.assemblyBegin()
            self.matrix.assemblyEnd()
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
            return _PETScMatrix(matrix=self.matrix - other.matrix)

    def __rsub__(self, other):
        return -self + other
    
    def __isub__(self, other):
        if other != 0:
            self.matrix.assemblyBegin()
            self.matrix.assemblyEnd()
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
            self.matrix = self.matrix - other.matrix
        return self

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=2)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PETScIdentityMatrix(size=3, bandwidth=3)
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
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp) ## The multiplication is broken. Numpy is calling __rmul__ for every element instead of with  the whole array.
            1

            
        """
        N = self._shape[1]

        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()

        if isinstance(other, _PETScMatrix):
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
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
                return result
            else:
                raise TypeError
            
    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            N = self._shape[1]
            x = PETSc.Vec().createMPI(N, comm=self.matrix.comm)
            y = x.duplicate()
            x[:] = other
            self.matrix.multTranspose(x, y)
            return numerix.asarray(y)
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
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=2)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print(L)
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
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
        row_ptr : array_like
            locations in the val vector that start a row, 
            terminated with len(val) + 1
        cols : array_like
            column indices
        val : array_like
            non-zero values
        """
#         from fipy.tools.debug import PRINT
#         PRINT("i:", i)
#         PRINT("j:", j)
#         PRINT("v:", v)
        
        i = numerix.asarray(i)
        j = numerix.asarray(j)
        v = numerix.asarray(v)
        start_row, end_row = self.matrix.getOwnershipRange()
#         PRINT("start_row:", start_row)
#         PRINT("end_row:", end_row)

        ix = numerix.lexsort([i, j])
        ix = ix[(j[ix] >= start_row) & (j[ix] < end_row)]
        cols = i[ix]
        row_ptr = numerix.searchsorted(j[ix], 
#                                        numerix.arange(0, self._shape[1]+1))
                                       numerix.arange(start_row, end_row + 1))
        vals = v[ix]
        
#         PRINT("i[ix]:", i[ix])
#         PRINT("j[ix]:", j[ix])
#         PRINT("v[ix]:", v[ix])

        
        
#         PRINT("size:", self._shape)
#         PRINT("shape:", self._shape[0])
#         
#         PRINT("row_ptr:", row_ptr)
#         PRINT("cols:", cols)
#         PRINT("vals:", vals)
#         
#         PRINT("size:", self.matrix.sizes)
#         PRINT("owns-rows:", self.matrix.getOwnershipRange())
#         PRINT("owns-cols:", self.matrix.getOwnershipRangesColumn())
        
        # note: PETSc (at least via pip) only seems to handle 32 bit addressing
        return row_ptr.astype('int32'), cols.astype('int32'), vals

    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=1)
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
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._shape[0])
            tmp = numerix.zeros((self._shape[0],), 'd')
            tmp[:] = vector
            self.put(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.put(vector, ids, ids)

    def take(self, id1, id2):
        # FIXME: This is a terrible implementation
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()
        vector = [self.matrix[i, j] for i, j in zip(id1, id2)]
        vector = numerix.array(vector, 'd')
        return vector

    def takeDiagonal(self):
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()

        return self.matrix.getDiagonal().array

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
            >>> print(L)
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        self.matrix.setValuesCSR(*self._ijv2csr(id2, id1, vector),
                                 addv=True)

    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._shape[0])
            tmp = numerix.zeros((self._shape[0],), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)
            
    def _petsc2app(self, ids):
        return ids

    @property
    def _scipy_csr(self):
        """Return the PETSc-ordered CSR matrix
        """
        from scipy import sparse
        mpi4pycomm = self.matrix.comm.tompi4py()

        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()

        indptr, indices, data = self.matrix.getValuesCSR()

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
                                  (self._petsc2app(coo.row),
                                   self._petsc2app(coo.col))),
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
    
class _PETScMatrixFromShape(_PETScMatrix):
    
    def __init__(self, rows, cols, bandwidth=0, sizeHint=None, matrix=None, comm=PETSc.COMM_SELF):
        """Instantiates and wraps a PETSc `Mat` matrix

        :Parameters:
          - `rows`: The number of matrix rows
          - `cols`: The number of matrix columns
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `comm`: communicator
          
        """
        bandwidth = bandwidth 
        if (bandwidth == 0) and (sizeHint is not None):
            bandwidth = sizeHint // max(rows, cols)
        if matrix is None:
            matrix = PETSc.Mat()
            matrix.create(comm)
            # rows are owned per process
            # cols are owned by everyone
            matrix.setSizes([[rows, None], [cols, None]])
            matrix.setType('aij') # sparse
            matrix.setPreallocationNNZ(None) # FIXME: ??? #bandwidth)
                
        _PETScMatrix.__init__(self, matrix=matrix)

class _PETScMeshMatrix(_PETScMatrixFromShape):
    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, numberOfVariables=1, numberOfEquations=1):

        """Creates a `_PETScMatrixFromShape` associated with a `Mesh`. Allows for different number of equations and/or variables

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `numberOfVariables`: The columns of the matrix is determined by `numberOfVariables * self.mesh.numberOfCells`.
          - `numberOfEquations`: The rows of the matrix is determined by `numberOfEquations * self.mesh.numberOfCells`.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations

        _PETScMatrixFromShape.__init__(self, 
                                       rows=numberOfEquations * len(self.mesh._localNonOverlappingCellIDs), 
                                       cols=numberOfVariables * len(self.mesh._localNonOverlappingCellIDs), # self.mesh.globalNumberOfCells, # , # 
                                       bandwidth=bandwidth, 
                                       sizeHint=sizeHint, 
                                       matrix=matrix,
                                       comm=mesh.communicator.petsc4py_comm)
                                       
    
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

            fipyIDs = self._globalNonOverlappingColIDs
            N = len(fipyIDs)

            count = numerix.zeros((comm.Nproc,), dtype=int)
            count[comm.procID] = N
            comm.mpi4py_comm.Allreduce(sendbuf=MPI.IN_PLACE, recvbuf=count, op=MPI.MAX)
            
            petscIDs = numerix.arange(N) + numerix.sum(count[:comm.procID])
            
            self._ao_ = PETSc.AO().createBasic(petsc=petscIDs.astype('int32'), 
                                               app=fipyIDs.astype('int32'), 
                                               comm=comm.petsc4py_comm)
        return self._ao_
        
    def _petsc2app(self, ids):
        return self._ao.petsc2app(ids)

    def _cellIDsToGlobalRowIDs(self, IDs):
         N = len(IDs)
         M = self.numberOfEquations
         return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.globalNumberOfCells).flatten()

    def _cellIDsToGlobalColIDs(self, IDs):
         N = len(IDs)
         M = self.numberOfVariables
         return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.globalNumberOfCells).flatten()

    def _cellIDsToLocalRowIDs(self, IDs):
         M = self.numberOfEquations
         N = len(IDs)
         return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.numberOfCells).flatten()

    def _cellIDsToLocalColIDs(self, IDs):
         M = self.numberOfVariables
         N = len(IDs)
         return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.numberOfCells).flatten()

    @property
    def _globalNonOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def _globalNonOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def _globalOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def _globalCommonColIDs(self):
        return list(range(0, self.numberOfVariables, self.mesh.globalNumberOfCells))
                     
    @property
    def _globalOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def _localNonOverlappingRowIDs(self):
        return self._cellIDsToLocalRowIDs(self.mesh._localNonOverlappingCellIDs)

    @property
    def _localOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localOverlappingCellIDs)

    @property
    def _localNonOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localNonOverlappingCellIDs)

    def copy(self):
        tmp = _PETScMatrixFromShape.copy(self)
        copy = self.__class__(mesh=self.mesh) # FIXME: ??? , bandwidth=self.bandwidth)
        copy.matrix = tmp.matrix
        return copy

    def _getStencil(self, id1, id2):
        id1 = self._globalOverlappingRowIDs[id1]
        id2 = self._globalOverlappingColIDs[id2]
            
        mask = numerix.in1d(id1, self._globalNonOverlappingRowIDs) 
        id1 = self._ao.app2petsc(id1[mask].astype('int32'))
        id2 = self._ao.app2petsc(id2[mask].astype('int32'))
        
        return id1, id2, mask

    def _globalNonOverlapping(self, vector, id1, id2):
        """Transforms and subsets local overlapping values and coordinates to global non-overlapping
        
        :Parameters:
          - `vector`: The overlapping values to insert.
          - `id1`: The local overlapping row indices.
          - `id2`: The local overlapping column indices.
          
        :Returns: 
          Tuple of (non-overlapping vector, 
                    global non-overlapping row indices, 
                    global non-overlapping column indices)
        """
        id1, id2, mask = self._getStencil(id1, id2)
        vector = vector[mask]
        return (vector, id1, id2)

    def put(self, vector, id1, id2):
        vector, id1, id2 = self._globalNonOverlapping(vector, id1, id2)
        _PETScMatrixFromShape.put(self, vector=vector, id1=id1, id2=id2)

    def addAt(self, vector, id1, id2):
        vector, id1, id2 = self._globalNonOverlapping(vector, id1, id2)
        _PETScMatrixFromShape.addAt(self, vector=vector, id1=id1, id2=id2)
    
    @property
    def _bodies(self):
        if not hasattr(self, "_bodies_"):
            self._bodies_ = numerix.in1d(self.mesh._globalOverlappingCellIDs, 
                                         self.mesh._globalNonOverlappingCellIDs)
        return self._bodies_
        
    @property
    def _ghosts(self):
        if not hasattr(self, "_ghosts_"):
            self._ghosts_ = self.mesh._globalOverlappingCellIDs[~self._bodies]
            self._ghosts_ = self._cellIDsToGlobalRowIDs(self._ghosts_)
            self._ghosts_ = self._ao.app2petsc(self._ghosts_.astype('int32'))
            
        return self._ghosts_
        
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
        corporeal = numerix.asarray(var[..., self._bodies]).ravel()
        incorporeal = numerix.asarray(var[..., ~self._bodies]).ravel()
        array = numerix.concatenate([corporeal, incorporeal])

        comm = self.mesh.communicator.petsc4py_comm
        vec = PETSc.Vec().createGhostWithArray(ghosts=self._ghosts.astype('int32'),
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
        M = self.numberOfEquations
        var = numerix.empty((M, N))
        bodies = numerix.array(vec)
        if M > 1:
            bodies = numerix.reshape(bodies, (M, -1))
        var[..., self._bodies] = bodies
        vec.ghostUpdate()
        with vec.localForm() as lf:
            if len(self._ghosts) > 0:
                ids = numerix.arange(-len(self._ghosts), 0)
                ghosts = numerix.reshape(numerix.array(lf)[ids], (M, -1))
                var[..., ~self._bodies] = ghosts

        return var.flatten()

    def __mul__(self, other):
        """Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=2)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PETScIdentityMatrix(size=3, bandwidth=3)
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
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp) ## The multiplication is broken. Numpy is calling __rmul__ for every element instead of with  the whole array.
            1

            
        """
        N = self._shape[1]

        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()
        
        if isinstance(other, (_PETScMatrix, PETSc.Vec)):
            return _PETScMatrixFromShape.__mul__(self, other=other)
        else:
            shape = numerix.shape(other)

            if shape == ():
                result = self.copy()
                result.matrix = self.matrix * other
                return result
            else:
                x = other[self._localNonOverlappingColIDs]
                x = PETSc.Vec().createWithArray(x, comm=self.matrix.comm)

                y = PETSc.Vec().createGhost(ghosts=self._ghosts.astype('int32'),
                                            size=(len(self._localNonOverlappingColIDs), None),
                                            comm=self.matrix.comm)
                self.matrix.mult(x, y)
                return self._petsc2fipyGhost(vec=y)

    def takeDiagonal(self):
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()

        y = PETSc.Vec().createGhost(ghosts=self._ghosts.astype('int32'),
                                    size=(len(self._localNonOverlappingColIDs), None),
                                    comm=self.matrix.comm)
        self.matrix.getDiagonal(result=y)
        return self._petsc2fipyGhost(vec=y)

    def flush(self):
        """Deletes the copy of the PETSc matrix held.
        """
    
        if not getattr(self, 'cache', False):
            del self.matrix

    def _test(self):
        """Tests
        
        >>> m = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=1)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> m.matrix.assemblyBegin()
        >>> m.matrix.assemblyEnd()
        
        # FIXME: are these names even right? is this a good test?
        
        >>> col, row, val = m.matrix.getValuesCSR()
        >>> print(numerix.allequal(col, [0, 1, 2, 3]))
        True
        >>> print(numerix.allequal(row, [1, 0, 2]))
        True
        >>> print(numerix.allclose(val, [1., 2., 0.]))
        True
        """
        pass
        
class _PETScIdentityMatrix(_PETScMatrixFromShape):
    """Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size, bandwidth=1, comm=PETSc.COMM_SELF):
        """Create a sparse matrix with `1` in the diagonal
        
            >>> print(_PETScIdentityMatrix(size=3))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScMatrixFromShape.__init__(self, rows=size, cols=size, bandwidth=bandwidth, comm=comm)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)
        
class _PETScIdentityMeshMatrix(_PETScIdentityMatrix):
    def __init__(self, mesh, bandwidth=1):
        """Create a sparse matrix associated with a `Mesh` with `1` in the diagonal
        
            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print(_PETScIdentityMeshMatrix(mesh=mesh))
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScIdentityMatrix.__init__(self, size=mesh.numberOfCells, bandwidth=bandwidth, 
                                      comm=mesh.communicator.petsc4py_comm)

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test()
