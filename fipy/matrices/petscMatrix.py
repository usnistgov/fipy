#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "petscMatrix.py"
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

from petsc4py import PETSc

from fipy.tools import numerix
from fipy.matrices.sparseMatrix import _SparseMatrix

class _PETScMatrix(_SparseMatrix):
    
    def __init__(self, matrix,
                 rowMap=None, colMap=None):
        """Creates a wrapper for a PETSc matrix

        Allows basic python operations __add__, __sub__ etc.
        Facilitate matrix populating in an easy way.
        
        :Parameters:
          - `matrix`: The starting `PETSc.Mat` 
        """
        self.matrix = matrix
        
        self.rowMap = rowMap
        self.colMap = colMap
        
        if self.rowMap is not None:
            self.matrix.setLGMap(rmap=rowMap, cmap=colMap)
   
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

        from fipy.tools import parallelComm
        return ''.join(parallelComm.allgather(_SparseMatrix.__str__(self)))

    def __iadd__(self, other):
        if other != 0:
            self.matrix.assemblyBegin()
            self.matrix.assemblyEnd()
            other.matrix.assemblyBegin()
            other.matrix.assemblyEnd()
            self.matrix = self.matrix + other.matrix
        return self

    def __add__(self, other):
        """
        Add two sparse matrices
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L + _PETScIdentityMatrix(size=3)
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
            elif shape == (N,):
                y = PETSc.Vec().createWithArray(other, comm=PETSc.COMM_WORLD)
                result = y.duplicate()
                self.matrix.mult(y, result)
                return result
            else:
                raise TypeError
            
    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            N = self._shape[1]
            x = PETSc.Vec().createMPI(N, comm=PETSc.COMM_WORLD)
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
        return range(self._shape[1]), range(self._shape[0])
        
    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=2)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> print L
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
        # FIXME: This is a terrible implementation
        self.matrix.assemblyBegin()
        self.matrix.assemblyEnd()
        vector = [self.matrix[i, j] for i, j in zip(id1, id2)]
        vector = numerix.array(vector, 'd')
        return vector

    def takeDiagonal(self):
        ids = numerix.arange(self._shape[0])
        return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)
        
            >>> L = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=3)
            >>> L.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L.addAt([1.73,2.2,8.4,3.9,1.23], [1,2,0,0,1], [2,2,0,0,2])
            >>> print L
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
        viewer = PETSc.Viewer().createASCII(name=filename, mode='w', 
                                            format=PETSc.Viewer.Format.ASCII_MATRIXMARKET)
        viewer.view(obj=self.matrix)
        viewer.destroy()
    
class _PETScMatrixFromShape(_PETScMatrix):
    
    def __init__(self, rows, cols, bandwidth=0, sizeHint=None, matrix=None,
                 rowMap=None, colMap=None):
        """Instantiates and wraps a PETSc `Mat` matrix

        :Parameters:
          - `rows`: The number of matrix rows
          - `cols`: The number of matrix columns
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `rowMap`: the map for the rows that this processor holds
          - `colMap`: the map for the colums that ???
          
        """
        bandwidth = bandwidth 
        if (bandwidth == 0) and (sizeHint is not None):
            bandwidth = sizeHint / max(rows, cols)
        if matrix is None:
            matrix = PETSc.Mat()
            if rowMap is None:
                if colMap is None:
                    comm = PETSc.COMM_WORLD
                else:
                    comm = colMap.comm
            else:
                comm = rowMap.comm
            matrix.create(comm)
            # rows are owned per process
            # cols are owned by everyone
            matrix.setSizes([[rows, None], [cols, None]])
            matrix.setType('aij') # sparse
            matrix.setPreallocationNNZ(None) # FIXME: ??? #bandwidth)
                
        _PETScMatrix.__init__(self, matrix=matrix, rowMap=rowMap, colMap=colMap)

class _PETScMeshMatrix(_PETScMatrixFromShape):
    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, numberOfVariables=1, numberOfEquations=1):

        """Creates a `_PETScMatrixFromShape` associated with a `Mesh`. Allows for different number of equations and/or variables

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `matrix`: pre-assembled `ll_mat` to use for storage
          - `numberOfVariables`: The columns of the matrix is determined by numberOfVariables * self.mesh.numberOfCells.
          - `numberOfEquations`: The rows of the matrix is determined by numberOfEquations * self.mesh.numberOfCells.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations

        comm = mesh.communicator.petsc4py_comm
        rowMap = PETSc.LGMap().create(self._globalNonOverlappingRowIDs.astype('int32'), comm=comm)
        colMap = PETSc.LGMap().create(self._globalOverlappingColIDs.astype('int32'), comm=comm)
        
        _PETScMatrixFromShape.__init__(self, 
                                       rows=numberOfEquations * len(self.mesh._localNonOverlappingCellIDs), 
                                       cols=numberOfVariables * self.mesh.globalNumberOfCells, 
                                       bandwidth=bandwidth, 
                                       sizeHint=sizeHint, 
                                       matrix=matrix,
                                       rowMap=rowMap,
                                       colMap=colMap)

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
        return range(0, self.numberOfVariables, self.mesh.globalNumberOfCells)
                     
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
        id1 = id1[mask]
        id2 = id2[mask]
        
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
            return _PETScMatrixFromShape.__mul__(self, other=other)
        else:
            shape = numerix.shape(other)

            if shape == ():
                result = self.copy()
                result.matrix = self.matrix * other
            else:
                x = PETSc.Vec().createMPI(N, comm=PETSc.COMM_WORLD)
                if self.colMap is not None:
                    x.setLGMap(self.colMap)
                y = x.duplicate()
                localNonOverlappingColIDs = self._localNonOverlappingColIDs.astype('int32')
                x.setValuesLocal(localNonOverlappingColIDs, other[localNonOverlappingColIDs])
                self.matrix.mult(x, y)
                return y
        
    @property
    def numpyArray(self):
        return super(_PETScMeshMatrix, self).numpyArray

    def flush(self):
        """
        Deletes the copy of the PETSc matrix held.
        """
    
        if (not hasattr(self, 'cache')) or (self.cache is False):
            del self.matrix

    def _test(self):
        """
        Tests
        
        >>> m = _PETScMatrixFromShape(rows=3, cols=3, bandwidth=1)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> m.matrix.assemblyBegin()
        >>> m.matrix.assemblyEnd()
        
        # FIXME: are these names even right? is this a good test?
        
        >>> col, row, val = m.matrix.getValuesCSR()
        >>> print numerix.allequal(col, [0, 1, 2, 3])
        True
        >>> print numerix.allequal(row, [1, 0, 2])
        True
        >>> print numerix.allclose(val, [1., 2., 0.])
        True
        """
        pass
        
class _PETScIdentityMatrix(_PETScMatrixFromShape):
    """
    Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size, bandwidth=1):
        """Create a sparse matrix with '1' in the diagonal
        
            >>> print _PETScIdentityMatrix(size=3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScMatrixFromShape.__init__(self, rows=size, cols=size, bandwidth=bandwidth)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)
        
class _PETScIdentityMeshMatrix(_PETScIdentityMatrix):
    def __init__(self, mesh, bandwidth=1):
        """Create a sparse matrix associated with a `Mesh` with '1' in the diagonal
        
            >>> from fipy import Grid1D
            >>> from fipy.tools import serialComm
            >>> mesh = Grid1D(nx=3, communicator=serialComm)
            >>> print _PETScIdentityMeshMatrix(mesh=mesh)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PETScIdentityMatrix.__init__(self, size=mesh.numberOfCells, bandwidth=bandwidth)

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test()
