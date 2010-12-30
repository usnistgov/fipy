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

from pysparse import spmatrix
from fipy.tools import numerix


from fipy.matrices.sparseMatrix import _SparseMatrix

class _PysparseMatrixBase(_SparseMatrix):
    
    """
    _PysparseMatrix class wrapper for pysparse.
    _PysparseMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, matrix):
        """Creates a `_PysparseMatrixBase`.

        :Parameters:
          - `matrix`: The starting `spmatrix` 
        """
        self.matrix = matrix

    def _getMatrix(self):
        return self.matrix
    
    def copy(self):
        return _PysparseMatrixBase(matrix=self.matrix.copy())
        
    def __getitem__(self, index):
        m = self.matrix[index]
        if type(m) is type(0) or type(m) is type(0.):
            return m
        else:
            return _PysparseMatrixBase(matrix=m)

    def __iadd__(self, other):
            return self._iadd(self._getMatrix(), other)
        
    def _iadd(self, L, other, sign = 1):
        if other != 0:
            L.shift(sign, other._getMatrix())
        return self

    def _add(self, other, sign = 1):
        L = self.matrix.copy()
        self._iadd(L, other, sign)
        return _PysparseMatrixBase(matrix=L)

    def __add__(self, other):
        """
        Add two sparse matrices
        
            >>> L = _PysparseMatrix(size=3)
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
            AttributeError: 'int' object has no attribute '_getMatrix'
        """

        if other is 0:
            return self
        else:
            L = self.matrix.copy()
            L.shift(1, other._getMatrix())
            return _PysparseMatrixBase(matrix=L)
        
    def __sub__(self, other):

        if other is 0:
            return self
        else:
            L = self.matrix.copy()
            L.shift(-1, other._getMatrix())
            return _PysparseMatrixBase(matrix=L)

    def __isub__(self, other):
            return self._iadd(self._getMatrix(), other, -1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = _PysparseMatrix(size=3)
            >>> L1.put([3.,10.,numerix.pi,2.5], [0,0,1,2], [2,1,1,0])
            >>> L2 = _PysparseIdentityMatrix(size=3)
            >>> L2.put([4.38,12357.2,1.1], [2,1,0], [1,0,2])
            
            >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
            ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
            ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

            >>> numerix.allclose((L1 * L2).getNumpyArray(), tmp)
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

        if isinstance(other, _PysparseMatrixBase):
            return _PysparseMatrixBase(matrix=spmatrix.matrixmultiply(self.matrix, other._getMatrix()))
        else:
            shape = numerix.shape(other)
            if shape == ():
                L = spmatrix.ll_mat(N, N, N)
                L.put(other * numerix.ones(N))
                return _PysparseMatrixBase(matrix=spmatrix.matrixmultiply(self.matrix, L))
            elif shape == (N,):
                y = other.copy()
                self.matrix.matvec(other, y)
                return y
            else:
                raise TypeError
            
    def __rmul__(self, other):
        if type(numerix.ones(1)) == type(other):
            y = other.copy()
            self.matrix.matvec_transp(other, y)
            return y
        else:
            return self * other
            
    def _getShape(self):
        return self.matrix.shape
        
    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)
        
            >>> L = _PysparseMatrix(size=3)
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
        
            >>> L = _PysparseMatrix(size=3)
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
            ids = numerix.arange(self._getShape()[0])
            tmp = numerix.zeros((self._getShape()[0],), 'd')
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
        ids = numerix.arange(self._getShape()[0])
        return self.take(ids, ids)

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)
        
            >>> L = _PysparseMatrix(size=3)
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
            ids = numerix.arange(self._getShape()[0])
            tmp = numerix.zeros((self._getShape()[0],), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)

    def getNumpyArray(self):
        shape = self._getShape()
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
    
class _PysparseMatrix(_PysparseMatrixBase):
    
    """
    _PysparseMatrix class wrapper for pysparse.
    _PysparseMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, size, bandwidth=0, sizeHint=None, matrix=None, storeZeros=True):
        """Creates a `_PysparseMatrix`.

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `storeZeros`: Instructs pysparse to store zero values if possible.
          
        """
        sizeHint = sizeHint or size * bandwidth
        if matrix is None:
            tmpMatrix = spmatrix.ll_mat(1, 1, 1)
            if hasattr(tmpMatrix, 'storeZeros'):
                matrix = spmatrix.ll_mat(size, size, sizeHint, storeZeros)
            else:
                matrix = spmatrix.ll_mat(size, size, sizeHint)
                
        _PysparseMatrixBase.__init__(self, matrix=matrix)

class _PysparseMeshMatrix(_PysparseMatrix):
    
    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, storeZeros=True):
        """Creates a `_PysparseMatrix` associated with a `Mesh`.

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `storeZeros`: Instructs pysparse to store zero values if possible.
          
        """
        self.mesh = mesh
        _PysparseMatrix.__init__(self, size=mesh.getNumberOfCells(), bandwidth=bandwidth, sizeHint=sizeHint, matrix=matrix, storeZeros=storeZeros)
        
    def __mul__(self, other):
        if isinstance(other, _PysparseMeshMatrix):
            return _PysparseMeshMatrix(mesh=self.mesh, 
                                       matrix=spmatrix.matrixmultiply(self.matrix, other._getMatrix()))
        else:
            return _PysparseMatrix.__mul__(self, other)
        
    def asTrilinosMeshMatrix(self):
        """Transforms a pysparse matrix into a trilinos matrix and maintains the
        trilinos matrix as an attribute.
        
        :Returns: 
          The trilinos matrix.

        """

        A = self.matrix.copy()
        values, irow, jcol = A.find()
        
        if not hasattr(self, 'trilinosMatrix'):
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrixKeepStencil
            if A.shape[0] == 0:
                bandwidth = 0
            else:
                bandwidth = int(numerix.ceil(float(len(values)) / float(A.shape[0])))
            self.trilinosMatrix = _TrilinosMeshMatrixKeepStencil(mesh=self.mesh, bandwidth=bandwidth) 

        self.trilinosMatrix.addAt(values, irow, jcol)
        self.trilinosMatrix.finalize()

        return self.trilinosMatrix

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
        
        >>> m = _PysparseMatrix(size=3, storeZeros=True)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(m.matrix.keys(), [(0, 1), (1, 0), (2, 2)])
        True
        >>> print not hasattr(m.matrix, 'storeZeros') or numerix.allequal(m.matrix.values(), [1., 2., 0.]) 
        True
        >>> m = _PysparseMatrix(size=3, storeZeros=False)
        >>> m.addAt((1., 0., 2.), (0, 2, 1), (1, 2, 0))
        >>> print numerix.allequal(m.matrix.keys(), [(0, 1), (1, 0)])
        True
        >>> print numerix.allequal(m.matrix.values(), numerix.array([1.0, 2.0]))
        True
        
        """
        pass
        
class _PysparseIdentityMatrix(_PysparseMatrix):
    """
    Represents a sparse identity matrix for pysparse.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal
        
            >>> print _PysparseIdentityMatrix(size=3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PysparseMatrix.__init__(self, size=size, bandwidth = 1)
        ids = numerix.arange(size)
        self.put(numerix.ones(size, 'd'), ids, ids)
        
class _PysparseIdentityMeshMatrix(_PysparseIdentityMatrix):
    def __init__(self, mesh):
        """
        Create a sparse matrix associated with a `Mesh` with '1' in the diagonal
        
            >>> from fipy import Grid1D
            >>> from fipy.tools import serial
            >>> mesh = Grid1D(nx=3, communicator=serial)
            >>> print _PysparseIdentityMeshMatrix(mesh=mesh)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _PysparseIdentityMatrix.__init__(self, size=mesh.getNumberOfCells())


def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
