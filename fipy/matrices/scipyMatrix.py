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

import scipy.sparse as sp
from fipy.tools import numerix

from fipy.matrices.sparseMatrix import _SparseMatrix

class _ScipyMatrixBase(_SparseMatrix):
    
    """
    _PysparseMatrix class wrapper for pysparse.
    _PysparseMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, matrix):
        """Creates a `_ScipyMatrixBase`.

        :Parameters:
          - `matrix`: The starting `spmatrix` 
        """
        self.matrix = matrix
   
    @property
    def _shape(self):
        return self.matrix.shape

    def _getRange(self):
        return range(self._shape[1]), range(self._shape[0])
        
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
        assert(len(id1) == len(id2) == len(vector))

        for v, i1, i2 in zip(vector, id1, id2):
            self.matrix[i1, i2] = v

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
                self.matrix.diagonal()
        """
        if type(vector) in [int, float]:
            vector = numerix.repeat(vector, self._shape[0])

        self.matrix.setdiag(vector)
            
    def take(self, id1, id2):
        assert(len(id1) == len(id2))

        return self.matrix[zip(id1, id2)] # fancy indexing ftw

    def takeDiagonal(self):
        self.matrix.diagonal()

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
        assert(len(id1) == len(id2) == len(vector))

        for v, i1, i2 in zip(vector, id1, id2):
            self.matrix[i1, i2] += v

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

class _PysparseMatrix(_ScipyMatrixBase):
    
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
                
        _ScipyMatrixBase.__init__(self, matrix=matrix)

class _PysparseMeshMatrix(_PysparseMatrix):
    
    def __init__(self, mesh, bandwidth=0, sizeHint=None, matrix=None, numberOfVariables=1, storeZeros=True):

        """Creates a `_PysparseMatrix` associated with a `Mesh`.

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `numberOfVariables`: The size of the matrix is determined by numberOfVariables * self.mesh.numberOfCells.
          - `storeZeros`: Instructs pysparse to store zero values if possible.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        _PysparseMatrix.__init__(self, size=self.numberOfVariables * self.mesh.numberOfCells, bandwidth=bandwidth, sizeHint=sizeHint, matrix=matrix, storeZeros=storeZeros)

    def __mul__(self, other):
        if isinstance(other, _PysparseMeshMatrix):
            return _PysparseMeshMatrix(mesh=self.mesh, 
                                       matrix=spmatrix.matrixmultiply(self.matrix, other.matrix))
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
            if A.shape[0] == 0:
                bandwidth = 0
            else:
                bandwidth = int(numerix.ceil(float(len(values)) / float(A.shape[0])))
            from fipy.matrices.trilinosMatrix import _TrilinosMeshMatrixKeepStencil
            self.trilinosMatrix = _TrilinosMeshMatrixKeepStencil(mesh=self.mesh, bandwidth=bandwidth, numberOfVariables=self.numberOfVariables)

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
        _PysparseIdentityMatrix.__init__(self, size=mesh.numberOfCells)

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test()   
