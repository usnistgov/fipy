#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "trilinosMatrix.py"
 #                                    created: 06/08/07
 #                                last update: 06/11/07
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Maxsim Gibiansky <maxsim.gibiansky@nist.gov>
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2007-06-11 MLG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

#EpetraExt necessary for matrix operations (Add, Multiply)
from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.tools import numerix


# List of current hacks - things that work in special cases, but 
# will not work in general.
# 1) Adding matrices - the matrix with fewer nonzeros gets added
# into the one that has more; this works as long as it's nonzero
# entries are a subset of the larger one's nonzero entries.
# Should be true for all cases in fipy.
# 2) put - not guaranteed to work on not-filled matrices. It does
# a replaceGlobalValues to replace them with 0, and then inserts 
# the new values, but this could go wrong if there have already been
# things "put" there. Fortunately, it's not used much. Maybe start 
# by FillCompleteing the matrix?
class _SparseMatrix:
    
    """
    _SparseMatrix class wrapper for a PyTrilinos Epetra.FECrsMatrix.
    _SparseMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, size = None, bandwidth = 0, matrix = None, sizeHint = None):
        """
        Creates a `_SparseMatrix`.

        :Parameters:
          - `size`: The size N for an N by N matrix.
          - `bandwidth`: The proposed band width of the matrix.
          - `matrix`: The starting `spmatrix` id there is one.
          -

        """
        if matrix != None:
            self.matrix = matrix
        else:
            if sizeHint is not None and bandwidth == 0:
                bandwidth = (sizeHint + size - 1)/size 
            self.comm = Epetra.PyComm()
            self.map = Epetra.Map(size, 0, self.comm)
            self.matrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, bandwidth)

    __array_priority__ = 100.0    

    def __array_wrap__(self, arr, context=None):
        if context is None:
            return arr
        else:
            return NotImplemented

    def _getMatrix(self):
        return self.matrix
    
    # All operations that require getting data out of the matrix need to muddle around
    # with FillComplete to make sure they work. Warnings will mark the spot. 
    def copy(self):
        if self._getMatrix().Filled():
            return _SparseMatrix(matrix = Epetra.FECrsMatrix(self.matrix))
        else: 
            #warnings.warn("""Matrix should have FillComplete called on it before being 
            #                 copied. FillComplete will now be called on this matrix.""",
            #                 UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()
            return _SparseMatrix(matrix = Epetra.FECrsMatrix(self.matrix))
            
        
    def __getitem__(self, index):
        if self.matrix.Filled():
            return self.matrix[index]
        else:
            #import warnings
            #warnings.warn("""Matrix should have FillComplete called on it before being 
            #                 accessed. FillComplete will be called on the matrix.""",
            #                 UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()
            return self.matrix[index]
        
    def __str__(self):
        if not self.matrix.Filled():
            self.matrix.FillComplete()
        s = ""
        cellWidth = 11
        shape = self._getShape()
        for j in range(shape[1]):
            for i in range(shape[0]):
                v = self[j,i]
                if v == 0:
                    s += "---".center(cellWidth)
                else:
                    exp = numerix.log(abs(v))
                    if abs(exp) <= 4:
                        if exp < 0:
                            s += ("%9.6f" % v).ljust(cellWidth)
                        else:
                            s += ("%9.*f" % (6,v)).ljust(cellWidth)
                    else:
                        s += ("%9.2e" % v).ljust(cellWidth)
            s += "\n"
        return s[:-1]
            
    def __repr__(self):
 	    return self.matrix.__repr__()
        
    def __setitem__(self, index, value):
        self.matrix[index] = value
        

    # Addition is tricky. 
    # Trilinos interface is as such: A can be added into B, where B is a
    # matrix, but this will automatically FillComplete() B; A has to be
    # Filled() beforehand.  If B is filled beforehand, this may or may not
    # crash, depending on whether things are being added into spots in B that
    # were not there before.  Have to really finesse it to make it look like
    # two things can (sometimes) be just added. Can be aided by changing order
    # of operands.

    def _iadd(self, L, other, sign = 1):
        if other != 0:
            if not other._getMatrix().Filled():
                #import warnings
                #warnings.warn("Matrix must be Filled() before being added",
                #                 UserWarning, stacklevel=2)
                other._getMatrix().FillComplete()
            
            # Depending on which one is more filled, pick the order of operations 
            if L.Filled() and other._getMatrix().NumGlobalNonzeros() \
                              > self._getMatrix().NumGlobalNonzeros():
                tempMatrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, (other._getMatrix().NumGlobalNonzeros()/L.NumGlobalRows())+1)
                if EpetraExt.Add(other._getMatrix(), False, sign, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in _iadd 1",
                                   UserWarning, stacklevel=2)

                if EpetraExt.Add(L, False, 1, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in _iadd 2",
                                   UserWarning, stacklevel=2)

                L = tempMatrix
            else:
                if EpetraExt.Add(other._getMatrix(), False,sign,L,1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in _iadd",
                                   UserWarning, stacklevel=2)
        return self

   
    # To add two things while modifying neither, both must be FillCompleted
    def _add(self, other, sign = 1):
        if not self._getMatrix().Filled():
            #import warnings
            #warnings.warn("Matrix must be FillComplete()d before being added",
            #               UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()
        if not other._getMatrix().Filled():
            #import warnings
            #warnings.warn("Matrix must be FillComplete()d before being added",
            #               UserWarning, stacklevel=2)
            other._getMatrix().FillComplete()
        
        # make the one with more nonzeros the right-hand operand
        # is likely to succeed
        if self._getMatrix().NumGlobalNonzeros() > other._getMatrix().NumGlobalNonzeros():
            L = Epetra.FECrsMatrix(self._getMatrix())
            other._iadd(L, other, sign)
        else:
            L = Epetra.FECrsMatrix(other._getMatrix())
            self._iadd(L, self, sign)
            
        return _SparseMatrix(matrix = L)

    def __add__(self, other):
        """
        Add two sparse matrices. The nonempty spots of one of them must be a 
        subset of the nonempty spots of the other one.
        
            >>> L = _SparseMatrix(size = 3)
            >>> L.addAt((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt([0,0,0], [0,1,2], [0,1,2])
            >>> print L + _SparseIdentityMatrix(3)
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
            return self._add(other)
        
    __radd__ = __add__

    def __sub__(self, other):
        if other is 0:
            return self
        else:
            return self._add(other, sign=-1)

    def __rsub__(self, other):
        if other is 0:
            return -self
        else:
            return other.add(self, sign=-1)
    # Check this, its different than the pysparse wrapper, pysparse didn't have the else
        
    def __iadd__(self, other):
        return self._iadd(self._getMatrix(), other)
        
    def __isub__(self, other):
        return self._iadd(self._getMatrix(), other, -1)
    # Check whether this should be _isub_ or __isub__ or _isub__
    # it's _isub__ in the pysparse wrapper. Does it matter? Is this ever used?
        
    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = _SparseMatrix(size = 3)
            >>> L1.addAt((3,10,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L2 = _SparseIdentityMatrix(size = 3)
            >>> L2.addAt((4.38,12357.2,1.1), (2,1,0), (1,0,2))
            
            >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
            ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
            ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

            >>> for i in range(0,3):
            ...     for j in range(0,3):
            ...         numerix.allclose(((L1*L2)[i,j],), tmp[i,j])
            True
            True
            True
            True
            True
            True
            True
            True
            True

        or a sparse matrix by a vector

            >>> tmp = numerix.array((29., 6.28318531, 2.5))       
            >>> numerix.allclose(L1 * numerix.array((1,2,3),'d'), tmp)
            1
            
        or a vector by a sparse matrix

            >>> tmp = numerix.array((7.5, 16.28318531,  3.))  
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp) 
            1

            
        """
        N = self._getMatrix().NumGlobalCols()

        if isinstance(other, _SparseMatrix):
            if isinstance(other._getMatrix(), Epetra.RowMatrix):
                if not self._getMatrix().Filled():
                    #import warnings
                    #warnings.warn("Matrix must be FillComplete()d before being multiplied",
                                   #UserWarning, stacklevel=2)
                    self._getMatrix().FillComplete()
                if not other._getMatrix().Filled():
                    #import warnings
                    #warnings.warn("Matrix must be FillComplete()d before being multiplied",
                                   #UserWarning, stacklevel=2)
                    other._getMatrix().FillComplete()

                result = Epetra.FECrsMatrix(Epetra.Copy, self.map, 0)

                EpetraExt.Multiply(self._getMatrix(), False, other._getMatrix(), False, result)
                return _SparseMatrix(matrix = result)
            else:
                raise TypeError
                
        else:
            shape = numerix.shape(other)
            if shape == ():
                result = self.copy()
                result._getMatrix().Scale(other)
                return result
            elif shape == (N,):
                if not self._getMatrix().Filled():
                    #import warnings
                    #warnings.warn("Matrix must be FillComplete()d before being multiplied",
                                   #UserWarning, stacklevel=2)
                    self._getMatrix().FillComplete()

                y = Epetra.Vector(other)
                result = Epetra.Vector(self.map)
                self._getMatrix().Multiply(False, y, result)
                return result
            else:
                raise TypeError
           
    def __rmul__(self, other):
        if type(numerix.ones(1)) == type(other):
            y = Epetra.Vector(other)
            result = Epetra.Vector(self.map)
            self._getMatrix().Multiply(True, y, result)
            return result
        else:
            return self * other
            
    def __neg__(self):
        """
        Negate a sparse matrix
        
            >>> print -_SparseIdentityMatrix(size = 3)
            -1.000000      ---        ---    
                ---    -1.000000      ---    
                ---        ---    -1.000000  
        """
        return self * -1
        
    def __pos__(self):
        return self
        
##     def __eq__(self,other):
## 	return self.matrix.__eq__(other._getMatrix())

    def _getShape(self):
        N = self._getMatrix().NumGlobalCols()
        return (N,N)
        
##     def transpose(self):
##         pass

    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)
        
            >>> L = _SparseMatrix(size = 3)
            >>> L.put((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> print L
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """
        if self._getMatrix().Filled():
            if self._getMatrix().ReplaceGlobalValues(id1, id2, vector) != 0:
                import warnings
                warnings.warn("ReplaceGlobalValues returned error code in put", 
                               UserWarning, stacklevel=2)

        else:
            #import warnings
            #warnings.warn("""Putting into an unfilled matrix can silently fail!""",
            #UserWarning, stackLevel=2)

            # This guarantees that it will actually replace the values that are there,
            # if there are any
            if self._getMatrix().NumGlobalNonzeros() == 0:
                self._getMatrix().InsertGlobalValues(id1, id2, vector)
            else:
                self._getMatrix().FillComplete()
                if self._getMatrix().ReplaceGlobalValues(id1, id2, vector) != 0:
                    import warnings
                    warnings.warn("ReplaceGlobalValues returned error code in put", 
                                   UserWarning, stacklevel=2)
            
                             


    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix
        
            >>> L = _SparseMatrix(size = 3)
            >>> L.putDiagonal((3.,10.,numerix.pi))
            >>> print L
             3.000000      ---        ---    
                ---    10.000000      ---    
                ---        ---     3.141593  
            >>> L.putDiagonal((10.,3.))
            >>> print L
            10.000000      ---        ---    
                ---     3.000000      ---    
                ---        ---     3.141593  
        """
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._getMatrix().NumGlobalRows())
            tmp = numerix.zeros((self._getMatrix().NumGlobalRows), 'd')
            tmp[:] = vector
            self.put(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.put(vector, ids, ids)

    def take(self, id1, id2):
        import warnings
        warnings.warn("""Trying to take from a Trilinos Matrix. That doesn't work.""",
                         UserWarning, stacklevel=2)
        raise TypeError
       #vector = numerix.zeros(len(id1), 'd')
       #self.matrix.take(vector, id1, id2)
       #return vector

    def takeDiagonal(self):
        if not self._getMatrix().Filled():
            #import warnings
            #warnings.warn("""Matrix should have FillComplete called on it before being 
                             #read. FillComplete will now be called on this matrix.""",
                             #UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()

        result = Epetra.Vector(self.map)
        self._getMatrix().ExtractDiagonalCopy(result)
        return result[:]
    
    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)
        
            >>> L = _SparseMatrix(size = 3)
            >>> L.addAt((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt((1.73,2.2,8.4,3.9,1.23), (1,2,0,0,1), (2,2,0,0,2))
            >>> print L
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """
        if not self._getMatrix().Filled():
            self._getMatrix().InsertGlobalValues(id1, id2, vector)
        else:
            #import warnings
            #warnings.warn("Adding into filled matrix", UserWarning, stacklevel=2)
            if self._getMatrix().SumIntoGlobalValues(id1, id2, vector) != 0:
                import warnings
                warnings.warn("Summing into filled matrix returned error code",
                               UserWarning, stacklevel=2)

    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._getMatrix().GetGlobalRows())
            tmp = numerix.zeros((self._getMatrix().GetGlobalRows(),), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)


    # These functions require take. This is a problem. 
    # If these are used a bunch, I might have to put in a very inefficient, 
    # python-loopy take.    

   #def getNumpyArray(self):
   #    shape = self._getShape()
   #    indices = numerix.indices(shape)
   #    numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
   #    return numerix.reshape(numMatrix, shape)
        
##     def __array__(self):
##      shape = self._getShape()
##      indices = numerix.indices(shape)
##         numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
##      return numerix.reshape(numMatrix, shape)

#   def matvec(self, x)
#       """
#       This method is required for scipy solvers.
#       """
#       return self * x
    

class _SparseIdentityMatrix(_SparseMatrix):
    """
    Represents a sparse identity matrix.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal
        
            >>> print _SparseIdentityMatrix(size = 3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _SparseMatrix.__init__(self, size = size, bandwidth = 1)
        ids = numerix.arange(size)
        self.addAt(numerix.ones(size), ids, ids)
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
