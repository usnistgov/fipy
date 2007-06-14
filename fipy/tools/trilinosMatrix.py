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

from fipy.tools.sparseMatrix import _SparseMatrix
from fipy.tools import numerix

# List of current hacks - things that work in special cases, but 
# will not work in general.
# 1) Adding matrices - the matrix with fewer nonzeros gets added
# into the one that has more; this works as long as it's nonzero
# entries are a subset of the larger one's nonzero entries.
# Should be true for all cases in fipy.

class _TrilinosMatrix(_SparseMatrix):
    
    """
    _TrilinosMatrix class wrapper for a PyTrilinos Epetra.FECrsMatrix.
    _TrilinosMatrix is always NxN.
    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """

    def __init__(self, size = None, bandwidth = 0, matrix = None, sizeHint = None):
        """
        Creates a `_TrilinosMatrix`.

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
            self.matrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, 2*bandwidth)
            # Leave twice the bandwidth, to handle multiple insertions into the
            # same spot. It's memory-inefficient, but it'll get cleaned up when
            # FillComplete is called, and according to the Trilinos devs the 
            # performance boost will be worth it.

    def _getMatrix(self):
        return self.matrix
    
    # All operations that require getting data out of the matrix need to muddle around
    # with FillComplete to make sure they work. 
    # There will be no warnings when FillComplete is implicitly called; there will only
    # be warnings when insertions fail.
    def copy(self):
        if not self._getMatrix().Filled():
            #import warnings
            #warnings.warn("""Matrix should have FillComplete called on it before being 
            #                 copied. FillComplete will now be called on this matrix.""",
            #                 UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()

        return _TrilinosMatrix(matrix = Epetra.FECrsMatrix(self.matrix))
            
        
    def __getitem__(self, index):
        if not self.matrix.Filled():
            #import warnings
            #warnings.warn("""Matrix should have FillComplete called on it before being 
            #                 accessed. FillComplete will be called on the matrix.""",
            #                 UserWarning, stacklevel=2)
            self._getMatrix().FillComplete()
        return self.matrix[index]
        
    def __str__(self):
        if not self.matrix.Filled():
            self.matrix.FillComplete()
        return _SparseMatrix.__str__(self)

    def __setitem__(self, index, value):
        self.matrix[index] = value
        

    # Addition is tricky. 
    # Trilinos interface is as such: A can be added into B, where B is a
    # matrix, but this will automatically FillComplete() B; A has to be
    # Filled() beforehand.  If B is filled beforehand, this may or may not
    # crash, depending on whether things are being added into spots in B that
    # were not there before.  Have put in some order-of-operands twiddling 
    # to make it look like two things can be added in any order. 

    # Addition not guaranteed to work for arbitrary matrices, but should
    # work for all those generated by FiPy and will give warnings if it 
    # encounters trouble.

    def __iadd__(self, other):
        if other != 0:
            if not other._getMatrix().Filled():
                #import warnings
                #warnings.warn("Matrix must be Filled() before being added",
                #                 UserWarning, stacklevel=2)
                other._getMatrix().FillComplete()
            
            # Depending on which one is more filled, pick the order of operations 
            if self._getMatrix().Filled() and other._getMatrix().NumGlobalNonzeros() \
                                            > self._getMatrix().NumGlobalNonzeros():
                tempMatrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, (other._getMatrix().NumGlobalNonzeros()/self._getMatrix().NumGlobalRows())+1)
                if EpetraExt.Add(other._getMatrix(), False, 1, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__ 1",
                                   UserWarning, stacklevel=2)

                if EpetraExt.Add(self._getMatrix(), False, 1, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__ 2",
                                   UserWarning, stacklevel=2)

                self.matrix = tempMatrix
            else:
                if EpetraExt.Add(other._getMatrix(), False,1,self._getMatrix(),1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__",
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
        # so addition is likely to succeed
        if self._getMatrix().NumGlobalNonzeros() > other._getMatrix().NumGlobalNonzeros():
            tempMatrix = self.copy()
            tempMatrix.__iadd__(other*sign)
        else:
            tempMatrix = other.copy()
            tempMatrix.__iadd__(self*sign)
            
        return tempMatrix

    def __add__(self, other):
        """
        Add two sparse matrices. The nonempty spots of one of them must be a 
        subset of the nonempty spots of the other one.
        
            >>> L = _TrilinosMatrix(size = 3)
            >>> L.addAt((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt([0,0,0], [0,1,2], [0,1,2])
            >>> print L + _TrilinosIdentityMatrix(3)
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
        
    def __sub__(self, other):
        if other is 0:
            return self
        else:
            return self._add(other, sign=-1)

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix
        
            >>> L1 = _TrilinosMatrix(size = 3)
            >>> L1.addAt((3,10,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L2 = _TrilinosIdentityMatrix(size = 3)
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

        if isinstance(other, _TrilinosMatrix):
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
                return _TrilinosMatrix(matrix = result)
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
            
    def _getShape(self):
        N = self._getMatrix().NumGlobalCols()
        return (N,N)
        
    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)
        
            >>> L = _TrilinosMatrix(size = 3)
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
                self._getMatrix().InsertGlobalValues(id1, id2, numerix.zeros(len(vector)))
                self._getMatrix().FillComplete()
                if self._getMatrix().ReplaceGlobalValues(id1, id2, vector) != 0:
                    import warnings
                    warnings.warn("ReplaceGlobalValues returned error code in put", 
                                   UserWarning, stacklevel=2)
            
                             


    def putDiagonal(self, vector):
        """
        Put elements of `vector` along diagonal of matrix
        
            >>> L = _TrilinosMatrix(size = 3)
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
        
            >>> L = _TrilinosMatrix(size = 3)
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

    def getNumpyArray(self):
        raise NotImplemented
        
class _TrilinosIdentityMatrix(_TrilinosMatrix):
    """
    Represents a sparse identity matrix for Trilinos.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal
        
            >>> print _TrilinosIdentityMatrix(size = 3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _TrilinosMatrix.__init__(self, size = size, bandwidth = 1)
        ids = numerix.arange(size)
        self.addAt(numerix.ones(size), ids, ids)
        
def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
