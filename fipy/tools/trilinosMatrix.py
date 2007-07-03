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

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.tools.sparseMatrix import _SparseMatrix
from fipy.tools import numerix

# This matrix class now passes the tests. 

# Current inadequacies:

# 1) Adding matrices - the matrix with fewer nonzeros gets added into the one
# that has more; this works as long as it's nonzero entries are a subset of the
# larger one's nonzero entries.  Should be true for all cases in fipy, but is
# not true in the general case.
#
# 2) addAt currently not guaranteed to work for fill-completed matrices, if
# elements are being added in new spots.
#
# 3) put currently not guaranteed to work for non-empty matrices that do not
# have all the target spots occupied. 
#
# None of these situations currently come up in FiPy; tests do not reveal any of 
# the warnings that guard for those. 

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

        """
        if matrix != None:
            self.matrix = matrix
        else:
            if sizeHint is not None and bandwidth == 0:
                bandwidth = (sizeHint + size - 1)/size 
            self.comm = Epetra.PyComm()
            if self.comm.NumProc() == 1:
                self.parallel = False
                self.map = Epetra.Map(size, 0, self.comm)
            else:
                self.parallel = True
                self.startRow = self.comm.MyPID()*size/self.comm.NumProc()
                self.endRow = (self.comm.MyPID()+1)*size/self.comm.NumProc()
                self.map = Epetra.Map(size, range(self.startRow, self.endRow), 0, self.comm)

            self.matrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, bandwidth*3/2)
            # Leave extra bandwidth, to handle multiple insertions into the
            # same spot. It's memory-inefficient, but it'll get cleaned up when
            # FillComplete is called, and according to the Trilinos devs the 
            # performance boost will be worth it.

    def _getMatrix(self):
        return self.matrix
    
    # All operations that require getting data out of the matrix need to muddle
    # around with FillComplete to make sure they work.  There will be no
    # warnings when FillComplete is implicitly called; there will only be
    # warnings when insertions fail.
    def copy(self):
        if not self._getMatrix().Filled():
            self._getMatrix().FillComplete()

        return _TrilinosMatrix(matrix = Epetra.FECrsMatrix(self.matrix))
            
        
    def __getitem__(self, index):
        if not self.matrix.Filled():
            self._getMatrix().FillComplete()

        return self.matrix[index]
        
    def __str__(self):
        if not self.matrix.Filled():
            self.matrix.FillComplete()
        return _SparseMatrix.__str__(self)

    def __setitem__(self, index, value):
        self.matrix[index] = value
        

    # Addition is tricky. 
    # Trilinos interface is as such: A can be added into B, but A has to be
    # Filled() beforehand. If B is filled beforehand, this may or may not
    # crash, depending on whether things are being added into spots in B that
    # were not there before.  Have put in some order-of-operands twiddling to
    # make it look like two things can be added in any order.

    # Though not guaranteed to work for arbitrary matrices, it should work for
    # all those generated by FiPy and will give warnings if it encounters
    # trouble (unless Trilinos runs into an error and aborts instead of
    # returning an error code)

    # With the changes to Trilinos that should be in the next version (Add no
    # longer FillComplete()s the destination matrix) it is possible to make
    # this work in the general case by adding both matrices into an empty
    # matrix. This has not yet been implemented. It is not strictly necessary
    # for FiPy currently, since all tests pass without it, but will be a good
    # idea to have so that the matrix class can be used more generally without
    # worrying so much about the underlying trilinos.

    def __iadd__(self, other):
        if other != 0:
            if not other._getMatrix().Filled():
                other._getMatrix().FillComplete()
            
            # Depending on which one is more filled, pick the order of operations 
            if self._getMatrix().Filled() and other._getMatrix().NumGlobalNonzeros() \
                                            > self._getMatrix().NumGlobalNonzeros():
                tempBandwidth = other._getMatrix().NumGlobalNonzeros() \
                                 /self._getMatrix().NumGlobalRows()+1
                tempMatrix = Epetra.FECrsMatrix(Epetra.Copy, self.map, tempBandwidth)
                
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
            self._getMatrix().FillComplete()
            
        if not other._getMatrix().Filled():
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
                    self._getMatrix().FillComplete()
                    
                if not other._getMatrix().Filled():
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
        
    def localizeToProcessor(self, row, col, val):
        if(self.parallel):
            tuples = zip(row, col, val)
            filtered = filter(lambda x: self.startRow <= x[0] and x[0] < self.endRow , tuples)
            if filtered is not None and len(filtered) > 0:
                row, col, val = zip( *filtered)
            else:
                row, col, val = (), (), ()
            return row, col, val
        else:
            return row, col, val
            
            
        
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

        id1, id2, vector = self.localizeToProcessor(id1, id2, vector)
        
        if self._getMatrix().Filled():
            if self._getMatrix().ReplaceGlobalValues(id1, id2, vector) != 0:
                import warnings
                warnings.warn("ReplaceGlobalValues returned error code in put", 
                               UserWarning, stacklevel=2)
                # Possible different algorithm, to guarantee success:
                # 
                # Make a new matrix, 
                # Use addAt to put the values in it, 
                # Use replaceGlobalValues in the original matrix to zero out the terms 
                # And add the old one into the new one, 
                # Replace the old one.
                #
                # Would incur performance costs, and since FiPy does not use 
                # this function in such a way as would generate these errors,
                # I have not implemented the change.

        else:

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
                    # Possible different algorithm, to guarantee that it does not fail:
                    # 
                    # Make a new matrix, 
                    # Use addAt to put the values in it, 
                    # Use replaceGlobalValues in the original matrix to zero out the terms 
                    # And add the old one into the new one, 
                    # Replace the old one.
                    #
                    # Would incur performance costs, and since FiPy does not use 
                    # this function in such a way as would generate these errors,
                    # I have not implemented the change.
            
                             


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

    def takeDiagonal(self):
        if not self._getMatrix().Filled():
            self._getMatrix().FillComplete()

        result = Epetra.Vector(self.map)
        self._getMatrix().ExtractDiagonalCopy(result)
        return numerix.array(result)
    
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
        id1, id2, vector = self.localizeToProcessor(id1, id2, vector)
        if not self._getMatrix().Filled():
            self._getMatrix().InsertGlobalValues(id1, id2, vector)
        else:
            if self._getMatrix().SumIntoGlobalValues(id1, id2, vector) != 0:
                import warnings
                warnings.warn("Summing into filled matrix returned error code",
                               UserWarning, stacklevel=2)
                # Possible change to this part of the code to do the following:
                #
                # Make a new matrix, 
                # Use addAt to put the values in it
                # Add the old one into the new one
                # Replace the old one. 
                #
                # Would incur performance costs, and since FiPy does not use 
                # this function in such a way as would generate these errors,
                # I have not implemented the change.


    def addAtDiagonal(self, vector):
        if type(vector) in [type(1), type(1.)]:
            ids = numerix.arange(self._getMatrix().GetGlobalRows())
            tmp = numerix.zeros((self._getMatrix().GetGlobalRows(),), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)


    # This function requires take. This is a problem. 
    # However, this does not seem to ever be called.

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
