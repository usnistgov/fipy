#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "trilinosMatrix.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from PyTrilinos import Epetra
from PyTrilinos import EpetraExt

from fipy.matrices.sparseMatrix import _SparseMatrix
from fipy.tools import numerix, parallelComm

# Current inadequacies of the matrix class:

# 1) Adding matrices - the matrix with fewer nonzeros gets added into the one
# that has more; this works as long as it's nonzero entries are a subset of the
# larger one's nonzero entries. Is true for all cases in fipy, but is not true
# in the general case - this isn't a general matrix class like the pysparse
# matrix class is.
#
# 2) addAt currently not guaranteed to work for fill-completed matrices, if
# elements are being added in new spots.
#
# 3) put currently not guaranteed to work for non-empty matrices that do not
# have all the target spots occupied.
#
# None of these situations currently come up in FiPy; tests do not reveal any of
# the warnings that guard for those, and all tests pass. Because of the way
# FiPy constructs its matrices, I do not anticipate any of these occurring.

class _TrilinosMatrix(_SparseMatrix):
    """class wrapper for a PyTrilinos Epetra.CrsMatrix.

    Allows basic python operations __add__, __sub__ etc.
    Facilitate matrix populating in an easy way.
    """
    def __init__(self, matrix, bandwidth=None,
                 rowMap=None, colMap=None, domainMap=None):
        """
        :Parameters:
          - `matrix`: The starting `Epetra.CrsMatrix` if there is one.
          - `bandwidth`: The proposed band width of the matrix.
        """
        self.matrix = matrix

        self.rowMap = rowMap or matrix.RowMap()
        self.colMap = colMap or self.rowMap
        self.domainMap = domainMap or self.colMap
        self.rangeMap = self.rowMap

        self.comm = matrix.Comm()
        if bandwidth is None:
            self.bandwidth = ((matrix.NumGlobalNonzeros() + matrix.NumGlobalRows() -1 )
                              // matrix.NumGlobalRows())
        else:
            self.bandwidth = bandwidth

    def _setMatrix(self, m):
        self._matrix = m

    matrix = property(lambda self: self._matrix, _setMatrix)

    # All operations that require getting data out of the matrix may need to
    # call FillComplete to make sure they work.  There will be no warnings when
    # FillComplete is implicitly called; there will only be warnings when
    # insertions fail.
    def copy(self):
        self.fillComplete()

        return _TrilinosMatrix(matrix=Epetra.CrsMatrix(self.matrix))


    def __getitem__(self, index):
        self.fillComplete()

        return self.matrix[index]

    def __str__(self):
        self.fillComplete()

        s = _SparseMatrix.__str__(self)

        comm = self.matrix.Map().Comm()
        if comm.NumProc() > 1:
            from fipy.tools import parallelComm
            return ''.join(parallelComm.allgather(s))
        else:
            return s

    @property
    def _range(self):
        return (range(self.rowMap.NumGlobalElements()), self.rowMap.MyGlobalElements())

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

    def __iadd__(self, other):
        if other != 0:
            other.fillComplete()

            # Depending on which one is more filled, pick the order of operations
            if self.matrix.Filled() and other.matrix.NumGlobalNonzeros() \
                                            > self.matrix.NumGlobalNonzeros():
                tempBandwidth = other.matrix.NumGlobalNonzeros() \
                                 /self.matrix.NumGlobalRows()+1

                tempMatrix = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, tempBandwidth)

                if EpetraExt.Add(other.matrix, False, 1, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__, 1",
                                   UserWarning, stacklevel=2)

                if EpetraExt.Add(self.matrix, False, 1, tempMatrix, 1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__, 2",
                                   UserWarning, stacklevel=2)

                self.matrix = tempMatrix

            else:
                if EpetraExt.Add(other.matrix, False,1,self.matrix,1) != 0:
                    import warnings
                    warnings.warn("EpetraExt.Add returned error code in __iadd__",
                                   UserWarning, stacklevel=2)

        return self


    # To add two things while modifying neither, both must be FillCompleted
    def _add(self, other, sign = 1):
        self.fillComplete()
        other.fillComplete()

        # make the one with more nonzeros the right-hand operand
        # so addition is likely to succeed
        if self.matrix.NumGlobalNonzeros() > other.matrix.NumGlobalNonzeros():
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

            >>> L = _TrilinosMatrixFromShape(rows=3, cols=3)
            >>> L.addAt((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt([0,0,0], [0,1,2], [0,1,2])
            >>> print L + _TrilinosIdentityMatrix(size=3)
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
            AttributeError: 'int' object has no attribute 'fillComplete'
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

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix.

        >>> L1 = _TrilinosMatrixFromShape(rows=3, cols=3)
        >>> L1.addAt((3,10,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
        >>> L2 = _TrilinosIdentityMatrix(size=3)
        >>> L2.addAt((4.38,12357.2,1.1), (2,1,0), (1,0,2))

        >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
        ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
        ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

        >>> L = (L1 * L2).numpyArray

        >>> print numerix.allclose(tmp, L) # doctest: +SERIAL
        True

        or a sparse matrix by a vector

        >>> tmp = numerix.array((29., 6.28318531, 2.5))
        >>> print numerix.allclose(L1 * numerix.array((1,2,3),'d'), tmp) # doctest: +SERIAL
        True

        or a vector by a sparse matrix

        >>> tmp = numerix.array((7.5, 16.28318531,  3.))
        >>> print numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp)  # doctest: +SERIAL
        True


        """
        N = self.matrix.NumMyCols()

        if isinstance(other, _TrilinosMatrix):
            if isinstance(other.matrix, Epetra.RowMatrix):
                self.fillComplete()
                other.fillComplete()

                result = Epetra.CrsMatrix(Epetra.Copy, self.rowMap, 0)

                EpetraExt.Multiply(self.matrix, False, other.matrix, False, result)
                copy = self.copy()
                copy.matrix = result
                return copy
            else:
                raise TypeError

        else:
            shape = numerix.shape(other)
            if shape == ():
                result = self.copy()
                result.matrix.Scale(other)
                return result
            elif shape == (N,):
                self.fillComplete()

                y = Epetra.Vector(self.domainMap, other)
                result = Epetra.Vector(self.rangeMap)
                self.matrix.Multiply(False, y, result)
                return numerix.array(result)
            else:
                raise TypeError

    def __rmul__(self, other):
        if type(numerix.ones(1, 'l')) == type(other):
            self.fillComplete()

            y = Epetra.Vector(self.rangeMap, other)
            result = Epetra.Vector(self.domainMap)
            self.matrix.Multiply(True, y, result)
            return numerix.array(result)
        else:
            return self * other

    @property
    def _shape(self):
        N = self.matrix.NumGlobalRows()
        return (N,N)



    def put(self, vector, id1, id2):
        """
        Put elements of `vector` at positions of the matrix corresponding to (`id1`, `id2`)

            >>> L = _TrilinosMatrixFromShape(rows=3, cols=3)
            >>> L.put((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> print L
                ---    10.000000   3.000000  
                ---     3.141593      ---    
             2.500000      ---        ---    
        """

        if hasattr(id1, 'dtype') and id1.dtype.name == 'int64':
            id1 = id1.astype('int32')
        if hasattr(id2, 'dtype') and id2.dtype.name == 'int64':
            id2 = id2.astype('int32')

        if self.matrix.Filled():
            if self.matrix.ReplaceGlobalValues(id1, id2, vector) != 0:
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
            if self.matrix.NumGlobalNonzeros() == 0:
                self.matrix.InsertGlobalValues(id1, id2, vector)
            else:
                self.matrix.InsertGlobalValues(id1, id2, numerix.zeros(len(vector), 'l'))
                self.fillComplete()
                if self.matrix.ReplaceGlobalValues(id1, id2, vector) != 0:
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

            >>> L = _TrilinosMatrixFromShape(rows=3, cols=3)
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
            ids = numerix.arange(self.matrix.NumGlobalRows())
            tmp = numerix.zeros((self.matrix.NumGlobalRows), 'd')
            tmp[:] = vector
            if ids.dtype.name == 'int64':
                ids = ids.astype('int32')
            self.put(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            if ids.dtype.name == 'int64':
                ids = ids.astype('int32')
            self.put(vector, ids, ids)

    def take(self, id1, id2):
        import warnings
        warnings.warn("""Trying to take from a Trilinos Matrix. That doesn't work.""",
                         UserWarning, stacklevel=2)
        raise TypeError

    def takeDiagonal(self):
        self.fillComplete()

        result = Epetra.Vector(self.rangeMap)
        self.matrix.ExtractDiagonalCopy(result)

        return result

    def addAt(self, vector, id1, id2):
        """
        Add elements of `vector` to the positions in the matrix corresponding to (`id1`,`id2`)

            >>> L = _TrilinosMatrixFromShape(rows=3, cols=3)
            >>> L.addAt((3.,10.,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L.addAt((1.73,2.2,8.4,3.9,1.23), (1,2,0,0,1), (2,2,0,0,2))
            >>> print L
            12.300000  10.000000   3.000000  
                ---     3.141593   2.960000  
             2.500000      ---     2.200000  
        """

        ## This was added as it seems that trilinos does not like int64 arrays
        if hasattr(id1, 'astype') and id1.dtype.name == 'int64':
            id1 = id1.astype('int32')
        if hasattr(id2, 'astype') and id2.dtype.name == 'int64':
            id2 = id2.astype('int32')

        if not self.matrix.Filled():
            err = self.matrix.InsertGlobalValues(id1, id2, vector)
            if err < 0:
                raise RuntimeError, "Processor %d, error code %d" \
                  % (self.comm.MyPID(), err)
        else:
            if self.matrix.SumIntoGlobalValues(id1, id2, vector) != 0:
                import warnings
                warnings.warn("Summing into unfilled matrix returned error code",
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

            if hasattr(self.matrix, 'GetMyRows'):
                Nrows = self.matrix.GetMyRows()
            else:
                Nrows = self.matrix.NumMyRows()

            ids = numerix.arange(Nrows)
            tmp = numerix.zeros((Nrows,), 'd')
            tmp[:] = vector
            self.addAt(tmp, ids, ids)
        else:
            ids = numerix.arange(len(vector))
            self.addAt(vector, ids, ids)

    def exportMmf(self, filename):
        """
        Exports the matrix to a Matrix Market file of the given filename.
        """
        self.fillComplete()
        EpetraExt.RowMatrixToMatrixMarketFile(filename, self.matrix)

    @property
    def numpyArray(self):
        import tempfile
        import os
        from scipy.io import mmio
        from fipy.tools import parallelComm

        if parallelComm.procID == 0:
            (f, mtxName) = tempfile.mkstemp(suffix='.mtx')
        else:
            mtxName = None

        mtxName = parallelComm.bcast(mtxName)

        self.exportMmf(mtxName)

        parallelComm.Barrier()
        mtx = mmio.mmread(mtxName)
        parallelComm.Barrier()

        if parallelComm.procID == 0:
            os.remove(mtxName)

        coo = mtx.tocoo()
        trilinosMatrix = self.matrix
        numpyArray = numerix.zeros((trilinosMatrix.NumGlobalRows(), trilinosMatrix.NumGlobalCols()), 'd')
        numpyArray[coo.row, coo.col] = coo.data
        return numpyArray

    def _getDistributedMatrix(self):
        """
        Returns an equivalent Trilinos matrix, but redistributed evenly over
        all processors.
        """
        if self.comm.NumProc() == 1:
            return self.matrix
            # No redistribution necessary in serial mode
        else:
##            self._matrix.GlobalAssemble()
            totalElements = self.matrix.NumGlobalRows()

            DistributedMap = Epetra.Map(totalElements, 0, self.comm)
            RootToDist = Epetra.Import(DistributedMap, self.rangeMap)

            DistMatrix = Epetra.CrsMatrix(Epetra.Copy, DistributedMap, self.bandwidth*3/2)

            DistMatrix.Import(self.matrix, RootToDist, Epetra.Insert)

            return DistMatrix

    def fillComplete(self):
        if not self.matrix.Filled():
            self.matrix.FillComplete(self.domainMap, self.rangeMap)

    def finalize(self):
        self.fillComplete()
        self.matrix.OptimizeStorage()

class _TrilinosMatrixFromShape(_TrilinosMatrix):
    def __init__(self, rows, cols, bandwidth=1, sizeHint=None,
                 rowMap=None, colMap=None, domainMap=None):
        """Instantiants and wraps an Epetra.CrsMatrix

        :Parameters:
          - `rows`: The number of matrix rows
          - `cols`: The number of matrix columns
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `map`: The Epetra `Map` for the rows that this processor holds
        """
        size = max(rows, cols)
        if sizeHint is not None and bandwidth == 0:
            bandwidth = (sizeHint + size - 1) / (size or 1)
        else:
            bandwidth = bandwidth

        if rowMap is None:
            comm = Epetra.SerialComm()
            # Matrix building gets done on all processors
            rowMap = Epetra.Map(rows, rows, 0, comm)
        else:
            comm = rowMap.Comm()

        if colMap is None:
            colMap = Epetra.Map(cols, cols, 0, comm)

        matrix = Epetra.CrsMatrix(Epetra.Copy, rowMap, (bandwidth*3)//2)

        # Leave extra bandwidth, to handle multiple insertions into the
        # same spot. It's memory-inefficient, but it'll get cleaned up when
        # FillComplete is called, and according to the Trilinos devs the
        # performance boost will be worth it.

        _TrilinosMatrix.__init__(self,
                                     matrix=matrix,
                                     rowMap=rowMap,
                                     colMap=colMap,
                                     domainMap=domainMap,
                                     bandwidth=bandwidth)

class _TrilinosMeshMatrix(_TrilinosMatrixFromShape):
    def __init__(self, mesh, bandwidth=0, sizeHint=None, numberOfVariables=1, numberOfEquations=1):
        """Creates a `_TrilinosMatrixFromShape` associated with a `Mesh`

        :Parameters:
          - `mesh`: The `Mesh` to assemble the matrix for.
          - `bandwidth`: The proposed band width of the matrix.
          - `sizeHint`: estimate of the number of non-zeros
          - `numberOfVariables`: The columns of the matrix is determined by numberOfVariables * self.mesh.globalNumberOfCells.
          - `numberOfEquations`: The rows of the matrix is determined by numberOfEquations * self.mesh.globalNumberOfCells.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations

        comm = mesh.communicator.epetra_comm
        rowMap = Epetra.Map(-1, list(self._globalNonOverlappingRowIDs), 0, comm)
        colMap = Epetra.Map(-1, list(self._globalOverlappingColIDs), 0, comm)
        domainMap = rowMap

        _TrilinosMatrixFromShape.__init__(self,
                                 rows=self.numberOfEquations * self.mesh.globalNumberOfCells,
                                 cols=self.numberOfVariables * self.mesh.globalNumberOfCells,
                                 bandwidth=bandwidth,
                                 sizeHint=sizeHint,
                                 rowMap=rowMap,
                                 colMap=colMap,
                                 domainMap=domainMap)

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
    def _localNonOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localNonOverlappingCellIDs)

    def copy(self):
        tmp = _TrilinosMatrixFromShape.copy(self)
        copy = self.__class__(mesh=self.mesh, bandwidth=self.bandwidth)
        copy.matrix = tmp._matrix
        return copy

    def asTrilinosMeshMatrix(self):
        self.finalize()
        return self

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

    def flush(self):
        pass

    def _getMatrixProperty(self):
        if not hasattr(self, '_matrix'):
            self._matrix = _TrilinosMeshMatrix(self.mesh,
                                               bandwidth=self.bandwidth,
                                               numberOfVariables=self.numberOfVariables,
                                               numberOfEquations=self.numberOfEquations).matrix
        return super(_TrilinosMeshMatrix, self).matrix

    matrix = property(_getMatrixProperty, _TrilinosMatrixFromShape._setMatrix)

    def put(self, vector, id1, id2):
        vector, id1, id2 = self._globalNonOverlapping(vector, id1, id2)
        _TrilinosMatrixFromShape.put(self, vector=vector, id1=id1, id2=id2)

    def addAt(self, vector, id1, id2):
        vector, id1, id2 = self._globalNonOverlapping(vector, id1, id2)
        _TrilinosMatrixFromShape.addAt(self, vector=vector, id1=id1, id2=id2)

    def takeDiagonal(self):
        nonoverlapping_result = _TrilinosMatrixFromShape.takeDiagonal(self)

        overlapping_result = Epetra.Vector(self.colMap)
        overlapping_result.Import(nonoverlapping_result,
                                  Epetra.Import(self.colMap,
                                                self.domainMap),
                                  Epetra.Insert)

        return overlapping_result

    def __mul__(self, other):
        """
        Multiply a sparse matrix by another sparse matrix.

            >>> L1 = _TrilinosMatrixFromShape(rows=3, cols=3)
            >>> L1.addAt((3,10,numerix.pi,2.5), (0,0,1,2), (2,1,1,0))
            >>> L2 = _TrilinosIdentityMatrix(size=3)
            >>> L2.addAt((4.38,12357.2,1.1), (2,1,0), (1,0,2))

            >>> tmp = numerix.array(((1.23572000e+05, 2.31400000e+01, 3.00000000e+00),
            ...                      (3.88212887e+04, 3.14159265e+00, 0.00000000e+00),
            ...                      (2.50000000e+00, 0.00000000e+00, 2.75000000e+00)))

            >>> L = (L1 * L2).numpyArray

            >>> print numerix.allclose(tmp, L)
            True

        or a sparse matrix by a vector

            >>> tmp = numerix.array((29., 6.28318531, 2.5))
            >>> print numerix.allclose(L1 * numerix.array((1,2,3),'d'), tmp) # doctest: +SERIAL
            True

        or a vector by a sparse matrix

            >>> tmp = numerix.array((7.5, 16.28318531,  3.))
            >>> numerix.allclose(numerix.array((1,2,3),'d') * L1, tmp) # doctest: +SERIAL
            True

        Should be able to multiply an overlapping value obtained from a
        CellVariable. This is required to make the '--no-pysparse' flag
        work correctly.

            >>> from fipy import *
            >>> m = Grid1D(nx=6)
            >>> v0 = CellVariable(mesh=m, value=numerix.arange(m.globalNumberOfCells, dtype=float))
            >>> v1 = CellVariable(mesh=m, value=_TrilinosIdentityMeshMatrix(mesh=m) * v0.value)
            >>> print numerix.allclose(v0, v1)
            True

        """
        self.fillComplete()

        N = self.matrix.NumMyCols()

        if isinstance(other, _TrilinosMatrix):
            return _TrilinosMatrixFromShape.__mul__(self, other=other)
        else:
            shape = numerix.shape(other)


            if shape == ():
                result = self.copy()
                result.matrix.Scale(other)
                return result
            else:

                if isinstance(other, Epetra.Vector):
                    other_map = other.Map()
                else:
                    other_map = self.colMap

                if other_map.SameAs(self.colMap):
                    localNonOverlappingColIDs = self._localNonOverlappingColIDs

                    other = Epetra.Vector(self.domainMap,
                                          other[localNonOverlappingColIDs])

                if other.Map().SameAs(self.matrix.DomainMap()):
                    nonoverlapping_result = Epetra.Vector(self.rangeMap)
                    self.matrix.Multiply(False, other, nonoverlapping_result)

                    if other_map.SameAs(self.colMap):
                        overlapping_result = Epetra.Vector(self.colMap)
                        overlapping_result.Import(nonoverlapping_result,
                                                  Epetra.Import(self.colMap,
                                                                self.domainMap),
                                                  Epetra.Insert)

                        return overlapping_result
                    else:
                        return nonoverlapping_result

                else:
                    raise TypeError("%s: %s != (%d,)" % (self.__class__, str(shape), N))

    def _test(self):
        """Tests

        >>> from fipy import *
        >>> matrix = _TrilinosMeshMatrix(mesh=Grid1D(nx=5), numberOfVariables=3, numberOfEquations=2)
        >>> GOC = matrix._globalOverlappingColIDs
        >>> GNOC = matrix._globalNonOverlappingColIDs
        >>> LNOC = matrix._localNonOverlappingColIDs
        >>> GOR = matrix._globalOverlappingRowIDs
        >>> GNOR = matrix._globalNonOverlappingRowIDs
        >>> LNOR = matrix._localNonOverlappingRowIDs

        5 cells, 3 variables, 1 processor

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _globalOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _globalNonOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _localOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _localNonOverlappingColIDs:0

        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])  # doctest: +SERIAL
        True
        >>> print numerix.allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]) # doctest: +SERIAL
        True


        5 cells, 2 equations, 1 processor

        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _globalOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _globalNonOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _localOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _localNonOverlappingRowIDs:0

        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) # doctest: +SERIAL
        True


        5 cells, 3 variables, 2 processors

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        0  1  2  3     0  1  2  3     0  1  2  3      _globalOverlappingCellIDs:0
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:1

        0  1  2  3     5  6  7  8    10 11 12 13      _globalOverlappingColIDs:0
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _globalOverlappingColIDs:1

        0  1           0  1           0  1            _globalNonOverlappingCellIDs:0
              2  3  4        2  3  4        2  3  4   _globalNonOverlappingCellIDs:1

        0  1           5  6          10 11            _globalNonOverlappingColIDs:0
              2  3  4        7  8  9       12 13 14   _globalNonOverlappingColIDs:1

        0  1  2  3     0  1  2  3     0  1  2  3      _localOverlappingCellIDs:0
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:1

        0  1  2  3     4  5  6  7     8  9 10 11      _localOverlappingColIDs:0
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _localOverlappingColIDs:1

        0  1           0  1           0  1            _localNonOverlappingCellIDs:0
              2  3  4        2  3  4        2  3  4   _localNonOverlappingCellIDs:1

        0  1           4  5           8  9            _localNonOverlappingColIDs:0
              2  3  4        7  8  9       12 13 14   _localNonOverlappingColIDs:1


        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(GNOC, [0, 1, 5, 6, 10, 11]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GNOC, [2, 3, 4, 7, 8, 9, 12, 13, 14]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(LNOC, [0, 1, 4, 5, 8, 9]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(LNOC, [2, 3, 4, 7, 8, 9, 12, 13, 14]) # doctest: +PROCESSOR_1_OF_2
        True


        5 cells, 2 equations, 2 processors

        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        0  1  2  3     0  1  2  3      _globalOverlappingCellIDs:0
        0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:1

        0  1  2  3     5  6  7  8      _globalOverlappingRowIDs:0
        0  1  2  3  4  5  6  7  8  9   _globalOverlappingRowIDs:1

        0  1           0  1            _globalNonOverlappingCellIDs:0
              2  3  4        2  3  4   _globalNonOverlappingCellIDs:1

        0  1           5  6            _globalNonOverlappingRowIDs:0
              2  3  4        7  8  9   _globalNonOverlappingRowIDs:1

        0  1  2  3     0  1  2  3      _localOverlappingCellIDs:0
        0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:1

        0  1  2  3     4  5  6  7      _localOverlappingRowIDs:0
        0  1  2  3  4  5  6  7  8  9   _localOverlappingRowIDs:1

        0  1           0  1            _localNonOverlappingCellIDs:0
              2  3  4        2  3  4   _localNonOverlappingCellIDs:1

        0  1           4  5            _localNonOverlappingRowIDs:0
              2  3  4        7  8  9   _localNonOverlappingRowIDs:1


        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 5, 6, 7, 8]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(GNOR, [0, 1, 5, 6]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GNOR, [2, 3, 4, 7, 8, 9]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(LNOR, [0, 1, 4, 5]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(LNOR, [2, 3, 4, 7, 8, 9]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> matrix = _TrilinosMeshMatrix(mesh=Grid1D(nx=5, communicator=serialComm), numberOfVariables=3, numberOfEquations=2)
        >>> GOC = matrix._globalOverlappingColIDs
        >>> GNOC = matrix._globalNonOverlappingColIDs
        >>> LNOC = matrix._localNonOverlappingColIDs
        >>> GOR = matrix._globalOverlappingRowIDs
        >>> GNOR = matrix._globalNonOverlappingRowIDs
        >>> LNOR = matrix._localNonOverlappingRowIDs

        5 cells, 3 variables, serial

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _globalOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _globalNonOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _localOverlappingColIDs:0

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   _localNonOverlappingColIDs:0

        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        True
        >>> print numerix.allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        True
        >>> print numerix.allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14])
        True


        5 cells, 2 equations, serial

        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        0  1  2  3  4  0  1  2  3  4   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _globalOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _globalNonOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _localOverlappingRowIDs:0

        0  1  2  3  4  0  1  2  3  4   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9   _localNonOverlappingRowIDs:0

        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        True
        >>> print numerix.allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        True
        >>> print numerix.allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        True

        >>> matrix = _TrilinosMeshMatrix(mesh=Grid1D(nx=7), numberOfVariables=3, numberOfEquations=2)
        >>> GOC = matrix._globalOverlappingColIDs
        >>> GNOC = matrix._globalNonOverlappingColIDs
        >>> LNOC = matrix._localNonOverlappingColIDs
        >>> GOR = matrix._globalOverlappingRowIDs
        >>> GNOR = matrix._globalNonOverlappingRowIDs
        >>> LNOR = matrix._localNonOverlappingRowIDs

        7 cells, 3 variables, 1 processor

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   _globalOverlappingColIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   _globalNonOverlappingColIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   _localOverlappingColIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   _localNonOverlappingColIDs:0

        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        ...                              11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        ...                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
        ...                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20]) # doctest: +SERIAL
        True


        7 cells, 2 equations, 1 processor

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   _globalOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   _globalOverlappingRowIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   _globalNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   _globalNonOverlappingRowIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   _localOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   _localOverlappingRowIDs:0

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   _localNonOverlappingCellIDs:0

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   _localNonOverlappingRowIDs:0

        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) # doctest: +SERIAL
        True
        >>> print numerix.allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]) # doctest: +SERIAL
        True


        7 cells, 3 variables, 2 processors

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        0  1  2  3  4        0  1  2  3  4        0  1  2  3  4         _globalOverlappingCellIDs:0
           1  2  3  4  5  6     1  2  3  4  5  6     1  2  3  4  5  6   _globalOverlappingCellIDs:1

        0  1  2  3  4        7  8  9 10 11       14 15 16 17 18         _globalOverlappingColIDs:0
           1  2  3  4  5  6     8  9 10 11 12 13    15 16 17 18 19 20   _globalOverlappingColIDs:1

        0  1  2              0  1  2              0  1  2               _globalNonOverlappingCellIDs:0
                 3  4  5  6           3  4  5  6           3  4  5  6   _globalNonOverlappingCellIDs:1

        0  1  2              7  8  9             14 15 16               _globalNonOverlappingColIDs:0
                 3  4  5  6          10 11 12 13          17 18 19 20   _globalNonOverlappingColIDs:1

        0  1  2  3  4        0  1  2  3  4        0  1  2  3  4         _localOverlappingCellIDs:0
           0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5   _localOverlappingCellIDs:1

        0  1  2  3  4        5  6  7  8  9       10 11 12 13 14         _localOverlappingColIDs:0
           0  1  2  3  4  5     6  7  8  9 10 11    12 13 14 15 16 17   _localOverlappingColIDs:1

        0  1  2              0  1  2              0  1  2               _localNonOverlappingCellIDs:0
                 2  3  4  5           2  3  4  5           2  3  4  5   _localNonOverlappingCellIDs:1

        0  1  2              5  6  7             10 11 12               _localNonOverlappingColIDs:0
                 2  3  4  5           8  9 10 11          14 15 16 17   _localNonOverlappingColIDs:1

        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 14, 15, 16, 17, 18]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GOC, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(GNOC, [0, 1, 2, 7, 8, 9, 14, 15, 16]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GNOC, [3, 4, 5, 6, 10, 11, 12, 13, 17, 18, 19, 20]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(LNOC, [0, 1, 2, 5, 6, 7, 10, 11, 12]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(LNOC, [2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 16, 17]) # doctest: +PROCESSOR_1_OF_2
        True

        7 cells, 2 equations, 2 processors

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        0  1  2  3  4        0  1  2  3  4         _globalOverlappingCellIDs:0
           1  2  3  4  5  6     1  2  3  4  5  6   _globalOverlappingCellIDs:1

        0  1  2  3  4        7  8  9 10 11         _globalOverlappingRowIDs:0
           1  2  3  4  5  6     8  9 10 11 12 13   _globalOverlappingRowIDs:1

        0  1  2              0  1  2               _globalNonOverlappingCellIDs:0
                 3  4  5  6           3  4  5  6   _globalNonOverlappingCellIDs:1

        0  1  2              7  8  9               _globalNonOverlappingRowIDs:0
                 3  4  5  6          10 11 12 13   _globalNonOverlappingRowIDs:1

        0  1  2  3  4        0  1  2  3  4         _localOverlappingCellIDs:0
           0  1  2  3  4  5     0  1  2  3  4  5   _localOverlappingCellIDs:1

        0  1  2  3  4        5  6  7  8  9         _localOverlappingRowIDs:0
           0  1  2  3  4  5     6  7  8  9 10 11   _localOverlappingRowIDs:1

        0  1  2              0  1  2               _localNonOverlappingCellIDs:0
                 2  3  4  5           2  3  4  5   _localNonOverlappingCellIDs:1

        0  1  2              5  6  7               _localNonOverlappingRowIDs:0
                 2  3  4  5           8  9 10 11   _localNonOverlappingRowIDs:1

        >>> print numerix.allequal(GOR, [0, 1, 2, 3, 4, 7, 8, 9, 10, 11]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GOR, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(GNOR, [0, 1, 2, 7, 8, 9]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(GNOR, [3, 4, 5, 6, 10, 11, 12, 13]) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print numerix.allequal(LNOR, [0, 1, 2, 5, 6, 7]) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print numerix.allequal(LNOR, [2, 3, 4, 5, 8, 9, 10, 11]) # doctest: +PROCESSOR_1_OF_2
        True


        7 cells, 3 variables, 3 processors

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        0  1  2  3           0  1  2  3           0  1  2  3            _globalOverlappingCellIDs:0
        0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5      _globalOverlappingCellIDs:1
              2  3  4  5  6        2  3  4  5  6        2  3  4  5  6   _globalOverlappingCellIDs:2

        0  1  2  3           7  8  9 10          14 15 16 17            _globalOverlappingColIDs:0
        0  1  2  3  4  5     7  8  9 10 11 12    14 15 16 17 18 19      _globalOverlappingColIDs:1
              2  3  4  5  6        9 10 11 12 13       16 17 18 19 20   _globalOverlappingColIDs:2

        0  1                 0  1                 0  1                  _globalNonOverlappingCellIDs:0
              2  3                 2  3                 2  3            _globalNonOverlappingCellIDs:1
                    4  5  6              4  5  6              4  5  6   _globalNonOverlappingCellIDs:2

        0  1                 7  8                14 15                  _globalNonOverlappingColIDs:0
              2  3                 9 10                16 17            _globalNonOverlappingColIDs:1
                    4  5  6             11 12 13             18 19 20   _globalNonOverlappingColIDs:2

        0  1  2  3           0  1  2  3           0  1  2  3            _localOverlappingCellIDs:0
        0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5      _localOverlappingCellIDs:1
              0  1  2  3  4        0  1  2  3  4        0  1  2  3  4   _localOverlappingCellIDs:2

        0  1  2  3           4  5  6  7           8  9 10 11            _localOverlappingColIDs:0
        0  1  2  3  4  5     6  7  8  9 10 11    12 13 14 15 16 17      _localOverlappingColIDs:1
              0  1  2  3  4        5  6  7  8  9       10 11 12 13 14   _localOverlappingColIDs:2

        0  1                 0  1                 0  1                  _localNonOverlappingCellIDs:0
              2  3                 2  3                 2  3            _localNonOverlappingCellIDs:1
                    2  3  4              2  3  4              2  3  4   _localNonOverlappingCellIDs:2

        0  1                 4  5                 8  9                  _localNonOverlappingColIDs:0
              2  3                 8  9                14 15            _localNonOverlappingColIDs:1
                    2  3  4              7  8  9             12 13 14   _localNonOverlappingColIDs:2

        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 7, 8, 9, 10, 14, 15, 16, 17]) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print numerix.allequal(GOC, [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print numerix.allequal(GOC, [2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20]) # doctest: +PROCESSOR_2_OF_3
        True

        >>> print numerix.allequal(GNOC, [0, 1, 7, 8, 14, 15]) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print numerix.allequal(GNOC, [2, 3, 9, 10, 16, 17]) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print numerix.allequal(GNOC, [4, 5, 6, 11, 12, 13, 18, 19, 20]) # doctest: +PROCESSOR_2_OF_3
        True

        >>> print numerix.allequal(LNOC, [0, 1, 4, 5, 8, 9]) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print numerix.allequal(LNOC, [2, 3, 8, 9, 14, 15]) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print numerix.allequal(LNOC, [2, 3, 4, 7, 8, 9, 12, 13, 14]) # doctest: +PROCESSOR_2_OF_3
        True


        7 cells, 2 equations, 3 processors

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        0  1  2  3           0  1  2  3            _globalOverlappingCellIDs:0
        0  1  2  3  4  5     0  1  2  3  4  5      _globalOverlappingCellIDs:1
              2  3  4  5  6        2  3  4  5  6   _globalOverlappingCellIDs:2

        0  1  2  3           7  8  9 10            _globalOverlappingRowIDs:0
        0  1  2  3  4  5     7  8  9 10 11 12      _globalOverlappingRowIDs:1
              2  3  4  5  6        9 10 11 12 13   _globalOverlappingRowIDs:2

        0  1                 0  1                  _globalNonOverlappingCellIDs:0
              2  3                 2  3            _globalNonOverlappingCellIDs:1
                    4  5  6              4  5  6   _globalNonOverlappingCellIDs:2

        0  1                 7  8                  _globalNonOverlappingRowIDs:0
              2  3                 9 10            _globalNonOverlappingRowIDs:1
                    4  5  6             11 12 13   _globalNonOverlappingRowIDs:2

        0  1  2  3           0  1  2  3            _localOverlappingCellIDs:0
        0  1  2  3  4  5     0  1  2  3  4  5      _localOverlappingCellIDs:1
              0  1  2  3  4        0  1  2  3  4   _localOverlappingCellIDs:2

        0  1  2  3           4  5  6  7            _localOverlappingRowIDs:0
        0  1  2  3  4  5     6  7  8  9 10 11      _localOverlappingRowIDs:1
              0  1  2  3  4        5  6  7  8  9   _localOverlappingRowIDs:2

        0  1                 0  1                  _localNonOverlappingCellIDs:0
              2  3                 2  3            _localNonOverlappingCellIDs:1
                    2  3  4              2  3  4   _localNonOverlappingCellIDs:2

        0  1                 4  5                  _localNonOverlappingRowIDs:0
              2  3                 8  9            _localNonOverlappingRowIDs:1
                    2  3  4              7  8  9   _localNonOverlappingRowIDs:2

        """
        pass


class _TrilinosIdentityMatrix(_TrilinosMatrixFromShape):
    """
    Represents a sparse identity matrix for Trilinos.
    """
    def __init__(self, size):
        """
        Create a sparse matrix with '1' in the diagonal

            >>> print _TrilinosIdentityMatrix(size=3)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _TrilinosMatrixFromShape.__init__(self, rows=size, cols=size, bandwidth=1)
        ids = numerix.arange(size)
        self.addAt(numerix.ones(size, 'l'), ids, ids)

class _TrilinosIdentityMeshMatrix(_TrilinosMeshMatrix):
    def __init__(self, mesh):
        """
        Create a sparse matrix associated with a `Mesh` with '1' in the diagonal

            >>> from fipy import Grid1D
            >>> mesh = Grid1D(nx=3)
            >>> print _TrilinosIdentityMeshMatrix(mesh=mesh)
             1.000000      ---        ---    
                ---     1.000000      ---    
                ---        ---     1.000000  
        """
        _TrilinosMeshMatrix.__init__(self, mesh=mesh, bandwidth=1)
        size = mesh.numberOfCells
        ids = numerix.arange(size)
        self.addAt(numerix.ones(size, 'l'), ids, ids)

class _TrilinosMeshMatrixKeepStencil(_TrilinosMeshMatrix):

    def _getStencil(self, id1, id2):
        if not hasattr(self, 'stencil'):
            self.stencil = _TrilinosMeshMatrix._getStencil(self, id1, id2)

        return self.stencil

    def flush(self, cacheStencil=False):
        """Deletes the matrix but maintains the stencil used
        `_globalNonOverlapping()` in as it can be expensive to construct.

        :Parameters:
          - `cacheStencil`: Boolean value to determine whether to keep the stencil (tuple of IDs and a mask) even after deleting the matrix.

        """

        del self._matrix
        if not cacheStencil:
            del self.stencil

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
