from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

class _SparseMatrix(object):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    matrix = None
    numpyArray = property()
    _shape = property()

    def __init__(self):
        pass

    __array_priority__ = 100.0

    def __array_wrap(self, arr, context=None):
        if context is None:
            return arr
        else:
            return NotImplemented

    def copy(self):
        raise NotImplementedError

    def __getitem__(self, index):
        raise NotImplementedError

    @property
    def _range(self):
        raise NotImplementedError

    def __str__(self):
        s = ''
        cellWidth = 11
        Irange, Jrange = self._range

        for j in Jrange:
            for i in Irange:
                v = self[j, i]
                if v == 0:
                    s += "---".center(cellWidth)
                else:
                    exp = numerix.log(abs(v))
                    if abs(exp) <= 4:
                        if exp < 0:
                            s += ("%9.6f" % v).ljust(cellWidth)
                        else:
                            s += ("%9.*f" % (6, v)).ljust(cellWidth)
                    else:
                        s += ("%9.2e" % v).ljust(cellWidth)
            s += "\n"
        return s[:-1]

    def __repr__(self):
        return repr(self.matrix)

    def __setitem__(self, index, value):
        raise NotImplementedError

    def __add__(self, other):
        raise NotImplementedError

    __radd__ = __add__

    def __iadd__(self, other):
        raise NotImplementedError

    def __sub__(self, other):
        raise NotImplementedError

    # Ask about this rsub
    def __rsub__(self, other):
        return -(self.__sub__(other))


    def __isub__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        raise NotImplementedError

    def __neg__(self):
        return self * -1

    def __pos__(self):
        return self

##     def __eq__(self,other):
## 	return self.matrix.__eq__(other.matrix)

    def put(self, vector, id1, id2):
        raise NotImplementedError

    def putDiagonal(self, vector):
        raise NotImplementedError

    def take(self, id1, id2):
        raise NotImplementedError

    def takeDiagonal(self):
        raise NotImplementedError

    def addAt(self, vector, id1, id2):
        raise NotImplementedError

    def addAtDiagonal(self, vector):
        raise NotImplementedError

    def exportMmf(self, filename):
        raise NotImplementedError

    @property
    def CSR(self):
        raise NotImplementedError

    @property
    def LIL(self):
        raise NotImplementedError

    @property
    def T(self):
        raise NotImplementedError

    def _matrix2mesh(self, ids):
        """Convert matrix row indices to mesh cell indices
        """
        return ids

    def _mesh2matrix(self, ids):
        """Convert mesh cell indices to matrix row indices
        """
        return ids

##     def __array__(self):
##      shape = self._shape
##      indices = numerix.indices(shape)
##         numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
##      return numerix.reshape(numMatrix, shape)

class _Mesh2Matrix(object):
    _bodies = None
    _ghosts = None

    def __init__(self, mesh, matrix,
                 numberOfEquations=1, numberOfVariables=1):
        """Creates a mapping between mesh cells and matrix rows and columns

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to be mapped to a `_SparseMatrix`
        matrix : ~fipy.matrices.sparseMatrix.SparseMatrix
            The `_SparseMatrix` to be mapped to a `Mesh`
        numberOfVariables : int
            The local columns of the matrix are determined by
            `numberOfVariables * len(mesh._localNonOverlappingCellIDs)`.
        numberOfEquations : int
            The local rows of the matrix are determined by
            `numberOfEquations * len(mesh._localNonOverlappingCellIDs)`.
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations

        # Note: storing matrix directly results in a reference cycle
        # between _PETScBaseMeshMatrix and _Mesh2Matrix, specifically
        # between _PETScBaseMeshMatrix.matrix and _PETScBaseMeshMatrix._ao.
        #
        # See https://lists.mcs.anl.gov/pipermail/petsc-users/2020-October/042652.html
        # and https://github.com/usnistgov/fipy/pull/761
        import weakref
        self.matrix = weakref.ref(matrix)

    @staticmethod
    def _cellIDsToGlobalIDs(IDs, M, L):
        N = len(IDs)
        return (numerix.vstack([IDs] * M) + numerix.indices((M, N))[0] * L).flatten()

    def _cellIDsToGlobalRowIDs(self, IDs):
        return self._cellIDsToGlobalIDs(IDs, M=self.numberOfEquations,
                                        L=self.mesh.globalNumberOfCells)

    def _cellIDsToLocalRowIDs(self, IDs):
        return self._cellIDsToGlobalIDs(IDs, M=self.numberOfEquations,
                                        L=self.mesh.numberOfCells)

    @property
    def globalNonOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def globalOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def localNonOverlappingRowIDs(self):
        return self._cellIDsToLocalRowIDs(self.mesh._localNonOverlappingCellIDs)

    def _cellIDsToGlobalColIDs(self, IDs):
        return self._cellIDsToGlobalIDs(IDs, M=self.numberOfVariables,
                                        L=self.mesh.globalNumberOfCells)

    def _cellIDsToLocalColIDs(self, IDs):
        return self._cellIDsToGlobalIDs(IDs, M=self.numberOfVariables,
                                        L=self.mesh.numberOfCells)

    @property
    def globalNonOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def globalOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def localOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localOverlappingCellIDs)

    @property
    def localNonOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localNonOverlappingCellIDs)

    def _getStencil_(self, id1, id2,
                     globalOverlappihgIDs, globalNonOverlappihgIDs,
                     overlapping=False):
        id1 = globalOverlappihgIDs[id1]

        if overlapping:
            mask = numerix.ones(id1.shape, dtype=bool)
        else:
            mask = numerix.isin(id1, globalNonOverlappihgIDs)

        id1 = self.matrix()._mesh2matrix(id1[mask])
        id2 = numerix.asarray(id2)[mask]

        return id1, id2, mask

    def _getStencil(self, id1, id2, overlapping=False):
        raise NotImplementedError

    def globalVectorAndIDs(self, vector, id1, id2, overlapping=False):
        """Transforms local overlapping values and coordinates to global

        Parameters
        ----------
        vector : array_like
            The overlapping values to insert.
        id1 : array_like
            The local overlapping row indices.
        id2 : array_like
            The local overlapping column indices.
        overlapping : bool
            Whether to return ghosted values or not (default False)

        Returns
        -------
        tuple of ndarray
            (vector, global row indices, global column indices)
        """
        id1, id2, mask = self._getStencil(id1, id2, overlapping)
        vector = vector[mask]
        return (vector, id1, id2)

    @property
    def bodies(self):
        if self._bodies is None:
            self._bodies = numerix.isin(self.mesh._globalOverlappingCellIDs,
                                        self.mesh._globalNonOverlappingCellIDs)
        return self._bodies

    @property
    def ghosts(self):
        if self._ghosts is None:
            self._ghosts = self.mesh._globalOverlappingCellIDs[~self.bodies]
            self._ghosts = self._cellIDsToGlobalRowIDs(self._ghosts)
            self._ghosts = self.matrix()._mesh2matrix(self._ghosts)

        return self._ghosts

class _RowMesh2Matrix(_Mesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        return self._getStencil_(id1, id2,
                                 self.globalOverlappingRowIDs,
                                 self.globalNonOverlappingRowIDs,
                                 overlapping)

class _ColMesh2Matrix(_Mesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        id2, id1, mask = self._getStencil_(id2, id1,
                                           self.globalOverlappingColIDs,
                                           self.globalNonOverlappingColIDs,
                                           overlapping)

        return id1, id2, mask

class _RowColMesh2Matrix(_RowMesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        id2 = self.globalOverlappingColIDs[id2]

        id1, id2, mask = super(_RowColMesh2Matrix, self)._getStencil(id1, id2, overlapping)

        id2 = self.matrix()._mesh2matrix(id2)

        return id1, id2, mask

    def _test(self):
        """Tests

        >>> from fipy import Grid1D
        >>> from fipy.tools.numerix import allequal

        >>> m2m = _RowColMesh2Matrix(mesh=Grid1D(nx=5),
        ...                          numberOfVariables=3,
        ...                          numberOfEquations=2)
        >>> GOC = m2m.globalOverlappingColIDs
        >>> GNOC = m2m.globalNonOverlappingColIDs
        >>> LNOC = m2m.localNonOverlappingColIDs
        >>> GOR = m2m.globalOverlappingRowIDs
        >>> GNOR = m2m.globalNonOverlappingRowIDs
        >>> LNOR = m2m.localNonOverlappingRowIDs

        5 cells, 3 variables, 1 processor

        ```
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        globalOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        globalNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        localOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        localNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0
        ```

        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                      8, 9, 10, 11, 12, 13, 14]))  # doctest: +SERIAL
        True
        >>> print(allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                       8, 9, 10, 11, 12, 13, 14])) # doctest: +SERIAL
        True
        >>> print(allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                       8, 9, 10, 11, 12, 13, 14])) # doctest: +SERIAL
        True


        5 cells, 2 equations, 1 processor

        ```
        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        globalOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        globalNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        _localOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        localNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0
        ```

        >>> print(allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])) # doctest: +SERIAL
        True
        >>> print(allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])) # doctest: +SERIAL
        True
        >>> print(allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])) # doctest: +SERIAL
        True


        5 cells, 3 variables, 2 processors

        ```
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        _globalOverlappingCellIDs

        0  1  2  3     0  1  2  3     0  1  2  3      proc 0
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc  1

        globalOverlappingColIDs

        0  1  2  3     5  6  7  8    10 11 12 13      proc 0
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc  1

        _globalNonOverlappingCellIDs

        0  1           0  1           0  1            proc 0
              2  3  4        2  3  4        2  3  4   proc  1

        globalNonOverlappingColIDs

        0  1           5  6          10 11            proc 0
              2  3  4        7  8  9       12 13 14   proc  1

        _localOverlappingCellIDs

        0  1  2  3     0  1  2  3     0  1  2  3      proc 0
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc  1

        localOverlappingColIDs

        0  1  2  3     4  5  6  7     8  9 10 11      proc 0
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc  1

        _localNonOverlappingCellIDs

        0  1           0  1           0  1            proc 0
              2  3  4        2  3  4        2  3  4   proc  1

        localNonOverlappingColIDs

        0  1           4  5           8  9            proc 0
              2  3  4        7  8  9       12 13 14   proc  1
        ```


        >>> print(allequal(GOC, [0, 1, 2, 3, 5, 6, 7,
        ...                      8, 10, 11, 12, 13])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                      8, 9, 10, 11, 12, 13, 14])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(GNOC, [0, 1, 5, 6, 10, 11])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GNOC, [2, 3, 4, 7, 8, 9,
        ...                       12, 13, 14])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(LNOC, [0, 1, 4, 5, 8, 9])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(LNOC, [2, 3, 4, 7, 8, 9,
        ...                       12, 13, 14])) # doctest: +PROCESSOR_1_OF_2
        True


        5 cells, 2 equations, 2 processors

        ```
        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        _globalOverlappingCellIDs

        0  1  2  3     0  1  2  3      proc 0
        0  1  2  3  4  0  1  2  3  4   proc  1

        globalOverlappingRowIDs

        0  1  2  3     5  6  7  8      proc 0
        0  1  2  3  4  5  6  7  8  9   proc  1

        _globalNonOverlappingCellIDs

        0  1           0  1            proc 0
              2  3  4        2  3  4   proc  1

        globalNonOverlappingRowIDs

        0  1           5  6            proc 0
              2  3  4        7  8  9   proc  1

        _localOverlappingCellIDs

        0  1  2  3     0  1  2  3      proc 0
        0  1  2  3  4  0  1  2  3  4   proc  1

        _localOverlappingRowIDs

        0  1  2  3     4  5  6  7      proc 0
        0  1  2  3  4  5  6  7  8  9   proc  1

        _localNonOverlappingCellIDs

        0  1           0  1            proc 0
              2  3  4        2  3  4   proc  1

        localNonOverlappingRowIDs

        0  1           4  5            proc 0
              2  3  4        7  8  9   proc  1
        ```


        >>> print(allequal(GOR, [0, 1, 2, 3, 5, 6, 7, 8])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GOR, [0, 1, 2, 3, 4, 5,
        ...                      6, 7, 8, 9])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(GNOR, [0, 1, 5, 6])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GNOR, [2, 3, 4, 7, 8, 9])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(LNOR, [0, 1, 4, 5])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(LNOR, [2, 3, 4, 7, 8, 9])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> m2m = _RowColMesh2Matrix(mesh=Grid1D(nx=5, communicator=serialComm),
        ...                          numberOfVariables=3,
        ...                          numberOfEquations=2)
        >>> GOC = m2m.globalOverlappingColIDs
        >>> GNOC = m2m.globalNonOverlappingColIDs
        >>> LNOC = m2m.localNonOverlappingColIDs
        >>> GOR = m2m.globalOverlappingRowIDs
        >>> GNOR = m2m.globalNonOverlappingRowIDs
        >>> LNOR = m2m.localNonOverlappingRowIDs

        5 cells, 3 variables, serial

        ```
        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   column IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        globalOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        globalNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        localOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4  0  1  2  3  4   proc 0

        localNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14   proc 0
        ```

        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]))
        True
        >>> print(allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]))
        True
        >>> print(allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]))
        True


        5 cells, 2 equations, serial

        ```
        0  1  2  3  4  0  1  2  3  4   cell IDs
        0  1  2  3  4  5  6  7  8  9   row IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        globalOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        globalNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        _localOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  0  1  2  3  4   proc 0

        localNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9   proc 0
        ```

        >>> print(allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        True
        >>> print(allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        True
        >>> print(allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))
        True

        >>> m2m = _RowColMesh2Matrix(mesh=Grid1D(nx=7),
        ...                          numberOfVariables=3,
        ...                          numberOfEquations=2)
        >>> GOC = m2m.globalOverlappingColIDs
        >>> GNOC = m2m.globalNonOverlappingColIDs
        >>> LNOC = m2m.localNonOverlappingColIDs
        >>> GOR = m2m.globalOverlappingRowIDs
        >>> GNOR = m2m.globalNonOverlappingRowIDs
        >>> LNOR = m2m.localNonOverlappingRowIDs

        7 cells, 3 variables, 1 processor

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        globalOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        globalNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        localOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        localNonOverlappingColIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   proc 0
        ```

        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        ...                      12, 13, 14, 15, 16, 17, 18, 19, 20])) # doctest: +SERIAL
        True
        >>> print(allequal(GNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        ...                       12, 13, 14, 15, 16, 17, 18, 19, 20])) # doctest: +SERIAL
        True
        >>> print(allequal(LNOC, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        ...                       12, 13, 14, 15, 16, 17, 18, 19, 20])) # doctest: +SERIAL
        True


        7 cells, 2 equations, 1 processor

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        _globalOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        globalOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   proc 0

        _globalNonOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        globalNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   proc 0

        _localOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        _localOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   proc 0

        _localNonOverlappingCellIDs

        0  1  2  3  4  5  6  0  1  2  3  4  5  6   proc 0

        localNonOverlappingRowIDs

        0  1  2  3  4  5  6  7  8  9 10 11 12 13   proc 0
        ```

        >>> print(allequal(GOR, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                      8, 9, 10, 11, 12, 13])) # doctest: +SERIAL
        True
        >>> print(allequal(GNOR, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                       8, 9, 10, 11, 12, 13])) # doctest: +SERIAL
        True
        >>> print(allequal(LNOR, [0, 1, 2, 3, 4, 5, 6, 7,
        ...                       8, 9, 10, 11, 12, 13])) # doctest: +SERIAL
        True


        7 cells, 3 variables, 2 processors

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        _globalOverlappingCellIDs

        0  1  2  3  4        0  1  2  3  4        0  1  2  3  4         proc 0
           1  2  3  4  5  6     1  2  3  4  5  6     1  2  3  4  5  6   proc  1

        globalOverlappingColIDs

        0  1  2  3  4        7  8  9 10 11       14 15 16 17 18         proc 0
           1  2  3  4  5  6     8  9 10 11 12 13    15 16 17 18 19 20   proc  1

        _globalNonOverlappingCellIDs

        0  1  2              0  1  2              0  1  2               proc 0
                 3  4  5  6           3  4  5  6           3  4  5  6   proc  1

        globalNonOverlappingColIDs

        0  1  2              7  8  9             14 15 16               proc 0
                 3  4  5  6          10 11 12 13          17 18 19 20   proc  1

        _localOverlappingCellIDs

        0  1  2  3  4        0  1  2  3  4        0  1  2  3  4         proc 0
           0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5   proc  1

        localOverlappingColIDs

        0  1  2  3  4        5  6  7  8  9       10 11 12 13 14         proc 0
           0  1  2  3  4  5     6  7  8  9 10 11    12 13 14 15 16 17   proc  1

        _localNonOverlappingCellIDs

        0  1  2              0  1  2              0  1  2               proc 0
                 2  3  4  5           2  3  4  5           2  3  4  5   proc  1

        localNonOverlappingColIDs

        0  1  2              5  6  7             10 11 12               proc 0
                 2  3  4  5           8  9 10 11          14 15 16 17   proc  1
        ```

        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 7, 8, 9, 10,
        ...                      11, 14, 15, 16, 17, 18])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GOC, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12,
        ...                      13, 15, 16, 17, 18, 19, 20])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(GNOC, [0, 1, 2, 7, 8,
        ...                       9, 14, 15, 16])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GNOC, [3, 4, 5, 6, 10, 11, 12,
        ...                       13, 17, 18, 19, 20])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(LNOC, [0, 1, 2, 5, 6,
        ...                       7, 10, 11, 12])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(LNOC, [2, 3, 4, 5, 8, 9, 10,
        ...                       11, 14, 15, 16, 17])) # doctest: +PROCESSOR_1_OF_2
        True

        7 cells, 2 equations, 2 processors

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        _globalOverlappingCellIDs

        0  1  2  3  4        0  1  2  3  4         proc 0
           1  2  3  4  5  6     1  2  3  4  5  6   proc  1

        globalOverlappingRowIDs

        0  1  2  3  4        7  8  9 10 11         proc 0
           1  2  3  4  5  6     8  9 10 11 12 13   proc  1

        _globalNonOverlappingCellIDs

        0  1  2              0  1  2               proc 0
                 3  4  5  6           3  4  5  6   proc  1

        globalNonOverlappingRowIDs

        0  1  2              7  8  9               proc 0
                 3  4  5  6          10 11 12 13   proc  1

        _localOverlappingCellIDs

        0  1  2  3  4        0  1  2  3  4         proc 0
           0  1  2  3  4  5     0  1  2  3  4  5   proc  1

        _localOverlappingRowIDs

        0  1  2  3  4        5  6  7  8  9         proc 0
           0  1  2  3  4  5     6  7  8  9 10 11   proc  1

        _localNonOverlappingCellIDs

        0  1  2              0  1  2               proc 0
                 2  3  4  5           2  3  4  5   proc  1

        localNonOverlappingRowIDs

        0  1  2              5  6  7               proc 0
                 2  3  4  5           8  9 10 11   proc  1
        ```

        >>> print(allequal(GOR, [0, 1, 2, 3, 4, 7,
        ...                      8, 9, 10, 11])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GOR, [1, 2, 3, 4, 5, 6, 8, 9,
        ...                      10, 11, 12, 13])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(GNOR, [0, 1, 2, 7, 8, 9])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(GNOR, [3, 4, 5, 6, 10, 11, 12, 13])) # doctest: +PROCESSOR_1_OF_2
        True

        >>> print(allequal(LNOR, [0, 1, 2, 5, 6, 7])) # doctest: +PROCESSOR_0_OF_2
        True
        >>> print(allequal(LNOR, [2, 3, 4, 5, 8, 9, 10, 11])) # doctest: +PROCESSOR_1_OF_2
        True


        7 cells, 3 variables, 3 processors

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20   column IDs

        _globalOverlappingCellIDs

        0  1  2  3           0  1  2  3           0  1  2  3            proc 0
        0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5      proc  1
              2  3  4  5  6        2  3  4  5  6        2  3  4  5  6   proc   2

        globalOverlappingColIDs

        0  1  2  3           7  8  9 10          14 15 16 17            proc 0
        0  1  2  3  4  5     7  8  9 10 11 12    14 15 16 17 18 19      proc  1
              2  3  4  5  6        9 10 11 12 13       16 17 18 19 20   proc   2

        _globalNonOverlappingCellIDs

        0  1                 0  1                 0  1                  proc 0
              2  3                 2  3                 2  3            proc  1
                    4  5  6              4  5  6              4  5  6   proc   2

        globalNonOverlappingColIDs

        0  1                 7  8                14 15                  proc 0
              2  3                 9 10                16 17            proc  1
                    4  5  6             11 12 13             18 19 20   proc   2

        _localOverlappingCellIDs

        0  1  2  3           0  1  2  3           0  1  2  3            proc 0
        0  1  2  3  4  5     0  1  2  3  4  5     0  1  2  3  4  5      proc  1
              0  1  2  3  4        0  1  2  3  4        0  1  2  3  4   proc   2

        localOverlappingColIDs

        0  1  2  3           4  5  6  7           8  9 10 11            proc 0
        0  1  2  3  4  5     6  7  8  9 10 11    12 13 14 15 16 17      proc  1
              0  1  2  3  4        5  6  7  8  9       10 11 12 13 14   proc   2

        _localNonOverlappingCellIDs

        0  1                 0  1                 0  1                  proc 0
              2  3                 2  3                 2  3            proc  1
                    2  3  4              2  3  4              2  3  4   proc   2

        localNonOverlappingColIDs

        0  1                 4  5                 8  9                  proc 0
              2  3                 8  9                14 15            proc  1
                    2  3  4              7  8  9             12 13 14   proc   2
        ```

        >>> print(allequal(GOC, [0, 1, 2, 3, 7, 8, 9,
        ...                      10, 14, 15, 16, 17])) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print(allequal(GOC, [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11,
        ...                      12, 14, 15, 16, 17, 18, 19])) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print(allequal(GOC, [2, 3, 4, 5, 6, 9, 10, 11, 12,
        ...                      13, 16, 17, 18, 19, 20])) # doctest: +PROCESSOR_2_OF_3
        True

        >>> print(allequal(GNOC, [0, 1, 7, 8, 14, 15])) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print(allequal(GNOC, [2, 3, 9, 10, 16, 17])) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print(allequal(GNOC, [4, 5, 6, 11, 12, 13,
        ...                       18, 19, 20])) # doctest: +PROCESSOR_2_OF_3
        True

        >>> print(allequal(LNOC, [0, 1, 4, 5, 8, 9])) # doctest: +PROCESSOR_0_OF_3
        True
        >>> print(allequal(LNOC, [2, 3, 8, 9, 14, 15])) # doctest: +PROCESSOR_1_OF_3
        True
        >>> print(allequal(LNOC, [2, 3, 4, 7, 8, 9,
        ...                       12, 13, 14])) # doctest: +PROCESSOR_2_OF_3
        True


        7 cells, 2 equations, 3 processors

        ```
        0  1  2  3  4  5  6  0  1  2  3  4  5  6   cell IDs
        0  1  2  3  4  5  6  7  8  9 10 11 12 13   row IDs

        _globalOverlappingCellIDs

        0  1  2  3           0  1  2  3            proc 0
        0  1  2  3  4  5     0  1  2  3  4  5      proc  1
              2  3  4  5  6        2  3  4  5  6   proc   2

        globalOverlappingRowIDs

        0  1  2  3           7  8  9 10            proc 0
        0  1  2  3  4  5     7  8  9 10 11 12      proc  1
              2  3  4  5  6        9 10 11 12 13   proc   2

        _globalNonOverlappingCellIDs

        0  1                 0  1                  proc 0
              2  3                 2  3            proc  1
                    4  5  6              4  5  6   proc   2

        globalNonOverlappingRowIDs

        0  1                 7  8                  proc 0
              2  3                 9 10            proc  1
                    4  5  6             11 12 13   proc   2

        _localOverlappingCellIDs

        0  1  2  3           0  1  2  3            proc 0
        0  1  2  3  4  5     0  1  2  3  4  5      proc  1
              0  1  2  3  4        0  1  2  3  4   proc   2

        _localOverlappingRowIDs

        0  1  2  3           4  5  6  7            proc 0
        0  1  2  3  4  5     6  7  8  9 10 11      proc  1
              0  1  2  3  4        5  6  7  8  9   proc   2

        _localNonOverlappingCellIDs

        0  1                 0  1                  proc 0
              2  3                 2  3            proc  1
                    2  3  4              2  3  4   proc   2

        localNonOverlappingRowIDs

        0  1                 4  5                  proc 0
              2  3                 8  9            proc  1
                    2  3  4              7  8  9   proc   2
        ```
        """
        pass
