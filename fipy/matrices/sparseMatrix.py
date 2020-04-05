from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

class _SparseMatrix(object):

    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, mesh=None, bandwidth=0, matrix=None, sizeHint=None):
        pass

    matrix     = None
    numpyArray = property()
    _shape     = property()

    __array_priority__ = 100.0

    def __array_wrap(self, arr, context=None):
        if context is None:
            return arr
        else:
            return NotImplemented

    def copy(self):
        pass

    def __getitem__(self, index):
        pass

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
        pass

    def __add__(self, other):
        pass

    __radd__ = __add__

    def __iadd__(self, other):
        pass

    def __sub__(self, other):
        pass

    # Ask about this rsub
    def __rsub__(self, other):
        return -(__sub__(self, other))


    def __isub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __neg__(self):
        return self * -1

    def __pos__(self):
        return self

##     def __eq__(self,other):
## 	return self.matrix.__eq__(other.matrix)

##     def transpose(self):
##         pass

    def put(self, vector, id1, id2):
        pass

    def putDiagonal(self, vector):
        pass

    def take(self, id1, id2):
        pass

    def takeDiagonal(self):
        pass

    def addAt(self, vector, id1, id2):
        pass

    def addAtDiagonal(self, vector):
        pass

    def exportMmf(self, filename):
        pass

    @property
    def CSR(self):
        pass

    @property
    def LIL(self):
        pass

    @property
    def T(self):
        pass

##     def __array__(self):
##      shape = self._shape
##      indices = numerix.indices(shape)
##         numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
##      return numerix.reshape(numMatrix, shape)

class _Mesh2Matrix(object):
    def __init__(self, mesh, numberOfEquations=1, numberOfVariables=1,
                 orderer=lambda IDs: IDs):
        """Creates a mapping between mesh cells and matrix rows and columns

        Parameters
        ----------
        mesh : ~fipy.meshes.mesh.Mesh
            The `Mesh` to be mapped to a matrix
        numberOfVariables : int
            The local columns of the matrix are determined by
            `numberOfVariables * len(mesh._localNonOverlappingCellIDs)`.
        numberOfEquations : int
            The local rows of the matrix are determined by
            `numberOfEquations * len(mesh._localNonOverlappingCellIDs)`.
        orderer : fn
            Function that changes the order of IDs to satisfy matrix
            construction requirements (Default returns IDs unchanged).
        """
        self.mesh = mesh
        self.numberOfVariables = numberOfVariables
        self.numberOfEquations = numberOfEquations
        self.orderer = orderer

    def _cellIDsToGlobalRowIDs(self, IDs):
        N = len(IDs)
        M = self.numberOfEquations
        return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.globalNumberOfCells).flatten()

    def _cellIDsToLocalRowIDs(self, IDs):
        M = self.numberOfEquations
        N = len(IDs)
        return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.numberOfCells).flatten()

    @property
    def _globalNonOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def _globalOverlappingRowIDs(self):
        return self._cellIDsToGlobalRowIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def _localNonOverlappingRowIDs(self):
        return self._cellIDsToLocalRowIDs(self.mesh._localNonOverlappingCellIDs)

    def _cellIDsToGlobalColIDs(self, IDs):
        N = len(IDs)
        M = self.numberOfVariables
        return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.globalNumberOfCells).flatten()

    def _cellIDsToLocalColIDs(self, IDs):
        M = self.numberOfVariables
        N = len(IDs)
        return (numerix.vstack([IDs] * M) + numerix.indices((M,N))[0] * self.mesh.numberOfCells).flatten()

    @property
    def _globalNonOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalNonOverlappingCellIDs)

    @property
    def _globalCommonColIDs(self):
        return list(range(0, self.numberOfVariables, self.mesh.globalNumberOfCells))

    @property
    def _globalOverlappingColIDs(self):
        return self._cellIDsToGlobalColIDs(self.mesh._globalOverlappingCellIDs)

    @property
    def _localOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localOverlappingCellIDs)

    @property
    def _localNonOverlappingColIDs(self):
        return self._cellIDsToLocalColIDs(self.mesh._localNonOverlappingCellIDs)

    def _getStencil(self, id1, id2, overlapping=False):
        raise NotImplementedError

    def _globalVectorAndIDs(self, vector, id1, id2, overlapping=False):
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
            self._ghosts_ = self.orderer(self._ghosts_)

        return self._ghosts_

class _RowMesh2Matrix(_Mesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        id1 = self._globalOverlappingRowIDs[id1]

        if overlapping:
            mask = numerix.ones(id1.shape, dtype=bool)
        else:
            mask = numerix.in1d(id1, self._globalNonOverlappingRowIDs)

        id1 = self.orderer(id1[mask])
        id2 = numerix.asarray(id2)[mask]

        return id1, id2, mask

class _ColMesh2Matrix(_Mesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        id2 = self._globalOverlappingColIDs[id2]

        if overlapping:
            mask = numerix.ones(id2.shape, dtype=bool)
        else:
            mask = numerix.in1d(id2, self._globalNonOverlappingColIDs)

        id1 = numerix.asarray(id1)[mask]
        id2 = self.orderer(id2[mask])

        return id1, id2, mask

class _RowColMesh2Matrix(_RowMesh2Matrix):
    def _getStencil(self, id1, id2, overlapping=False):
        id2 = self._globalOverlappingColIDs[id2]

        id1, id2, mask = super(_RowColMesh2Matrix, self)._getStencil(id1, id2, overlapping)

        id2 = self.orderer(id2)

        return id1, id2, mask
