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

##     def __array__(self):
##      shape = self._shape
##      indices = numerix.indices(shape)
##         numMatrix = self.take(indices[0].ravel(), indices[1].ravel())
##      return numerix.reshape(numMatrix, shape)

class _Mesh2Matrix(object):
    _bodies = None
    _ghosts = None

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
            mask = numerix.in1d(id1, globalNonOverlappihgIDs)

        id1 = self.orderer(id1[mask])
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
            self._bodies = numerix.in1d(self.mesh._globalOverlappingCellIDs,
                                        self.mesh._globalNonOverlappingCellIDs)
        return self._bodies

    @property
    def ghosts(self):
        if self._ghosts is None:
            self._ghosts = self.mesh._globalOverlappingCellIDs[~self.bodies]
            self._ghosts = self._cellIDsToGlobalRowIDs(self._ghosts)
            self._ghosts = self.orderer(self._ghosts)

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

        id2 = self.orderer(id2)

        return id1, id2, mask
