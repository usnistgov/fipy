from __future__ import unicode_literals
from fipy.tools import numerix

__all__ = ["OffsetSparseMatrix"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

def OffsetSparseMatrix(SparseMatrix, numberOfVariables, numberOfEquations):
    """
    Used in binary terms. `equationIndex` and `varIndex` need to be set statically before instantiation.
    """

    class OffsetSparseMatrixClass(SparseMatrix):
        equationIndex = 0
        varIndex = 0

        def __init__(self, mesh, bandwidth=0, sizeHint=None,
                     numberOfVariables=numberOfVariables, numberOfEquations=numberOfEquations):
            SparseMatrix.__init__(self, mesh=mesh, bandwidth=bandwidth, sizeHint=sizeHint,
                                  numberOfVariables=numberOfVariables, numberOfEquations=numberOfEquations)

        def put(self, vector, id1, id2):
            SparseMatrix.put(self, vector, id1 + self.mesh.numberOfCells * self.equationIndex, id2 + self.mesh.numberOfCells * self.varIndex)

        def addAt(self, vector, id1, id2):
            SparseMatrix.addAt(self, vector, id1 + self.mesh.numberOfCells * self.equationIndex, id2 + self.mesh.numberOfCells * self.varIndex)

        def addAtDiagonal(self, vector):
            if type(vector) in [type(1), type(1.)]:
                tmp = numerix.zeros((self.mesh.numberOfCells,), 'd')
                tmp[:] = vector
                SparseMatrix.addAtDiagonal(self, tmp)
            else:
                SparseMatrix.addAtDiagonal(self, vector)

    return OffsetSparseMatrixClass
