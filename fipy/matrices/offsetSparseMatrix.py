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

        def __init__(self, mesh, nonZerosPerRow=1, exactNonZeros=False,
                     numberOfVariables=numberOfVariables,
                     numberOfEquations=numberOfEquations):
            if hasattr(nonZerosPerRow, "__iter__"):
                # nonZerosPerRow is an iterable for each row.
                # need to pad rows for other equations with zeros.
                # can't compare to collections.abc.Iterable because PySparse.
                before = self.equationIndex
                after = numberOfEquations - self.equationIndex - 1
                N = len(nonZerosPerRow)
                nonZerosPerRow = numerix.concatenate([[0] * N * before,
                                                      nonZerosPerRow,
                                                      [0] * N * after])
            else:
                nonZerosPerRow //= numberOfEquations
            SparseMatrix.__init__(self,
                                  mesh=mesh,
                                  nonZerosPerRow=nonZerosPerRow,
                                  exactNonZeros=exactNonZeros,
                                  numberOfVariables=numberOfVariables,
                                  numberOfEquations=numberOfEquations)

        def put(self, vector, id1, id2):
            N = self.mesh.numberOfCells
            SparseMatrix.put(self,
                             vector=vector,
                             id1=id1 + N * self.equationIndex,
                             id2=id2 + N * self.varIndex)

        def addAt(self, vector, id1, id2):
            N = self.mesh.numberOfCells
            SparseMatrix.addAt(self,
                               vector=vector,
                               id1=id1 + N * self.equationIndex,
                               id2=id2 + N * self.varIndex)

        def addAtDiagonal(self, vector):
            if type(vector) in [type(1), type(1.)]:
                tmp = numerix.zeros((self.mesh.numberOfCells,), 'd')
                tmp[:] = vector
                SparseMatrix.addAtDiagonal(self, tmp)
            else:
                SparseMatrix.addAtDiagonal(self, vector)

    return OffsetSparseMatrixClass
