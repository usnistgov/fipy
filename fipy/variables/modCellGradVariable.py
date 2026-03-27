__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.gaussCellGradVariable import _GaussCellGradVariable
from fipy.tools import inline
from fipy.tools import numerix

class _ModCellGradVariable(_GaussCellGradVariable):
    def __init__(self, var, modIn, modPy):
        _GaussCellGradVariable.__init__(self, var)
        self.modIn = modIn
        self.modPy = modPy


    def _calcValueInline(self, N, M, ids, orientations, volumes):
        val = self._array.copy()

        inline._runIterateElementInline(self.modIn + """
            ITEM(val, i, vec) = 0.;

            int k;
            for (k = 0; k < M; k++) {
                int id = ITEM(ids, i, &k);
                ITEM(val, i, vec) += ITEM(orientations, i, &k) * ITEM(areaProj, id, vec) * ITEM(faceValues, id, NULL);
            }

            ITEM(val, i, vec) = mod(ITEM(val, i, vec) * ITEM(avgFaceSize, i, vec) /  * ITEM(avgFaceSize, i, vec);
            ITEM(val, i, vec) /= ITEM(volumes, i, NULL);
        """, val = val,
            ids = numerix.array(ids),
            orientations = numerix.array(orientations),
            volumes = numerix.array(volumes),
            areaProj = numerix.array(self.mesh._areaProjections),
            faceValues = numerix.array(self.var.arithmeticFaceValue),
            M = M,
            ni = N,
            avgFaceSize = self.mesh._averageFaceSizePerCell,
            shape=numerix.array(numerix.shape(val)))

        return self._makeValue(value = val)

    def _gradAreaPerCell(self, N, M, ids, orientations, volumes):
        gradArea = super(_ModCellGradVariable, self)._gradAreaPerCell(N, M, ids, orientations, volumes)
        avgFaceSize = self.mesh._averageFaceSizePerCell
        return self.modPy(gradArea / avgFaceSize) * avgFaceSize
