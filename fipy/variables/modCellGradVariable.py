from __future__ import division
from __future__ import unicode_literals
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

            ITEM(val, i, vec) /= ITEM(volumes, i, NULL);
            ITEM(val, i, vec) = mod(ITEM(val, i, vec) * gridSpacing[vec[0]]) /  gridSpacing[vec[0]];
        """, val = val,
            ids = numerix.array(ids),
            orientations = numerix.array(orientations),
            volumes = numerix.array(volumes),
            areaProj = numerix.array(self.mesh._areaProjections),
            faceValues = numerix.array(self.var.arithmeticFaceValue),
            M = M,
            ni = N,
            gridSpacing = numerix.array(self.mesh._meshSpacing),
            shape=numerix.array(numerix.shape(val)))

        return self._makeValue(value = val)

    def _calcValueNoInline(self, N, M, ids, orientations, volumes):
        value = _GaussCellGradVariable._calcValueNoInline(self, N, M, ids, orientations, volumes)
        gridSpacing = self.mesh._meshSpacing
        return self.modPy(value * gridSpacing) / gridSpacing
