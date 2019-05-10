from __future__ import division
from __future__ import unicode_literals
__all__ = []

from fipy.tools import numerix
from fipy.tools import inline
from fipy.variables.cellVariable import CellVariable

class _AddOverFacesVariable(CellVariable):
    r"""surface integral of `self.faceVariable`, :math:`\phi_f`

    .. math:: \int_S \phi_f\,dS \approx \frac{\sum_f \phi_f A_f}{V_P}

    Returns
    -------
    integral : CellVariable
        volume-weighted sum
    """

    def __init__(self, faceVariable, mesh = None):
        if not mesh:
            mesh = faceVariable.mesh

        CellVariable.__init__(self, mesh, hasOld = 0, elementshape=faceVariable.shape[:-1])
        self.faceVariable = self._requires(faceVariable)

    def _calcValue(self):
        if inline.doInline and self.faceVariable.rank < 2:
            return self._calcValueInline()
        else:
            return self._calcValueNoInline()

    def _calcValueInline(self):

        NCells = self.mesh.numberOfCells
        ids = self.mesh.cellFaceIDs

        val = self._array.copy()

        inline._runInline("""
        int i;

        for(i = 0; i < numberOfCells; i++)
          {
          int j;
          value[i] = 0.;
          for(j = 0; j < numberOfCellFaces; j++)
            {
              // cellFaceIDs can be masked, which caused subtle and
              // unreproducible problems on OS X (who knows why not elsewhere)
              long id = ids[i + j * numberOfCells];
              if (id >= 0) {
                  value[i] += orientations[i + j * numberOfCells] * faceVariable[id];
              }
            }
            value[i] = value[i] / cellVolume[i];
          }
          """,
                          numberOfCellFaces = self.mesh._maxFacesPerCell,
                          numberOfCells = NCells,
                          faceVariable = self.faceVariable.numericValue,
                          ids = numerix.array(ids),
                          value = val,
                          orientations = numerix.array(self.mesh._cellToFaceOrientations),
                          cellVolume = numerix.array(self.mesh.cellVolumes))

        return self._makeValue(value = val)

    def _calcValueNoInline(self):
        ids = self.mesh.cellFaceIDs

        contributions = numerix.take(self.faceVariable, ids, axis=-1)

        # FIXME: numerix.MA.filled casts away dimensions
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0, None, None),) + (slice(0, None, None),)

        faceContributions = contributions * self.mesh._cellToFaceOrientations[s]

        return numerix.tensordot(numerix.ones(faceContributions.shape[-2], 'd'),
                                 numerix.MA.filled(faceContributions, 0.), (0, -2)) / self.mesh.cellVolumes
