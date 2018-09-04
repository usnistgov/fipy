#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "addOverFacesVariable.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

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
        s = (numerix.newaxis,) * (len(contributions.shape) - 2) + (slice(0,None,None),) + (slice(0,None,None),)

        faceContributions = contributions * self.mesh._cellToFaceOrientations[s]

        return numerix.tensordot(numerix.ones(faceContributions.shape[-2], 'd'),
                                 numerix.MA.filled(faceContributions, 0.), (0, -2)) / self.mesh.cellVolumes
