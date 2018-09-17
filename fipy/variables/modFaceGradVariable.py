#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "faceGradVariable.py"
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

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.variables.faceGradVariable import _FaceGradVariable
from fipy.tools import inline
from fipy.tools import numerix

class _ModFaceGradVariable(_FaceGradVariable):
    def __init__(self, var, modIn):
        _FaceGradVariable.__init__(self, var)
        self.modIn = modIn

    if inline.doInline:
        def _calcValue(self):

            id1, id2 = self.mesh._adjacentCellIDs

            tangents1 = self.mesh._faceTangents1
            tangents2 = self.mesh._faceTangents2

            val = self._array.copy()

            inline._runIterateElementInline(self.modIn + """
            int j;
            double t1grad1, t1grad2, t2grad1, t2grad2, N;
            int ID1 = ITEM(id1, i, NULL);
            int ID2 = ITEM(id2, i, NULL);
            N = mod(ITEM(var, ID2, NULL) - ITEM(var, ID1, NULL)) / ITEM(dAP, i, NULL);

            t1grad1 = t1grad2 = t2grad1 = t2grad2 = 0.;

            t1grad1 += ITEM(tangents1, i, vec) * ITEM(cellGrad, ID1, vec);
            t1grad2 += ITEM(tangents1, i, vec) * ITEM(cellGrad, ID2, vec);
            t2grad1 += ITEM(tangents2, i, vec) * ITEM(cellGrad, ID1, vec);
            t2grad2 += ITEM(tangents2, i, vec) * ITEM(cellGrad, ID2, vec);

            ITEM(val, i, vec) =  ITEM(normals, i, vec) * N;
            ITEM(val, i, vec) += ITEM(tangents1, i, vec) * (t1grad1 + t1grad2) / 2.;
            ITEM(val, i, vec) += ITEM(tangents2, i, vec) * (t2grad1 + t2grad2) / 2.;

            """,tangents1 = tangents1,
                tangents2 = tangents2,
                cellGrad = self.var.grad.numericValue,
                normals = numerix.array(self.mesh._orientedFaceNormals),
                id1 = numerix.array(id1),
                id2 = numerix.array(id2),
                dAP = numerix.array(self.mesh._cellDistances),
                var = self.var.numericValue,
                val = val,
                ni = tangents1.shape[1],
                shape=numerix.array(numerix.shape(val)))

            return self._makeValue(value = val)
    else:
        def _calcValue(self):
            dAP = self.mesh._cellDistances
            id1, id2 = self.mesh._adjacentCellIDs
            N = (numerix.take(self.var,id2) - numerix.take(self.var,id1)) / dAP
            normals = self.mesh._orientedFaceNormals

            tangents1 = self.mesh._faceTangents1
            tangents2 = self.mesh._faceTangents2
            cellGrad = self.var.grad.numericValue

            grad1 = numerix.take(cellGrad, id1, axis=1)
            grad2 = numerix.take(cellGrad, id2, axis=1)
            t1grad1 = numerix.sum(tangents1*grad1,0)
            t1grad2 = numerix.sum(tangents1*grad2,0)
            t2grad1 = numerix.sum(tangents2*grad1,0)
            t2grad2 = numerix.sum(tangents2*grad2,0)

            T1 = (t1grad1 + t1grad2) / 2.
            T2 = (t2grad1 + t2grad2) / 2.

            return normals * N + tangents1 * T1 + tangents2 * T2
