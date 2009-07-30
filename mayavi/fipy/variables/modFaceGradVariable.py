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
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

from fipy.variables.faceGradVariable import _FaceGradVariable
from fipy.tools import inline
from fipy.tools import numerix

class _ModFaceGradVariable(_FaceGradVariable):
    def __init__(self, var, modIn):
        _FaceGradVariable.__init__(self, var)
        self.modIn = modIn
        
    def _calcValueInline(self):

        id1, id2 = self.mesh._getAdjacentCellIDs()
        
        tangents1 = self.mesh._getFaceTangents1()
        tangents2 = self.mesh._getFaceTangents2()
 
        val = self._getArray().copy()

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
            cellGrad = self.var.getGrad().getNumericValue(),
            normals = numerix.array(self.mesh._getOrientedFaceNormals()),
            id1 = numerix.array(id1),
            id2 = numerix.array(id2),
            dAP = numerix.array(self.mesh._getCellDistances()),
            var = self.var.getNumericValue(),
            val = val,
            ni = tangents1.shape[1],
            shape=numerix.array(numerix.shape(val)))
            
        return self._makeValue(value = val)
##         return self._makeValue(value = val, unit = self.getUnit())
