#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceGradVariable.py"
 #                                    created: 12/18/03 {2:52:12 PM} 
 #                                last update: 1/3/07 {3:26:15 PM}
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

from fipy.variables.faceVariable import FaceVariable
from fipy.tools import numerix
from fipy.tools.inline import inline

class _FaceGradVariable(FaceVariable):
    def __init__(self, var):
        FaceVariable.__init__(self, mesh=var.getMesh(), rank=var.getRank() + 1)
        self.var = self._requires(var)

    def _calcValue(self):        
        return inline._optionalInline(self._calcValueInline, self._calcValuePy)
    
    def _calcValuePy(self):
        dAP = self.mesh._getCellDistances()
        id1, id2 = self.mesh._getAdjacentCellIDs()
##      N = self.mod(numerix.take(self.var,id2) - numerix.take(self.var,id1)) / dAP
        N = (numerix.take(self.var,id2) - numerix.take(self.var,id1)) / dAP
        normals = self.mesh._getOrientedFaceNormals()
        
        tangents1 = self.mesh._getFaceTangents1()
        tangents2 = self.mesh._getFaceTangents2()
        cellGrad = self.var.getGrad().getNumericValue()
        
        grad1 = numerix.take(cellGrad, id1, axis=1)
        grad2 = numerix.take(cellGrad, id2, axis=1)
        t1grad1 = numerix.sum(tangents1*grad1,0)
        t1grad2 = numerix.sum(tangents1*grad2,0)
        t2grad1 = numerix.sum(tangents2*grad1,0)
        t2grad2 = numerix.sum(tangents2*grad2,0)
        
        T1 = (t1grad1 + t1grad2) / 2.
        T2 = (t2grad1 + t2grad2) / 2.
        
        return normals * N + tangents1 * T1 + tangents2 * T2

    def _calcValueInline(self):

        id1, id2 = self.mesh._getAdjacentCellIDs()
        
        tangents1 = self.mesh._getFaceTangents1()
        tangents2 = self.mesh._getFaceTangents2()
 
        val = self._getArray().copy()

        inline._runInline("""
            int j;
            double t1grad1, t1grad2, t2grad1, t2grad2, N;
            int ID1 = id1[i];
            int ID2 = id2[i];
            
            N = (var[ID2] - var[ID1]) / dAP[i];

            t1grad1 = t1grad2 = t2grad1 = t2grad2 = 0.;
            
            for (j = 0; j < NJ; j++) {
                int arrayID = i * NJ + j;
                t1grad1 += tangents1[arrayID] * cellGrad[ID1 * NJ + j];
                t1grad2 += tangents1[arrayID] * cellGrad[ID2 * NJ + j];
                t2grad1 += tangents2[arrayID] * cellGrad[ID1 * NJ + j];
                t2grad2 += tangents2[arrayID] * cellGrad[ID2 * NJ + j];
            }
            
            for (j = 0; j < NJ; j++) {
                int arrayID = i * NJ + j;
                val[arrayID] = normals[arrayID] * N;
                val[arrayID] += tangents1[arrayID] * (t1grad1 + t1grad2) / 2.;
                val[arrayID] += tangents2[arrayID] * (t2grad1 + t2grad2) / 2.;
            }
        """,tangents1 = tangents1,
            tangents2 = tangents2,
            cellGrad = self.var.getGrad().getNumericValue(),
            normals = self.mesh._getOrientedFaceNormals(),
            id1 = id1,
            id2 = id2,
            dAP = numerix.array(self.mesh._getCellDistances()),
            var = self.var.getNumericValue(),
            val = val,
            ni = tangents1.shape[0],
            NJ = tangents1.shape[1])
            
        return self._makeValue(value = val)
##         return self._makeValue(value = val, unit = self.getUnit())

    
