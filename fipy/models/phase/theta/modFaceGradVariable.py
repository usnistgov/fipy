#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceGradVariable.py"
 #                                    created: 12/18/03 {2:52:12 PM} 
 #                                last update: 9/3/04 {10:41:39 PM}
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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

import Numeric

from fipy.variables.faceGradVariable import FaceGradVariable
from fipy.tools.inline import inline

class ModFaceGradVariable(FaceGradVariable):
    def __init__(self, var, modIn):
	FaceGradVariable.__init__(self, var)
        self.modIn = modIn
        
    def _calcValueInline(self):

	id1, id2 = self.mesh.getAdjacentCellIDs()
	
	tangents1 = self.mesh.getFaceTangents1()
	tangents2 = self.mesh.getFaceTangents2()

	inline.runInline(self.modIn + """
        int i;
        for(i = 0; i < ni; i++)
        {
            int j;
            double t1grad1, t1grad2, t2grad1, t2grad2, N;

            N = mod(var(id2(i)) - var(id1(i))) / dAP(i);

	    t1grad1 = t1grad2 = t2grad1 = t2grad2 = 0.;
            
	    for (j = 0; j < nj; j++) {
		t1grad1 += tangents1(i,j) * cellGrad(id1(i),j);
		t1grad2 += tangents1(i,j) * cellGrad(id2(i),j);
		t2grad1 += tangents2(i,j) * cellGrad(id1(i),j);
		t2grad2 += tangents2(i,j) * cellGrad(id2(i),j);
	    }
	    
	    for (j = 0; j < nj; j++) {
		val(i,j) = normals(i,j) * N;
		val(i,j) += tangents1(i,j) * (t1grad1 + t1grad2) / 2.;
		val(i,j) += tangents2(i,j) * (t2grad1 + t2grad2) / 2.;
	    }
        }
        """,tangents1 = tangents1,
            tangents2 = tangents2,
            cellGrad = self.var.getGrad().getNumericValue(),
            normals = Numeric.array(self.mesh.getOrientedFaceNormals()),
            id1 = Numeric.array(id1),
            id2 = Numeric.array(id2),
            dAP = Numeric.array(self.mesh.getCellDistances()),
            var = self.var.getNumericValue(),
            val = self._getArray(),
            ni = tangents1.shape[0],
            nj = tangents1.shape[1])
