#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "faceGradVariable.py"
 #                                    created: 12/18/03 {2:52:12 PM} 
 #                                last update: 7/12/05 {1:08:19 PM}
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

import Numeric

from fipy.variables.vectorFaceVariable import VectorFaceVariable
from fipy.tools import numerix
from fipy.tools.inline import inline

class _FaceGradVariable(VectorFaceVariable):
    def __init__(self, var):
	VectorFaceVariable.__init__(self, var.getMesh())
	self.var = self._requires(var)

    def _calcValue(self):        
	inline._optionalInline(self._calcValueInline, self._calcValuePy)
    
    def _calcValuePy(self):
    
        dAP = self.mesh._getCellDistances()
	id1, id2 = self.mesh._getAdjacentCellIDs()
##	N = self.mod(numerix.take(self.var,id2) - numerix.take(self.var,id1)) / dAP
	N = (numerix.take(self.var,id2) - numerix.take(self.var,id1)) / dAP
	normals = self.mesh._getOrientedFaceNormals()
        
	tangents1 = self.mesh._getFaceTangents1()
	tangents2 = self.mesh._getFaceTangents2()
	cellGrad = self.var.getGrad().getNumericValue()
        
      
        
	grad1 = numerix.take(cellGrad,id1)
	grad2 = numerix.take(cellGrad,id2)
	t1grad1 = numerix.sum(tangents1*grad1,1)
	t1grad2 = numerix.sum(tangents1*grad2,1)
	t2grad1 = numerix.sum(tangents2*grad1,1)
	t2grad2 = numerix.sum(tangents2*grad2,1)
	
	T1 = (t1grad1 + t1grad2) / 2.
	T2 = (t2grad1 + t2grad2) / 2.
	
	N = N[:,Numeric.NewAxis]
	T1 = T1[:,Numeric.NewAxis]
	T2 = T2[:,Numeric.NewAxis]
        
	self.value = normals * N + tangents1 * T1 + tangents2 * T2

    def _calcValueInline(self):

	id1, id2 = self.mesh._getAdjacentCellIDs()
	
	tangents1 = self.mesh._getFaceTangents1()
	tangents2 = self.mesh._getFaceTangents2()

	inline._runInlineLoop1("""
            int j;
            double t1grad1, t1grad2, t2grad1, t2grad2, N;

	    N = (var(id2(i)) - var(id1(i))) / dAP(i);

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
        """,tangents1 = tangents1,
            tangents2 = tangents2,
            cellGrad = self.var.getGrad().getNumericValue(),
            normals = self.mesh._getOrientedFaceNormals(),
            id1 = id1,
            id2 = id2,
	    dAP = Numeric.array(self.mesh._getCellDistances()),
            var = self.var.getNumericValue(),
            val = self._getArray(),
            ni = tangents1.shape[0],
            nj = tangents1.shape[1])


    
