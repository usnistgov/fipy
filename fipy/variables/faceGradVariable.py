#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "faceGradVariable.py"
 #                                    created: 12/18/03 {2:52:12 PM} 
 #                                last update: 1/29/04 {6:52:32 PM}
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
 # protection and is in the public domain.  PFM is an experimental
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

from fivol.variables.vectorFaceVariable import VectorFaceVariable
import fivol.tools.array as array
import fivol.inline.inline as inline

class FaceGradVariable(VectorFaceVariable):
    def __init__(self, var, mod = None):
	VectorFaceVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)
	if mod is None:
	    self.mod = lambda argument: argument
	else:
	    self.mod = mod
	    
    def _calcValuePy(self, normals, cellGrad, id1, id2, tangents1, tangents2, N):
	grad1 = array.take(cellGrad,id1)
	grad2 = array.take(cellGrad,id2)
	t1grad1 = array.sum(tangents1*grad1,1)
	t1grad2 = array.sum(tangents1*grad2,1)
	t2grad1 = array.sum(tangents2*grad1,1)
	t2grad2 = array.sum(tangents2*grad2,1)
	
	T1 = (t1grad1 + t1grad2) / 2.
	T2 = (t2grad1 + t2grad2) / 2.
	
	N = N[:,Numeric.NewAxis]
	T1 = T1[:,Numeric.NewAxis]
	T2 = T2[:,Numeric.NewAxis]

	self.value = normals * N + tangents1 * T1 + tangents2 * T2

    def _calcValueIn(self, normals, cellGrad, id1, id2, tangents1, tangents2, N):
	inline.runInlineLoop1("""
	int j;
	double t1grad1, t1grad2, t2grad1, t2grad2;
	
	t1grad1 = t1grad2 = t2grad1 = t2grad2 = 0;
	for (j = 0; j < nj; j++) {
	    t1grad1 += tangents1(i,j) * cellGrad(id1(i),j);
	    t1grad2 += tangents1(i,j) * cellGrad(id2(i),j);
	    t2grad1 += tangents2(i,j) * cellGrad(id1(i),j);
	    t2grad2 += tangents2(i,j) * cellGrad(id2(i),j);
	}
	
	for (j = 0; j < nj; j++) {
	    val(i,j) = normals(i,j) * N(i);
	    val(i,j) += tangents1(i,j) * (t1grad1 + t1grad2) / 2.;
	    val(i,j) += tangents2(i,j) * (t2grad1 + t2grad2) / 2.;
	}
	""", 
	tangents1 = tangents1, tangents2 = tangents2,
	cellGrad = cellGrad, normals = normals, N = N, 
	id1 = id1, id2 = id2,
	val = self.value.value,
	ni = tangents1.shape[0], nj = tangents1.shape[1], nk = 0)

## 	t1grad1 = self.t1grad1, t1grad2 = self.t1grad2, 
## 	t2grad1 = self.t2grad1,  t2grad2 = self.t2grad2, 
	
    def calcValue(self):        
	dAP = self.mesh.getCellDistances()
	id1, id2 = self.mesh.getAdjacentCellIDs()
	N = self.mod(array.take(self.var,id2) - array.take(self.var,id1))/dAP
	
	normals = self.mesh.getOrientedFaceNormals()
	
	tangents1 = self.mesh.getFaceTangents1()
	tangents2 = self.mesh.getFaceTangents2()
	cellGrad = self.var.getGrad().getNumericValue()
	
	inline.optionalInline(self._calcValueIn, self._calcValuePy, normals, cellGrad, id1, id2, tangents1, tangents2, N)
	
##    def calcValue(self):
##        id1, id2, N, normals = self.calcValue1()
##        T1, T2, tangents1, tangents2 = self.calcValue2(id1, id2)
##        self.calcValue3(N, normals, T1, T2, tangents1, tangents2)

##    def calcValue1(self):
##	dAP = self.mesh.getCellDistances()
##	id1, id2 = self.mesh.getAdjacentCellIDs()
##	value = self.var[:]
##	N = self.mod(Numeric.take(value, id2) - Numeric.take(value, id1))/dAP
##	normals = self.mesh.getOrientedFaceNormals()

##        return id1, id2, N, normals
    
##    def calcValue2(self, id1, id2):
##	tangents1 = self.mesh.getFaceTangents1()
##	tangents2 = self.mesh.getFaceTangents2()
##	cellGrad = self.var.getGrad()
##	grad1 = Numeric.take(cellGrad[:], id1)
##	grad2 = Numeric.take(cellGrad[:], id2)
##	t1grad1 = Numeric.sum(tangents1*grad1,1)
##	t1grad2 = Numeric.sum(tangents1*grad2,1)
##	t2grad1 = Numeric.sum(tangents2*grad1,1)
##	t2grad2 = Numeric.sum(tangents2*grad2,1)
##	T1 = (t1grad1 + t1grad2) / 2.
##	T2 = (t2grad1 + t2grad2) / 2.

##        return T1, T2, tangents1, tangents2

##    def calcValue3(self, N, normals, T1, T2, tangents1, tangents2):
	
##	N = Numeric.reshape(N, (len(normals),1)) 
##	T1 = Numeric.reshape(T1, (len(normals),1)) 
##	T2 = Numeric.reshape(T2, (len(normals),1)) 

##	self.value = normals * N + tangents1 * T1 + tangents2 * T2
