#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "cellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 2/3/04 {3:31:53 PM} 
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

from fivol.variables.vectorCellVariable import VectorCellVariable
from fivol.tools import array
from fivol.inline import inline
from fivol.variables.faceGradContributionsVariable import FaceGradContributions

class CellGradVariable(VectorCellVariable):
    def __init__(self, var):
	VectorCellVariable.__init__(self, var.getMesh())
	self.var = self.requires(var)
        self.faceGradientContributions = FaceGradContributions(self.var)
        
    def _calcValueIn(self, N, M, ids, orientations, volumes):
	inline.runInlineLoop2("""
	    val(i,j) = 0.;
	    
	    int k;
	    for (k = 0; k < nk; k++) {
		val(i,j) += orientations(i,k) * faceGradientContributions(ids(i,k), j);
	    }
		
	    val(i,j) /= volumes(i);
	""",
	val = self.value.value, ids = ids, orientations = orientations, volumes = volumes,
	faceGradientContributions = self.faceGradientContributions.getNumericValue(),
	ni = N, nj = self.mesh.getDim(), nk = M
	)
	    
    def _calcValuePy(self, N, M, ids, orientations, volumes):
## 	print "self.var:",self.var.__class__
## 	print 'getFaceValue:',self.var.getFaceValue()[:]
## 	raw_input()

	contributions = array.take(self.faceGradientContributions[:],ids.flat)
	contributions = contributions.reshape((N,M,self.mesh.getDim()))

	grad = (orientations*contributions).sum(1)
	
	grad = grad / volumes[:,Numeric.NewAxis]

	self.value = grad
	    
    def calcValue(self):
	N = len(self.mesh.getCells())
	M = self.mesh.getMaxFacesPerCell()
	
	ids = self.mesh.getCellFaceIDs()

	orientations = self.mesh.getCellFaceOrientations()
	volumes = self.mesh.getCellVolumes()

	inline.optionalInline(self._calcValueIn, self._calcValuePy, N, M, ids, orientations, volumes)
