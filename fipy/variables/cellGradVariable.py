#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 7/24/04 {9:02:02 AM} 
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

from fipy.variables.vectorCellVariable import VectorCellVariable
import fipy.tools.array as array
from fipy.tools.inline import inline
from fipy.variables.faceGradContributionsVariable import FaceGradContributions

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
		val(i, j) += orientations(i, k) * areaProj(ids(i, k), j) * faceValues(ids(i, k));
	    }
		
	    val(i, j) /= volumes(i);
	""",val = self._getArray(),
            ids = Numeric.array(ids),
            orientations = Numeric.array(orientations),
            volumes = Numeric.array(volumes),
            areaProj = Numeric.array(self.mesh.getAreaProjections()),
            faceValues = Numeric.array(self.var.getArithmeticFaceValue()),
	    ni = N, nj = self.mesh.getDim(), nk = M)
        
##    def _calcValueIn(self, N, M, ids, orientations, volumes):
##	inline.runInlineLoop2("""
##	    val(i,j) = 0.;
	    
##	    int k;
##	    for (k = 0; k < nk; k++) {
##		val(i,j) += orientations(i,k) * faceGradientContributions(ids(i,k), j);
##	    }
		
##	    val(i,j) /= volumes(i);
##	""",
##	val = self._getArray(), ids = ids, orientations = orientations, volumes = volumes,
##	faceGradientContributions = self.faceGradientContributions.getNumericValue(),
##	ni = N, nj = self.mesh.getDim(), nk = M
##	)
	    
    def _calcValuePy(self, N, M, ids, orientations, volumes):
	contributions = array.take(self.faceGradientContributions[:],ids.flat)

##        contributions = contributions.reshape((N,M,self.mesh.getDim()))
        contributions = array.reshape(contributions, (N, M, self.mesh.getDim()))
        orientations = array.reshape(orientations, (N, M, 1))
	grad = array.sum(orientations * contributions, 1)
##        grad = (orientations*contributions).sum(1)

	grad = grad / volumes[:,Numeric.NewAxis]

	self.value = grad

    def _calcValue(self):
	N = self.mesh.getNumberOfCells()
	M = self.mesh.getMaxFacesPerCell()
	
	ids = self.mesh.getCellFaceIDs()

	orientations = self.mesh.getCellFaceOrientations()
	volumes = self.mesh.getCellVolumes()

	inline.optionalInline(self._calcValueIn, self._calcValuePy, N, M, ids, orientations, volumes)
