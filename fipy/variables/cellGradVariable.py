#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 4/2/05 {7:30:47 PM} 
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
 
import Numeric, MA

from fipy.variables.vectorCellVariable import VectorCellVariable
import fipy.tools.array as array
from fipy.tools.inline import inline
from fipy.variables.faceGradContributionsVariable import FaceGradContributions


class CellGradVariable(VectorCellVariable):
    def __init__(self, var, name = ''):
	VectorCellVariable.__init__(self, mesh = var.getMesh(), name = name)
	self.var = self._requires(var)
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
            ids = Numeric.array(MA.filled(ids, 0)),
            orientations = Numeric.array(MA.filled(orientations, 0)),
            volumes = Numeric.array(volumes),
            areaProj = Numeric.array(self.mesh._getAreaProjections()),
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
        from fipy.tools.array import MAtake
	contributions = MAtake(self.faceGradientContributions[:],ids.flat)

        contributions = array.reshape(contributions, (N, M, self.mesh.getDim()))
        orientations = array.reshape(orientations, (N, M, 1))
	grad = Numeric.array(array.sum(orientations * contributions, 1))

	grad = grad / volumes[:,Numeric.NewAxis]

	self.value = grad

    def _calcValue(self):
	N = self.mesh.getNumberOfCells()
	M = self.mesh._getMaxFacesPerCell()
	
	ids = self.mesh._getCellFaceIDs()

	orientations = self.mesh._getCellFaceOrientations()
	volumes = self.mesh.getCellVolumes()

	inline.optionalInline(self._calcValueIn, self._calcValuePy, N, M, ids, orientations, volumes)
