#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modCellGradVariable.py"
 #                                    created: 12/18/03 {2:28:00 PM} 
 #                                last update: 9/3/04 {10:33:20 PM} 
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

from fipy.variables.cellGradVariable import _CellGradVariable
from fipy.tools.inline import inline

class _ModCellGradVariable(_CellGradVariable):
    def __init__(self, var, modIn, modPy):
        _CellGradVariable.__init__(self, var)
        self.modIn = modIn
        self.modPy = modPy
        
    def _calcValueIn(self, N, M, ids, orientations, volumes):
        
	inline._runInlineLoop2(self.modIn + """
	    val(i,j) = 0.;
	    
	    int k;
            
	    for (k = 0; k < nk; k++) {
		val(i, j) += orientations(i, k) * areaProj(ids(i, k), j) * faceValues(ids(i, k));
	    }
		
	    val(i, j) /= volumes(i);
            val(i, j) = mod(val(i,j) * gridSpacing(j)) /  gridSpacing(j);
	""",
	val = self._getArray(),
        ids = Numeric.array(ids),
        orientations = Numeric.array(orientations),
        volumes = Numeric.array(volumes),
        areaProj = Numeric.array(self.mesh._getAreaProjections()),
        faceValues = Numeric.array(self.var.getArithmeticFaceValue()),
	ni = N, nj = self.mesh.getDim(), nk = M,
        gridSpacing = Numeric.array(self.mesh._getMeshSpacing()))

    def _calcValuePy(self, N, M, ids, orientations, volumes):
        _CellGradVariable._calcValuePy(self, N, M, ids, orientations, volumes)
        gridSpacing = self.mesh._getMeshSpacing()
	self.value = self.modPy(self.value * gridSpacing) / gridSpacing

