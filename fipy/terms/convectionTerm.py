#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "convectionTerm.py"
 #                                    created: 11/13/03 {11:39:03 AM} 
 #                                last update: 1/16/04 {11:47:14 AM} 
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-13 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fivol.terms.faceTerm import FaceTerm
from fivol.variables.vectorFaceVariable import VectorFaceVariable

class ConvectionTerm(FaceTerm):
    def __init__(self, convCoeff, mesh, boundaryConditions, diffusionTerm = None):

	if not isinstance(convCoeff, VectorFaceVariable):
	    convCoeff = VectorFaceVariable(mesh = mesh, value = convCoeff)

	self.diffusionTerm = diffusionTerm
	
	self.projectedCoefficients = convCoeff * mesh.getOrientedAreaProjections()
	
	self.coeff = self.projectedCoefficients.sum(1)
	
	if self.diffusionTerm == None:
	    diffCoeff = 1e-20
	else:
	    diffCoeff = self.diffusionTerm.getCoeff()
	    
	P = -self.coeff / diffCoeff
	
	alpha = self.Alpha(P)

	weight = {
	    'implicit':{
		'cell 1 diag':    -alpha,
		'cell 1 offdiag': -(1-alpha),
		'cell 2 diag':     (1-alpha),
		'cell 2 offdiag':  alpha
	    }
	}
	FaceTerm.__init__(self,weight,mesh,boundaryConditions)
