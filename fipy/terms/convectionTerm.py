"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "convectionTerm.py"
                                   created: 11/13/03 {11:39:03 AM} 
                               last update: 12/19/03 {4:09:33 PM} 
 Author: Jonathan Guyer
 E-mail: guyer@nist.gov
 Author: Daniel Wheeler
 E-mail: daniel.wheeler@nist.gov
   mail: NIST
    www: http://ctcms.nist.gov
 
========================================================================
This software was developed at the National Institute of Standards
and Technology by employees of the Federal Government in the course
of their official duties.  Pursuant to title 17 Section 105 of the
United States Code this software is not subject to copyright
protection and is in the public domain.  PFM is an experimental
system.  NIST assumes no responsibility whatsoever for its use by
other parties, and makes no guarantees, expressed or implied, about
its quality, reliability, or any other characteristic.  We would
appreciate acknowledgement if the software is used.

This software can be redistributed and/or modified freely
provided that any derivative works bear some notice that they are
derived from it, and any modified versions bear some notice that
they have been modified.
========================================================================
 
 Description: 

 History

 modified   by  rev reason
 ---------- --- --- -----------
 2003-11-13 JEG 1.0 original
###################################################################
"""

from faceTerm import FaceTerm
import Numeric

class ConvectionTerm(FaceTerm):
    def __init__(self, convCoeff, mesh, boundaryConditions, diffusionTerm = None):
	weight = {
	    'implicit':{
		'cell 1 diag':    -0.5, 
		'cell 1 offdiag':  0.5, 
		'cell 2 diag':     0.5, 
		'cell 2 offdiag': -0.5
	    }
	}
	FaceTerm.__init__(self,weight,mesh,boundaryConditions)
	if type(convCoeff) in [type(1), type(1.), type(Numeric.array((1,)))]:
	    self.convCoeff = Numeric.array(convCoeff)[Numeric.NewAxis,:]
	else:
	    self.convCoeff = convCoeff
	self.diffusionTerm = diffusionTerm
	
    def calculateAlpha(self, P):
	pass
	
    def calculateCoeffGeom(self,dt):
	areas = self.mesh.getOrientedAreaProjections()
	self.coeff = Numeric.sum((self.convCoeff * areas)[:], 1)
	if self.diffusionTerm == None:
	    diffCoeff = 1e-20
	else:
	    diffCoeff = self.diffusionTerm.getCoeff()
	
	P = -self.coeff / diffCoeff
	
	alpha = self.calculateAlpha(P)
	
	self.weight['implicit']['cell 1 diag'] = -alpha
	self.weight['implicit']['cell 1 offdiag'] = -(1-alpha)
	self.weight['implicit']['cell 2 diag'] = (1-alpha)
	self.weight['implicit']['cell 2 offdiag'] = alpha
