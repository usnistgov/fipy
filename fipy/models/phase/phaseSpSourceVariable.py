## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "phaseSourceField.py"
 #                                    created: 12/8/03 {4:44:40 PM} 
 #                                last update: 12/9/03 {2:24:47 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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

from variables.cellVariable import CellVariable

class PhaseSpSourceVariable(CellVariable):
    def __init__(self,mesh,parameters):
	self.parameters = parameters
	CellVariable.__init__(self,name = "SpSource", mesh = mesh, hasOld = False)
	
    def getValue(self):
	phi = self.parameters['phi']
	m = self.parameters['mPhi']

	## driving force double well

	sp = m * (phi - (m < 0.))
    
	## theta source terms

	thetaMag = self.parameters['theta'].getOld().getGradMag()
	s = self.parameters['s']
	epsilon = self.parameters['epsilon']

	sp += (2*s + epsilon**2 * thetaMag) * thetaMag
	
	return sp
