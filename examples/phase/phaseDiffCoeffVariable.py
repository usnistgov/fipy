## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "phaseDiffCoeffVariable.py"
 #                                    created: 12/8/03 {6:05:05 PM} 
 #                                last update: 12/11/03 {12:16:34 AM} 
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

from variables.faceVariable import FaceVariable
import Numeric

class PhaseDiffCoeffVariable(FaceVariable):
    def __init__(self,mesh, parameters):
	FaceVariable.__init__(self, mesh)
	self.parameters = parameters
	self.phi = self.requires(parameters['phi'])
	self.m = self.requires(parameters['mPhi'])
	
    def calcValue(self):
# 	phi = self.parameters['phi']
# 	m = self.parameters['mPhi']
	phi = self.phi
	m = self.m

	alpha = self.parameters['alpha']
	N = self.parameters['symmetry']
	c2 = self.parameters['anisotropy']

	dphi = phi.getFaceGrad()
	
	z = Numeric.arctan2(dphi[:,1],dphi[:,0])
	z = N * (z - self.parameters['theta'].getFaceValue())
	z = Numeric.tan(z / 2.)
	z = z * z
	z = (1. - z) / (1. + z)
	z = (1.+ c2 * z)
	
	self.value = alpha**2 * z * z
	
