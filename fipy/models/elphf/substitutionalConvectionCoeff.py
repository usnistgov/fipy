## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "substitutionalConvectionCoeff.py"
 #                                    created: 12/9/03 {3:32:09 PM} 
 #                                last update: 12/9/03 {4:12:44 PM} 
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

from variables.vectorFaceVariable import VectorFaceVariable
from substitutionalSumVariable import SubstitutionalSumVariable
import Numeric

class SubstitutionalConvectionCoeff(VectorFaceVariable):
    def __init__(self,mesh,diffusivity,Cj,substitutionals):
	VectorFaceVariable.__init__(self,mesh,Cj.name + "_convection")
	self.Dj = diffusivity
	self.substitutionalSum = SubstitutionalSumVariable(
	    mesh = mesh, 
	    Cj = Cj, 
	    substitutionals = substitutionals)
	    
    def getValue(self):
	num = self.Dj * self.substitutionalSum.getFaceGrad() 
	den = (1. - self.substitutionalSum.getFaceValue())
	return num / Numeric.reshape(den,(len(num),1))
