"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "phaseScSourceVariable.py"
 #                                    created: 12/8/03 {4:44:40 PM} 
 #                                last update: 12/19/03 {12:30:42 PM} 
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
"""

from variables.cellVariable import CellVariable

class PhaseScSourceVariable(CellVariable):
    def __init__(self,mesh,parameters):
	CellVariable.__init__(self,name = "ScSource", mesh = mesh, hasOld = 0)
	self.phi = self.requires(parameters['phi'])
	self.m = self.requires(parameters['mPhi'])
	
    def calcValue(self):
	phi = self.phi
	m = self.m
	
	## driving force double well

	sc = (m > 0.) * m * phi
    
	## theta source terms

#         ## anisotropy
# 
# ##        z = Numeric.atan2(self.dphi[:,1],self.dphi[:,0]);
# ##        z = N * (z-self.theta);
# ##        z = tan(0.5 * z);
# ##        zsq = z * z;
# ##        b = (1-zsq) / (1+zsq);
# ##        db = -N * 2. *z / (1+zsq);
# ##        ff = alphasq * c2 * (1.+c2 * b) * db
#         
# ##        sc + = phaseTools.add_over_faces_inline(self.ff,-self.dphi[:,1],self.dphi[:,0],mesh)

	self.value = sc.getValue()
