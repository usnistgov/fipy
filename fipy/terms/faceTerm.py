"""
#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "faceTerm.py"
#                                   created: 11/12/03 {11:01:56 AM} 
#                               last update: 11/13/03 {11:36:08 AM} 
# Author: Jonathan Guyer
# Author: Daniel Wheeler
# E-mail: guyer@nist.gov
#   mail: NIST
#    www: http://ctcms.nist.gov/
# 
#========================================================================
#This software was developed at the National Institute of Standards
#and Technology by employees of the Federal Government in the course
#of their official duties.  Pursuant to title 17 Section 105 of the
#United States Code this software is not subject to copyright
#protection and is in the public domain.  PFM is an experimental
#system.  NIST assumes no responsibility whatsoever for its use by
#other parties, and makes no guarantees, expressed or implied, about
#its quality, reliability, or any other characteristic.  We would
#appreciate acknowledgement if the software is used.
#
#This software can be redistributed and/or modified freely
#provided that any derivative works bear some notice that they are
#derived from it, and any modified versions bear some notice that
#they have been modified.
#========================================================================
# See the file "license.terms" for information on usage and  redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
# 
####################################################################
#----
"""

import term

class FaceTerm(term.Term):
    def __init__(self,stencil,equation):
	"""
	stencil = [phi_adj, phi]
	"""
	term.Term.__init__(self,stencil,equation)
	
    def buildMatrix(self):
	var = self.equation.var()
	N = var.size()
	
	for face in var.mesh().faces():
	    if len(face.cells()) == 2:
		id1 = face.cells()[0].id()
		id2 = face.cells()[1].id()
		self.equation.L()[id1,id1]+=self.coeff * self.stencil[1]
		self.equation.L()[id1,id2]-=self.coeff * self.stencil[0]
		self.equation.L()[id2,id1]-=self.coeff * self.stencil[0]
		self.equation.L()[id2,id2]+=self.coeff * self.stencil[1]
	    else if len(face.cells()) == 1:
# 		do boundary conditions
# 		for the moment this is zero flux
		pass
	    else:
# 		cause catastrophic failure
		pass
		    
		    