"""
#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "transientTerm.py"
#                                   created: 11/12/03 {11:35:45 AM} 
#                               last update: 11/14/03 {5:04:34 PM} 
# Author: Jonathan Guyer
# E-mail: guyer@nist.gov
#   mail: NIST
#    www: http://ctcms.nist.gov/
# 
#========================================================================
#This software was developed at the National Institute of Standards
#and Technology by employees of the Federal Government in the course
#of their official duties.  Pursuant to title 17 Section 105 of the
#United States Code this software is not subject to copyright
#protection and is in the public domain.  PFM is an experimentalœf
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

import cellTerm

class TransientTerm(cellTerm.CellTerm):
    def __init__(self,equation,tranCoeff):
	cellTerm.CellTerm.__init__(self, stencil = (0,1,1), equation) 
	self.tranCoeff = tranCoeff
	    
    def updateCoeff(self,dt):
	cells = self.equation.mesh().cells()
	self.coeff = Numeric.zeroes([len(cells),'d')
	for cell in cells:
	    self.coeff[cell.id()] = self.tranCoeff * cell.volume() / dt
	

