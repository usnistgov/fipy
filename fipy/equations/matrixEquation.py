"""
#-----*-Pyth-*-
####################################################################
# PFM - Python-based phase field solver
#
# FILE: "matrixEquation.py"
#                                   created: 11/12/03 {10:41:06 AM} 
#                               last update: 11/14/03 {4:52:24 PM} 
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

import equation
import Numeric
import spmatrix


class MatrixEquation(equation.Equation):
    bandwidth = 5
    
    def __init__(self,name,mesh,terms,solver):
	self.mesh = mesh
	equation.Equation.__init__(
	    self,
	    name,
	    var = Numeric.zeroes([len(mesh.cells())],'d'),
	    terms,
	    solver)
	
    def L(self):
	return self.L
	
    def b(self):
	return self.b
	    
    def solve(self):
	N = var.size()
	self.L = spmatrix.ll_mat_sym(N,self.bandwidth*N)
	self.b = Numeric.zeros((N),'d')
	for term in self.terms:
	    term.buildMatrix()
	self.solver.solve(self.L,self.var,self.b)
	