"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "matrixEquation.py"
                                   created: 11/12/03 {10:41:06 AM} 
                               last update: 11/17/03 {4:13:57 PM} 
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
 2003-11-12 JEG 1.0 original
###################################################################
"""

import equation
import Numeric
import spmatrix


class MatrixEquation(equation.Equation):
    bandwidth = 5
    
    def __init__(self,
        name,
        var,
        terms,
        solver):

	equation.Equation.__init__(self,
	    name,
	    var,
	    terms,
	    solver)

    def getVar(self):
        return self.var    
	
    def getL(self):
	return self.L
	
    def getB(self):
	return self.b

    def solve(self,dt):
        array=self.var.getArray()
	N = len(array)
	self.L = spmatrix.ll_mat(N,N,self.bandwidth)
	self.b = Numeric.zeros((N),'d')
	for term in self.terms:
	    term.updateCoeff(dt)
	    term.buildMatrix(self.L,array,self.b)
	self.solver.solve(self.L,array,self.b)
	
