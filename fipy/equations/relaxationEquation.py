#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "relaxationEquation.py"
 #                                    created: 11/12/03 {10:41:06 AM} 
 #                                last update: 9/3/04 {10:35:34 PM} 
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
 # protection and is in the public domain.  FiPy is an experimental
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

from fipy.equations.matrixEquation import MatrixEquation

class RelaxationEquation(MatrixEquation):
    def __init__(
	self,
	var,
	terms,
	solver,
	solutionTolerance = 1e-4,
	relaxation = 1.):

	self.relaxation = relaxation
	MatrixEquation.__init__(self, var, terms, solver, solutionTolerance)
	
    def postSolve(self, residual):
## 	self.relaxation *= (self.solutionTolerance/abs(residual))**0.2
	self.relaxation = (self.solutionTolerance/abs(residual))**0.2
## 	self.relaxation = self.solutionTolerance/abs(residual)
	self.relaxation = Numeric.where(self.relaxation > 1., 1., self.relaxation)
	self.relaxation = Numeric.where(self.relaxation < 1e-12, 1e-12, self.relaxation)
## 	self.relaxation = min(1.5,self.relaxation)
## 	self.relaxation = max(1e-3,self.relaxation)
	
    def getRelaxation(self):
	return self.relaxation