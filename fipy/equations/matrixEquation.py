#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "matrixEquation.py"
 #                                    created: 11/12/03 {10:41:06 AM} 
 #                                last update: 2/20/04 {5:08:35 PM} 
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
import spmatrix

from fivol.equations.equation import Equation
import fivol.tools.vector as vector
from fivol.tools.dimensions.physicalField import PhysicalField

class MatrixEquation(Equation):
    bandwidth = 5
    
    def getVar(self):
        return self.var    
	
    def getL(self):
	return self.L
	
    def getB(self):
	return self.b

    def buildMatrix(self):
	N = len(self.array)
	self.L = spmatrix.ll_mat(N,N,self.bandwidth)
	self.b = Numeric.zeros((N),'d')
##	coeffScale = self.terms[0].getCoeffScale()
##        print self.__class__.__name__,coeffScale
##        print
	varScale = PhysicalField(1, self.var.getUnit())
	for term in self.terms:
	    term.buildMatrix(self.L,self.var.getOld().getValue(),self.b,self.terms[0].getCoeffScale(),varScale)
	    
    def postSolve(self, array):
	pass
	
    def getResidual(self):
	Lx = self.oldSweepArray.copy()
	self.L.matvec(self.oldSweepArray,Lx)
	
	residual = Lx - self.b
	
	denom = Numeric.where(Lx == 0, self.b, Lx)
	denom = Numeric.where(denom == 0, 1, denom)
	
	residual /= denom
		
	return abs(residual) + self.solutionTolerance * 1e-10
	    
    def getResidual2(self):
	Lx = self.oldSweepArray.copy()
	self.L.matvec(self.oldSweepArray,Lx)
	
	# prevent divide-by-zero
## 	epsilon = 1e-60
## 	b = Numeric.where(self.b == 0, epsilon, self.b)
## 	residual = Numeric.where(residual == 0, b, residual)
## 	b = Numeric.where(abs(self.b) < epsilon, epsilon, self.b)
## 	residual = Numeric.where(abs(residual) < epsilon, b, residual)
## 	residual = Numeric.where(residual == 0, epsilon, residual)
## 	b = Numeric.where(self.b == 0, epsilon, self.b)
## 	print self
## 	print "L.x:", residual 
## 	print "b:", b
## 	print "b/L.x:", b / residual
## 	residual = 1. - b / residual
## 	residual = Numeric.where(residual == 0, residual - self.b, 1. - self.b / residual)
## 	residual = Numeric.where(residual != epsilon, 1. - b / residual, -b)
	
	residual = Lx - self.b
	residual = vector.sqrtDot(residual,residual)
	
## 	denom = Numeric.where(Lx == 0, self.b, Lx)
## 	denom = Numeric.where(denom == 0, 1., denom)
## 	
## 	residual /= denom
## 	
## 	residual = vector.sqrtDot(residual,residual)
	
	Lx = vector.sqrtDot(Lx,Lx)
	if Lx != 0:
	    residual /= Lx
	else:
	    b = vector.sqrtDot(self.b,self.b)
	    if b != 0:
		residual /= b
	    else:
		residual = 0
		
	return residual

    def solve(self):
	self.array = self.var.getNumericValue()
	self.oldSweepArray = self.array.copy()
	self.buildMatrix()
	self.solver.solve(self.L,self.array,self.b)
	residual = self.getResidual()
	self.postSolve(residual)
	self.var[:] = self.array
	
	self.residual = residual
	self.converged = Numeric.alltrue(self.residual < self.solutionTolerance)
	
