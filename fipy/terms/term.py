#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "term.py"
 #                                    created: 11/12/03 {10:54:37 AM} 
 #                                last update: 12/7/04 {12:05:54 PM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
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

from fipy.variables.variable import Variable
from fipy.tools.dimensions.physicalField import PhysicalField

class Term:
    def __init__(self):
	self.coeff = None
        
    def buildMatrix(self, var, boundaryConditions, oldArray, dt):
	pass
	
    def getFigureOfMerit(self):
	return None

    def getResidual(self, matrix, var, RHSvector):
	Lx = matrix * Numeric.array(var[:])
	
	residual = Lx - RHSvector
	
	denom = max(abs(Lx))
	if denom == 0:
	    denom = max(abs(RHSvector))
	if denom == 0:
	    denom = 1.
	    
	residual /= denom
		
	return abs(residual)

    def isConverged(self):
	return self.converged

    def solve(self, var, solver = None, boundaryConditions = (), dt = 1., solutionTolerance = 1e-4):
 	matrix, RHSvector = self.buildMatrix(var, boundaryConditions, dt = dt)
        residual = self.getResidual(matrix, var, RHSvector)
	if solver is None:
	    from fipy.solvers.linearPCGSolver import LinearPCGSolver
	    solver = LinearPCGSolver()
	    
	array = var.getNumericValue()
	solver.solve(matrix, array, RHSvector)
	var[:] = array
	
	self.residual = residual
	self.converged = Numeric.alltrue(self.residual < solutionTolerance)

    def __add__(self, other):
	from fipy.terms.binaryTerm import AdditionTerm
	return AdditionTerm(term1 = self, term2 = other)
	    
    __radd__ = __add__
    
    def __sub__(self, other):
	from fipy.terms.binaryTerm import SubtractionTerm
	return SubtractionTerm(term1 = self, term2 = other)

    def __rsub__(self, other):
	from fipy.terms.binaryTerm import SubtractionTerm
	return SubtractionTerm(term1 = other, term2 = self)

    def calcCoeff(self, mesh):
	pass
	
    def getCoeff(self, mesh):
	if self.coeff is None:
	    self.calcCoeff(mesh)
	return self.coeff
	
    def getWeight(self, mesh):
	pass
	    
