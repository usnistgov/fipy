#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearPCGSolver.py"
 #                                    created: 11/14/03 {3:56:49 PM} 
 #                                last update: 9/3/04 {10:35:32 PM} 
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
 #  2003-11-14 JEG 1.0 original
 # ###################################################################
 ##

import sys

import precon
import itsolvers

from fipy.solvers.solver import Solver

class LinearPCGSolver(Solver):
    def __init__(self, tolerance, steps):
	Solver.__init__(self, tolerance, steps)
	
    def solve(self, L, x, b):
## 	print 'L:',L
## 	print 'x:',x
## 	print 'b:',b
## 	raw_input('end output')
    
	A = L.getMatrix().to_sss()

	Assor=precon.ssor(A)

 	info, iter, relres = itsolvers.pcg(A,b,x,self.tolerance,self.steps,Assor)
##        print info, iter, relres

	if (info != 0):
	    print >> sys.stderr, 'cg not converged'
