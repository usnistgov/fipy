"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "linearPCGSolver.py"
                                   created: 11/14/03 {3:56:49 PM} 
                               last update: 11/24/03 {7:12:04 PM} 
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
 2003-11-14 JEG 1.0 original
###################################################################
"""

from solver import Solver
import precon
import itsolvers
import sys

class LinearPCGSolver(Solver):
    def __init__(self, tolerance, steps):
	Solver.__init__(self, tolerance, steps)
	
    def solve(self, L, x, b):

	A = L.to_sss()
	
	Assor=precon.ssor(A)

	info, iter, relres = itsolvers.pcg(A,b,x,self.tolerance,self.steps,Assor)
        print info, iter, relres
	    
	if (info != 0):
	    print >> sys.stderr, 'cg not converged'
