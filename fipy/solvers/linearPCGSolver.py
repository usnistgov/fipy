"""-*-Pyth-*-
###################################################################
 Alpha - Core code development for Alpha

 FILE: "linearPCGSolver.py"
                                   created: 11/14/03 {3:56:49 PM} 
                               last update: 11/14/03 {4:08:41 PM} 
 Author: Jonathan Guyer
 E-mail: jguyer@his.com
   mail: Alpha Cabal
    www: http://alphatcl.sourceforge.net
 
Copyright (c) 2003  Jonathan Guyer

See the file "license.terms" for information on usage and redistribution
of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 See the file "license.terms" for information on usage and  redistribution
 of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 
###################################################################
"""

import solver

class linearPCGSolver(solver.Solver):
    def __init__(self, tolerance, steps):
	solver.Solver.__init__(self, tolerance, steps)
	
    def solve(self, L, x, b):
	A = self.L.to_sss()
	
	Assor=precon.ssor(A)
	
	info, iter, relres = itsolvers.pcg(A,b,x,self.tolerance,self.steps,Assor)
	
# 	print info, iter, relres
	    
	if (info != 0):
	    print >> sys.stderr, 'cg not converged'
	
