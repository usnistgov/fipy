#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "linearGMRESSolver.py"
 #                                    created: 11/14/03 {3:56:49 PM} 
 #                                last update: 12/6/04 {4:32:02 PM} 
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

__docformat__ = 'restructuredtext'

import sys

import precon
import itsolvers

from fipy.solvers.solver import Solver

class LinearGMRESSolver(Solver):
    """
    
    The `LinearGMRESSolver` solves a linear system of equations using the
    Generalised Minimal Residual method (GMRES) with jacobi
    preconditioning. GMRES solves systems with a general non-symmetric
    coefficient matrix.

    The `LinearGMRESSolver` is a wrapper class for the the PySparse_
    `itsolvers.gmres` and `precon.jacobi` methods. Usage:

    ::

        solver = LinearGMRESSolver(tolerance = 1e-10, steps = 1000)

    .. _PySparse: http://pysparse.sourceforge.net
    
    """
    
    def _solve(self, L, x, b):

## 	print "L: ", L
## 	print "b: ", b
## 	print "x: ", x
	
	A = L._getMatrix().to_csr()
        
        Assor=precon.jacobi(L._getMatrix())
        
        info, iter, relres = itsolvers.gmres(A,b,x,self.tolerance,self.steps,Assor)
        
## 	print info, iter, relres
	
## 	y = x.copy()
## 	L.matvec(x,y)
## 	print "L * x: ", y
	
	if (info != 0):
	    print >> sys.stderr, 'cg not converged'
