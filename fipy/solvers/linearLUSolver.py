"""
-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "linearLUSolver.py"
                                   created: 11/14/03 {3:56:49 PM} 
                               last update: 12/4/03 {10:23:29 PM} 
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
import superlu
import sys

class LinearLUSolver(Solver):
    def __init__(self):
	Solver.__init__(self, tolerance = 0., steps = 0)
	
    def solve(self, L, x, b):
        print L
        print b
        
        LU = superlu.factorize(L.to_csr(), diag_pivot_thresh = 0.)
        LU.solve(b, x)
        

