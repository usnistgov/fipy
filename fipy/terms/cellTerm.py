"""-*-Pyth-*-
###################################################################
 PFM - Python-based phase field solver

 FILE: "cellTerm.py"
                                   created: 11/12/03 {11:00:54 AM} 
                               last update: 11/17/03 {11:33:32 AM} 
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

import Numeric
import term

class CellTerm(term.Term):
    def __init__(self,stencil,equation):
	"""
	stencil = [b, phi, phi_old]
	"""
	term.Term.__init__(self,stencil,equation)
	
    def buildMatrix(self):
	var = self.equation.getVar()
	N = len(var)
        b = self.equation.getB()
        L = self.equation.getL()
        stencil = self.stencil
	b += var*self.coeff*self.stencil[2]
	b += Numeric.ones([N])*self.coeff*self.stencil[0]
	L.update_add_pyarray(Numeric.ones([N])*self.coeff*stencil[1])
	
		
		
		
