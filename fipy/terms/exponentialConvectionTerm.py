#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "exponentialConvectionTerm.py"
 #                                    created: 12/5/03 {2:50:05 PM} 
 #                                last update: 4/1/05 {11:02:53 AM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

import Numeric

from fipy.terms.convectionTerm import ConvectionTerm
from fipy.variables.faceVariable import FaceVariable

class ExponentialConvectionTerm(ConvectionTerm):
    class Alpha(FaceVariable):
	def __init__(self, P):
	    FaceVariable.__init__(self, P.getMesh())
	    self.P = self._requires(P)
	    
	def _calcValue(self):
	    eps = 1e-3
	    P  = self.P.getNumericValue()

	    P = Numeric.where(abs(P) < eps, eps, P)
	    alpha = Numeric.where(P > 101., (P - 1) / P, 0.5)
	    alpha = Numeric.where(abs(P) > eps and P <= 101., ((P - 1) * Numeric.exp(P) + 1) / (P * (Numeric.exp(P) - 1)), alpha)

	    self.value = alpha
