"""
## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "exponentialConvectionTerm.py"
 #                                    created: 12/5/03 {2:50:05 PM} 
 #                                last update: 12/22/03 {5:02:09 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##
"""

from convectionTerm import ConvectionTerm
from variables.faceVariable import FaceVariable
import Numeric

class PowerLawConvectionTerm(ConvectionTerm):
    class Alpha(FaceVariable):
	def __init__(self, P):
	    FaceVariable.__init__(self, P.getMesh())
	    self.P = self.requires(P)
	    
	def calcValue(self):
	    eps = 1e-3
	    P  = self.P[:]
	    P = Numeric.where(Numeric.absolute(P) < eps, eps, P)
	    
	    alpha = Numeric.where(                                   P > 10.,                 (P - 1.) / P,   0.5)

	    tmp = (1. - P/10.)
	    tmpSqr = tmp * tmp
	    alpha = Numeric.where(    Numeric.logical_and(10. >= P, P > eps), ((P-1.) + tmpSqr*tmpSqr*tmp)/P, alpha)

	    tmp = (1. + P/10.)
	    tmpSqr = tmp * tmp
	    alpha = Numeric.where( Numeric.logical_and(eps  >  P, P >= -10.),     (tmpSqr*tmpSqr*tmp - 1.)/P, alpha)
	    
	    self.value = alpha
