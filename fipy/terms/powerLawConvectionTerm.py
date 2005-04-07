#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "powerLawConvectionTerm.py"
 #                                    created: 12/5/03 {2:50:05 PM} 
 #                                last update: 4/1/05 {11:03:04 AM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
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

__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.convectionTerm import _ConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.tools.inline import inline

class PowerLawConvectionTerm(_ConvectionTerm):
    r"""
    The discretization for the `PowerLawConvectionTerm` is given by

    .. raw:: latex
    
       $$ \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f $$

       where $ \phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A $ and
       $\alpha_f$ is calculated using the power law scheme.
       For further details see Section 3.5 of the main \FiPy{}
       guide~\cite{FiPyGuide}.

    """    
    class Alpha(FaceVariable):
	def __init__(self, P):
	    FaceVariable.__init__(self, mesh = P.getMesh())
	    self.P = self._requires(P)
	    
	def _calcValuePy(self, eps, P):
	    P = Numeric.where(abs(P) < eps, eps, P)
	    
## 	    print "P:", P
	    
	    alpha = Numeric.where(                    P > 10.,                     (P - 1.) / P,   0.5)

	    tmp = (1. - P / 10.)
	    tmpSqr = tmp * tmp
	    alpha = Numeric.where(   (10. >= P) and (P > eps), ((P-1.) + tmpSqr*tmpSqr*tmp) / P, alpha)

	    tmp = (1. + P / 10.)
	    tmpSqr = tmp * tmp
	    alpha = Numeric.where((eps  >  P) and (P >= -10.),     (tmpSqr*tmpSqr*tmp - 1.) / P, alpha)
	    
	    alpha = Numeric.where(                   P < -10.,                          -1. / P, alpha)
	    
## 	    print "alpha:", alpha
## 	    raw_input()
	    
	    self.value = PhysicalField(value = alpha)

	def _calcValueIn(self, eps, P):
##            print P.shape
##            print len(self.mesh.getCells())
##            raw_input()
	    inline._runInlineLoop1("""
		if (fabs(P(i)) < eps) {
		    P(i) = eps;
		}
		
		alpha(i) = 0.5;
		
		if (P(i) > 10.) {
		    alpha(i) = (P(i) - 1.) / P(i);
		} else if (10. >= P(i) && P(i) > eps) {
		    double	tmp = (1. - P(i) / 10.);
		    double	tmpSqr = tmp * tmp;
		    alpha(i) = ((P(i) - 1.) + tmpSqr*tmpSqr*tmp) / P(i);
		} else if (eps >= P(i) && P(i) >= -10) {
		    double	tmp = (1. + P(i) / 10.);
		    double	tmpSqr = tmp * tmp;
		    alpha(i) = (tmpSqr*tmpSqr*tmp - 1.) / P(i);
		} else {	// P(i) < -10.
		    alpha(i) = -1. / P(i);
		}
	    """,
	    alpha = self._getArray(), eps = eps, P = P,
	    ni = len(self.mesh.getFaces())
	    )


	def _calcValue(self):	    
	    eps = 1e-3
	    P  = self.P.getNumericValue()
	    
	    inline._optionalInline(self._calcValueIn, self._calcValuePy, eps, P)
