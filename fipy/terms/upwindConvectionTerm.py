#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "upwindConvectionTerm.py"
 #                                    created: 12/5/03 {2:50:05 PM} 
 #                                last update: 12/22/05 {4:00:55 PM} 
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

from fipy.terms.convectionTerm import ConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.tools.inline import inline

class UpwindConvectionTerm(ConvectionTerm):
    r"""
    The discretization for the `UpwindConvectionTerm` is given by

    .. raw:: latex
    
       $$ \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f $$

       where $ \phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A $ and
       $\alpha_f$ is calculated using the upwind convection scheme.
       For further details see ``\nameref{FiPy-sec:NumericalSchemes}'' in the
       main \FiPy{} guide\cite[\S~\ref{FiPy-sec:NumericalSchemes}]{FiPyGuide}.
    """
    
    class _Alpha(FaceVariable):
	def __init__(self, P):
	    FaceVariable.__init__(self, mesh = P.getMesh())
	    self.P = self._requires(P)
	    
	def _calcValuePy(self, P):
	    alpha = Numeric.where(P > 0., 1., 0.)

	    return PhysicalField(value = alpha)

	def _calcValueIn(self, P):
            alpha = self._getArray().copy()
            
	    inline._runInlineLoop1("""
		alpha(i) = 0.5;
		
		if (P(i) > 0.) {
		    alpha(i) = 1.;
		} else {
		    alpha(i) = 0.;
		}
	    """,
	    alpha = alpha, P = P,
	    ni = len(self.mesh.getFaces())
	    )
     
            return self._makeValue(value = alpha)
##         return self._makeValue(value = alpha, unit = self.getUnit())


	def _calcValue(self):	    
	    P  = self.P.getNumericValue()
	    
	    return inline._optionalInline(self._calcValueIn, self._calcValuePy, P)
