#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "powerLawConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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

from fipy.terms.convectionTerm import ConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import numerix

class PowerLawConvectionTerm(ConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::
    
       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the power law scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """    
    class _Alpha(FaceVariable):
	def __init__(self, P):
	    FaceVariable.__init__(self, mesh = P.getMesh())
	    self.P = self._requires(P)
	    
	def _calcValuePy(self, eps, P):
            """

            Test case added because `and` was being used instead of bitwise `&`.

                >>> from fipy.meshes.grid1D import Grid1D
                >>> mesh = Grid1D(nx = 3)
                >>> from fipy.variables.faceVariable import FaceVariable
                >>> P = FaceVariable(mesh = mesh, value = (1e-3, 1e+71, 1e-3, 1e-3))
                >>> alpha = PowerLawConvectionTerm._Alpha(P)
                >>> print numerix.allclose(alpha, [ 0.5,  1.,   0.5 , 0.5])
                True
                
            """
            
	    P = numerix.where(abs(P) < eps, eps, P)
	    
	    alpha = numerix.where(                  P > 10.,                     (P - 1.) / P,   0.5)

	    tmp = (1. - P / 10.)
	    tmpSqr = tmp * tmp
	    alpha = numerix.where(  (10. >= P) & (P > eps), ((P-1.) + tmpSqr*tmpSqr*tmp) / P, alpha)

	    tmp = (1. + P / 10.)
	    tmpSqr = tmp * tmp
	    alpha = numerix.where((-eps >  P) & (P >= -10.),     (tmpSqr*tmpSqr*tmp - 1.) / P, alpha)

	    alpha = numerix.where(                 P < -10.,                          -1. / P, alpha)

	    return PhysicalField(value = alpha)

	def _calcValueIn(self, eps, P):
            alpha = self._getArray().copy()
            
	    inline._runInline("""
		if (fabs(P[i]) < eps) {
		    P[i] = eps;
		}
		
		alpha[i] = 0.5;
		
		if (P[i] > 10.) {
		    alpha[i] = (P[i] - 1.) / P[i];
		} else if (10. >= P[i] && P[i] > eps) {
		    double	tmp = (1. - P[i] / 10.);
		    double	tmpSqr = tmp * tmp;
		    alpha[i] = ((P[i] - 1.) + tmpSqr*tmpSqr*tmp) / P[i];
		} else if (-eps > P[i] && P[i] >= -10.) {
		    double	tmp = (1. + P[i] / 10.);
		    double	tmpSqr = tmp * tmp;
		    alpha[i] = (tmpSqr*tmpSqr*tmp - 1.) / P[i];
		} else if (P[i] < -10.) {
		    alpha[i] = -1. / P[i];
		}
	    """,
	    alpha = alpha, eps = eps, P = P,
	    ni = self.mesh._getNumberOfFaces()
	    )

            return self._makeValue(value = alpha)
##         return self._makeValue(value = alpha, unit = self.getUnit())


	def _calcValue(self):	    
	    eps = 1e-3
	    P  = self.P.getNumericValue()
	    
            return inline._optionalInline(self._calcValueIn, self._calcValuePy, eps, P)

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
