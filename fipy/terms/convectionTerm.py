#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "convectionTerm.py"
 #                                    created: 11/13/03 {11:39:03 AM} 
 #                                last update: 2/25/05 {5:35:41 PM} 
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
 #  2003-11-13 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import Numeric

from fipy.terms.faceTerm import FaceTerm
from fipy.variables.vectorFaceVariable import VectorFaceVariable

class ConvectionTerm(FaceTerm):
    def __init__(self, coeff = 1.0, diffusionTerm = None):
        """
        `Convection Term` should not be called directly.
        
        :Parameters:
          - `coeff` : The `term`'s coefficient value.
          - `diffusionTerm` : If a `DiffusionTerm` is given, the `ConvectionTerm` uses the diffusion coefficient to calculate the Peclet number.
        """
	self.diffusionTerm = diffusionTerm
        self.stencil = None
	FaceTerm.__init__(self, coeff = coeff)
	
    def __neg__(self):
        """
        Negate the term.

           >>> -ConvectionTerm(coeff = 1.0)
           ConvectionTerm(coeff = -1.0)
        """
        return self.__class__(coeff = -self.coeff, diffusionTerm = self.diffusionTerm)

    def _calcGeomCoeff(self, mesh):
	if not isinstance(self.coeff, VectorFaceVariable):
	    self.coeff = VectorFaceVariable(mesh = mesh, value = self.coeff)

	projectedCoefficients = self.coeff * mesh._getOrientedAreaProjections()
	
	self.geomCoeff = projectedCoefficients.sum(1)
	
    def _getWeight(self, mesh):

        if self.stencil is None:


            if self.diffusionTerm == None:
                diffCoeff = 1e-20
            else:
                diffCoeff = self.diffusionTerm._getGeomCoeff(mesh)
                if diffCoeff == 0.:
                    diffCoeff = 1e-20

            alpha = self.Alpha(-self._getGeomCoeff(mesh) / diffCoeff)
            
            self.stencil = {'implicit' : {'cell 1 diag'    : alpha,
                                          'cell 1 offdiag' : (1-alpha),
                                          'cell 2 diag'    : -(1-alpha),
                                          'cell 2 offdiag' : -alpha}}

	return self.stencil

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
