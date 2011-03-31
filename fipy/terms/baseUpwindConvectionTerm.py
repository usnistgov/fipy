#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "baseUpwindConvectionTerm.py"
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

class _BaseUpwindConvectionTerm(ConvectionTerm):

    class _Alpha(FaceVariable):
        def __init__(self, P):
            FaceVariable.__init__(self, mesh = P.mesh)
            self.P = self._requires(P)
            
        def _calcValuePy(self, P):
            alpha = numerix.where(P > 0., 1., 0.)
            return PhysicalField(value = alpha)

        def _calcValueIn(self, P):
            alpha = self._array.copy()
            inline._runInline("""
                alpha[i] = 0.5;
                
                if (P[i] > 0.) {
                    alpha[i] = 1.;
                } else {
                    alpha[i] = 0.;
                }
            """,
            alpha = alpha, P = P,
            ni = self.mesh._numberOfFaces
            )

            return self._makeValue(value = alpha)

        def _calcValue(self):
            P  = self.P.numericValue

            return inline._optionalInline(self._calcValueIn, self._calcValuePy, P)


def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
