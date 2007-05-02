#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "implicitSourceTerm.py"
 #                                    created: 11/28/03 {11:36:25 AM} 
 #                                last update: 3/29/07 {10:40:42 AM} 
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
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.sourceTerm import SourceTerm

class ImplicitSourceTerm(SourceTerm):
    r"""

    The `ImplicitSourceTerm` represents

    .. raw:: latex

       \[ \int_V \phi S \,dV \simeq \phi_P S_P V_P \] where $S$ is the

    `coeff` value.       
    """
    def _calcCoeffVectors(self, var, equation=None):
        coeff = self._getGeomCoeff(var.getMesh())
        from fipy.tools.numerix import sign
        combinedSign = self._diagonalSign * sign(coeff)
        self.coeffVectors = {
            'diagonal': coeff * (combinedSign >= 0),
            'old value': 0,
            'b vector': -coeff * var * (combinedSign < 0),
            'new value': 0
        }
