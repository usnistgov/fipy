#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "explicitSourceTerm.py"
 #
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.terms.sourceTerm import SourceTerm

class _ExplicitSourceTerm(SourceTerm):
    r"""

    The `_ExplicitSourceTerm` discretisation is given by

    .. math::

       \int_V S \,dV \simeq S_P V_P 
       
    where :math:`S` is the `coeff` value. This source is added to the RHS vector and
    does not contribute to the solution matrix.

    """
        
    def _getWeight(self, var, transientGeomCoeff=None, diffusionGeomCoeff=None):
        return {
            'b vector': -1, 
            'new value': 0, 
            'old value': 0, 
            'diagonal' : 0
        }
        
    def __repr__(self):
        return repr(self.coeff)

