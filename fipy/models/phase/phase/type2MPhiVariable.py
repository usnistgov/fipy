#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "type2MPhiVariable.py"
 #                                    created: 12/24/03 {10:39:23 AM} 
 #                                last update: 7/5/05 {3:41:45 PM} 
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

from mPhiVariable import _MPhiVariable
import Numeric
import fipy.tools.array as array

class Type2MPhiVariable(_MPhiVariable):
    r"""
    The `Type2MPhiVariable` is given by,

    .. raw:: latex

        $$\phi - \frac{1}{2} - \frac{ \kappa_1 }{ \pi } \arctan \left( \kappa_2 T \right). $$

    Usage ::

        Type2MPhiVariable(phase = <CellVariable>, 
                          temperature = <CellVariable>, 
                          parameters = {'kappa1' : X, 'kappa2' : Y})

    """
    def _calcValue(self):        
        self.value = self.phase[:] - 0.5 - self.parameters['kappa 1'] / Numeric.pi * array.arctan(self.parameters['kappa 2'] * self.temperature[:])


        
