#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "constraintMask.py"
 #
 # Author: Jonathan Guyer <guyer@nist.gov>
 # Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 # Author: James Warren   <jwarren@nist.gov>
 #   mail: NIST
 #    www: <http://www.ctcms.nist.gov/fipy/>
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
 # #############################################################################
 ##

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools import numerix

def _ConstraintMask(var):
    class _ConstraintMaskVariable(var._variableClass):
        def __init__(self, var):
            super(_ConstraintMaskVariable, self).__init__(mesh=var.mesh, rank=0, value=False)
            for constraint in var.constraints:
                self._requires(constraint.where)
            self.var = var

        def _calcValue(self):
            returnMask = numerix.zeros(self.shape, dtype=bool)    
            for constraint in self.var.constraints:
                returnMask = returnMask | numerix.asarray(constraint.where)
            return returnMask

    return _ConstraintMaskVariable(var)
