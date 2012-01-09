#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "sourceTerm.py"
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

from fipy.terms.cellTerm import CellTerm
from fipy.terms import AbstractBaseClassError
from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

__all__ = ["SourceTerm"]

class SourceTerm(CellTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff=0., var=None):
        if self.__class__ is SourceTerm:
            raise AbstractBaseClassError
        CellTerm.__init__(self, coeff=coeff, var=var) 
        
    def _calcGeomCoeff(self, var):
        self._checkCoeff(var)

        if self.coeff.shape != () and self.coeff.shape[-1] != len(var.mesh.cellVolumes):        
            return self.coeff[...,numerix.newaxis] * CellVariable(mesh=var.mesh, value=var.mesh.cellVolumes)
        else:
            return self.coeff * CellVariable(mesh=var.mesh, value=var.mesh.cellVolumes)

    def _checkDt(self, dt):
        return 1.
