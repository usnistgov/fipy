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

from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.variables.cellVariable import CellVariable

__all__ = ["ResidualTerm"]

class ResidualTerm(_ExplicitSourceTerm):
    r"""

    The `ResidualTerm` is a special form of explicit `SourceTerm` that adds the
    residual of one equation to another equation. Useful for Newton's method.
    """
    def __init__(self, equation, underRelaxation=1.):
        self.equation = equation
        self.underRelaxation = underRelaxation
        
        _ExplicitSourceTerm.__init__(self, var=None)
        
    def __repr__(self):
        return r"$\Delta$[" + repr(self.equation) + "]"

    def _getGeomCoeff(self, var):
        return self.coeff
        
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        vec = self.equation.justResidualVector(var=None, 
                                               boundaryConditions=boundaryConditions,
                                               dt=dt)

        self.coeff = CellVariable(mesh=var.mesh, value=vec * self.underRelaxation)
        self.geomCoeff = None
        self.coeffVectors = None
        
        return _ExplicitSourceTerm._buildMatrix(self, var=var, SparseMatrix=SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)


