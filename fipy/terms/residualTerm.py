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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
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
