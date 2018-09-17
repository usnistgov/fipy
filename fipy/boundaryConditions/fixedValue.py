#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "fixedValue.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James Warren <jwarren@nist.gov>
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

"""Fixed value (Dirichlet) boundary condition
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.tools import numerix
from fipy.tools import vector
from fipy.variables.variable import Variable

__all__ = ["FixedValue"]

class FixedValue(BoundaryCondition):
    r"""
    The `FixedValue` boundary condition adds a contribution, equivalent to a
    fixed value (Dirichlet condition), to the equation's RHS vector and
    coefficient matrix. The contributions are given by :math:`-\mathtt{value}
    \times G_{\text{face}}` for the RHS vector and :math:`G_{\text{face}}` for
    the coefficient matrix. The parameter :math:`G_{\text{face}}` represents the
    term's geometric coefficient, which depends on the type of term and the mesh
    geometry.

    Contributions are only added to entries corresponding to the
    specified faces.
    """


    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Set boundary equal to value.

        A `tuple` of (`LL`, `bb`) is calculated, to be added to the
        Term's (:math:`\mathsf{L}`, :math:`\mathsf{b}`) matrices.

        :Parameters:
          - `SparseMatrix`: Sparse matrix class to use
          - `Ncells`:       Size of matrices
          - `MaxFaces`:     bandwidth of :math:`\mathsf{L}`
          - `coeff`:        contribution to adjacent cell diagonal and
            :math:`\mathsf{b}`-vector by this exterior face
        """
        faces = self.faces.value

        LL = SparseMatrix(mesh=self.faces.mesh, sizeHint=len(self.faces), bandwidth=1)
        LL.addAt(coeff['cell 1 diag'][faces], self.adjacentCellIDs, self.adjacentCellIDs)

        ## The following has been commented out because
        ## FixedValue's _buildMatrix() method is called for
        ## each term in the equation. Thus minusCoeff can be different for each term.
        ##
        ## if not hasattr(self, 'minusCoeff'):
        ##     self.minusCoeff = -coeff['cell 1 offdiag']
        ##     self.minusCoeff.dontCacheMe()

        bb = numerix.zeros((Ncells,),'d')

        value = self.value
        if isinstance(value, Variable):
            value = value.value
        if value.shape == faces.shape:
            value = value[faces]

        vector.putAdd(bb, self.adjacentCellIDs, -coeff['cell 1 offdiag'].value[faces] * value)

        return (LL, bb)
