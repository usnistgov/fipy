#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "fixedFlux.py"
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

from fipy.tools import numerix

from fipy.boundaryConditions.boundaryCondition import BoundaryCondition
from fipy.boundaryConditions.fixedValue import FixedValue
from fipy.tools import vector

__all__ = ["FixedFlux"]

class FixedFlux(BoundaryCondition):
    r"""

    The `FixedFlux` boundary condition adds a contribution, equivalent to a
    fixed flux (Neumann condition), to the equation's RHS vector.  The
    contribution, given by `value`, is only added to entries corresponding to
    the specified `faces`, and is weighted by the face areas.

    """

    def __init__(self,faces,value):
        """
        Creates a `FixedFlux` object.

        :Parameters:
            - `faces`: A `list` or `tuple` of `Face` objects to which this condition applies.
            - `value`: The value to impose.

        """
        BoundaryCondition.__init__(self,faces,value)
        ## The extra index [self.faces.value] makes self.contribution the same length as self.adjacentCellIDs
        self.contribution = (self.value * self.faces.mesh._faceAreas)[self.faces.value]

    def _buildMatrix(self, SparseMatrix, Ncells, MaxFaces, coeff):
        """Leave **L** unchanged and add gradient to **b**

        :Parameters:
          - `SparseMatrix`: *unused* (empty matrix)
          - `Ncells`:       Size of **b**-vector
          - `MaxFaces`:     *unused*
          - `coeff`:        *unused*
        """

        bb = numerix.zeros((Ncells,),'d')

        if not self.boundaryConditionApplied:
            vector.putAdd(bb, self.adjacentCellIDs, -self.contribution)
            self.boundaryConditionApplied = True

        return (0, bb)

    def _getDerivative(self, order):
        if order == 1:
            return FixedValue(self.faces, self.value)
        else:
            return BoundaryCondition._getDerivative(self, order)

    def _test(self):
        """
        The following tests check that self.contributions is the same length as
        self.adjacentCellIDs.

           >>> from fipy import *
           >>> m = Grid1D(nx = 10)
           >>> v = FaceVariable(mesh=m)
           >>> bc = FixedFlux(value=v.globalValue[-1], faces=m.facesRight)
           >>> len(bc.contribution) == len(bc.adjacentCellIDs)
           True
           >>> bc = FixedFlux(value=v, faces=m.facesRight)
           >>> len(bc.contribution) == len(bc.adjacentCellIDs)
           True
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
