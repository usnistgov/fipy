#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "cylindricalUniformGrid1D.py"
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

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.uniformGrid1D import UniformGrid1D
from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools import parallelComm

__all__ = ["CylindricalUniformGrid1D"]

class CylindricalUniformGrid1D(UniformGrid1D):
    """
    Creates a 1D cylindrical grid mesh.

        >>> mesh = CylindricalUniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]

    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2, communicator=parallelComm, *args, **kwargs):
        UniformGrid1D.__init__(self, dx=dx, nx=nx, origin=origin, overlap=overlap, communicator=communicator, *args, **kwargs)

    def _translate(self, vector):
        return CylindricalUniformGrid1D(dx=self.args['dx'],
                                        nx=self.args['nx'],
                                        origin=self.args['origin'] + numerix.array(vector),
                                        overlap=self.args['overlap'])

    @property
    def _faceAreas(self):
        return self.faceCenters[0].value

    @property
    def _cellAreas(self):
        return numerix.array((self.faceAreas[:-1], self.faceAreas[1:]))

    @property
    def _cellAreaProjections(self):
        return MA.array(self.cellNormals) * self.cellAreas

    @property
    def _faceAspectRatios(self):
        return self._faceAreas / self._cellDistances

    @property
    def _areaProjections(self):
        return self.faceNormals * self._faceAreas

    @property
    def cellVolumes(self):
        return self.dx * self.cellCenters[0].value

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. The following was broken, now fixed.

            >>> from fipy import *
            >>> mesh = CylindricalUniformGrid1D(nx=3., dx=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        This test is for https://github.com/usnistgov/fipy/issues/372. Cell
        volumes were being returned as binOps rather than arrays.

            >>> m = CylindricalUniformGrid1D(dx=1., nx=4)
            >>> print isinstance(m.cellVolumes, numerix.ndarray)
            True
            >>> print isinstance(m._faceAreas, numerix.ndarray)
            True

        If the above types aren't correct, the divergence operator's value can be a binOp

            >>> print isinstance(CellVariable(mesh=m).arithmeticFaceValue.divergence.value, numerix.ndarray)
            True

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
