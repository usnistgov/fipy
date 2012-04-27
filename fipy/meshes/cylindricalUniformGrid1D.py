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

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.meshes.uniformGrid1D import UniformGrid1D
from fipy.tools import numerix
from fipy.tools import parallel

__all__ = ["CylindricalUniformGrid1D"]

class CylindricalUniformGrid1D(UniformGrid1D):
    """
    Creates a 1D cylindrical grid mesh.
    
        >>> mesh = CylindricalUniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
    """
    def __init__(self, dx=1., nx=1, origin=(0,), overlap=2, communicator=parallel, *args, **kwargs):
        UniformGrid1D.__init__(self, dx=dx, nx=nx, origin=origin, overlap=overlap, communicator=communicator, *args, **kwargs)
    
    def _translate(self, vector):
        return CylindricalUniformGrid1D(dx=self.args['dx'],
                                        nx=self.args['nx'],
                                        origin=self.args['origin'] + numerix.array(vector),
                                        overlap=self.args['overlap'])
                 
    @property
    def _faceAreas(self):
        return self.faceCenters[0]

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
        return self._faceNormals * self._faceAreas
 
    @property
    def cellVolumes(self):
        return self.dx * self.cellCenters[0]

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. The following was broken, now fixed.

            >>> from fipy import *
            >>> mesh = CylindricalUniformGrid1D(nx=3., dx=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
