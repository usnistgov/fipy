#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cylindricalGrid1D.py"
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

from fipy.tools import numerix
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.meshes.numMesh.grid1D import Grid1D

class CylindricalGrid1D(Grid1D):
    """
    Creates a 1D cylindrical grid mesh.
    
        >>> mesh = CylindricalGrid1D(nx = 3)
        >>> print mesh.getCellCenters()
        [[ 0.5  1.5  2.5]]
         
        >>> mesh = CylindricalGrid1D(dx = (1, 2, 3))
        >>> print mesh.getCellCenters()
        [[ 0.5  2.   4.5]]
         
        >>> mesh = CylindricalGrid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

        >>> mesh = CylindricalGrid1D(nx=2, dx=(1., 2.)) + ((1.,),)
        >>> print mesh.getCellCenters()
        [[ 1.5  3. ]]
        >>> from fipy.tools import parallel
        >>> print parallel.procID > 0  or numerix.allclose(mesh.getCellVolumes(), (1.5, 6))
        True
        
    """
    def __init__(self, dx=1., nx=None, origin=(0,), overlap=2):
        scale = PhysicalField(value=1, unit=PhysicalField(value=dx).getUnit())
        self.origin = PhysicalField(value=origin)
        self.origin /= scale
    
        Grid1D.__init__(self, dx=dx, nx=nx, overlap=overlap)

        self.args['origin'] = origin


    def _calcFaceCenters(self):
        Grid1D._calcFaceCenters(self)
        self.faceCenters += self.origin
        
    def _calcFaceAreas(self):
        self.faceAreas = self.getFaceCenters()[0]

    def _calcCellVolumes(self):
        Grid1D._calcCellVolumes(self)
        self.cellVolumes = self.cellVolumes / 2.

    def _translate(self, vector):
        return CylindricalGrid1D(dx=self.args['dx'], nx=self.args['nx'], 
                                 origin=numerix.array(self.args['origin']) + vector, overlap=self.args['overlap'])
                                     
    def __mul__(self, factor):
        return CylindricalGrid1D(dx=self.args['dx'] * factor, nx=self.args['nx'], 
                                 origin=numerix.array(self.args['origin']) * factor, overlap=self.args['overlap'])

    def getVertexCoords(self):
        return self.vertexCoords + self.origin

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when nx is a float.

            >>> from fipy import CellVariable, DiffusionTerm
            >>> mesh = CylindricalGrid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
