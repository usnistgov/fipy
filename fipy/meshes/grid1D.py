#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid1D.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
from fipy.tools.decorators import getsetDeprecated
from mesh1D import Mesh1D
from fipy.tools import parallel

from fipy.meshes.builders import NonuniformGrid1DBuilder
from fipy.meshes.abstractGrid import AbstractGrid1DFactory

class Grid1D(AbstractGrid1DFactory(Mesh1D)):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = Grid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
        >>> mesh = Grid1D(dx = (1, 2, 3))
        >>> print mesh.cellCenters
        [[ 0.5  2.   4.5]]
         
        >>> mesh = Grid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

    """
    def __init__(self, dx=1., nx=None, overlap=2, 
                 communicator=parallel,
                 BuilderClass=NonuniformGrid1DBuilder):

        builder = BuilderClass()

        self.args = {
            'dx': dx, 
            'nx': nx, 
            'overlap': overlap
        }

        builder.buildGridData([dx], [nx], overlap, communicator)

        ([self.dx],
         [self.nx],
         self.dim,
         scale,
         self.globalNumberOfCells,
         self.globalNumberOfFaces,
         self.overlap,
         self.offset,
         self.numberOfVertices,
         self.numberOfFaces,
         self.numberOfCells,
         self.shape,
         self.physicalShape,
         self._meshSpacing,
         self.occupiedNodes,
         vertices,
         faces,
         cells) = builder.gridData

        Mesh1D.__init__(self, vertices, faces, cells, communicator=communicator)
        
        self.scale = scale

    
## pickling

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when nx is a float.

            >>> # from fipy import *
            >>> from fipy import CellVariable, DiffusionTerm
            >>> mesh = Grid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
