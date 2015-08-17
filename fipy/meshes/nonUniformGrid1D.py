#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "nonUniformGrid1D.py"
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

from fipy.tools import parallelComm

from fipy.meshes.mesh1D import Mesh1D
from fipy.meshes.builders import _NonuniformGrid1DBuilder
from fipy.meshes.representations.gridRepresentation import _Grid1DRepresentation
from fipy.meshes.topologies.gridTopology import _Grid1DTopology

__all__ = ["NonUniformGrid1D"]

class NonUniformGrid1D(Mesh1D):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = NonUniformGrid1D(nx = 3)
        >>> print mesh.cellCenters
        [[ 0.5  1.5  2.5]]
         
        >>> mesh = NonUniformGrid1D(dx = (1, 2, 3))
        >>> print mesh.cellCenters
        [[ 0.5  2.   4.5]]
         
        >>> mesh = NonUniformGrid1D(nx = 2, dx = (1, 2, 3))
        Traceback (most recent call last):
        ...
        IndexError: nx != len(dx)

    """
    def __init__(self, dx=1., nx=None, overlap=2, 
                 communicator=parallelComm,
                 _BuilderClass=_NonuniformGrid1DBuilder,
                 _RepresentationClass=_Grid1DRepresentation,
                 _TopologyClass=_Grid1DTopology):

        builder = _BuilderClass()

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

        Mesh1D.__init__(self, vertices, faces, cells, communicator=communicator, 
                        _RepresentationClass=_RepresentationClass, _TopologyClass=_TopologyClass)
        
        self.scale = scale

## pickling

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. Fixed a bug where the following throws
        an error on solve() when nx is a float.

            >>> from fipy import *
            >>> mesh = NonUniformGrid1D(nx=3., dx=(1., 2., 3.))
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var, solver=DummySolver())

        Test for ticket https://github.com/usnistgov/fipy/issues/364.

            >>> from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D
            >>> m = NonUniformGrid1D(nx=9, overlap=1)
            >>> print min(m.x) == 0.5 # doctest: +SERIAL
            True
            >>> print min(m.x) == 3.5 # doctest: +PROCESSOR_1_OF_2
            True
            >>> print min(m.x) == 5.5 # doctest: +PROCESSOR_2_OF_3
            True
            
        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
