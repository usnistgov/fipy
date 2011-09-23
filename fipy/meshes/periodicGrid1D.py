#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "periodicGrid1D.py"
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
Peridoic 1D Mesh
"""
__docformat__ = 'restructuredtext'

from grid1D import Grid1D
from fipy.tools import numerix
from fipy.tools.decorators import getsetDeprecated
from fipy.meshes.builders import PeriodicGrid1DBuilder

class PeriodicGrid1D(Grid1D):
    """
    
        >>> from fipy.tools import parallel

    Creates a Periodic grid mesh.
        
        >>> mesh = PeriodicGrid1D(dx = (1, 2, 3))

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(numerix.nonzero(mesh.exteriorFaces)[0], 
        ...                            [3]))
        True

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh.faceCellIDs.filled(-999),
        ...                            [[2, 0, 1, 2],
        ...                             [0, 1, 2, -999]]))
        True

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._cellDistances,
        ...                            [ 2., 1.5, 2.5, 1.5]))
        True

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._cellToCellDistances,
        ...                            [[ 2.,   1.5,  2.5],
        ...                             [ 1.5,  2.5,  2. ]]))
        True
        
        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._faceNormals,
        ...                            [[ 1.,  1.,  1.,  1.]]))
        True

        >>> print (parallel.procID > 0 
        ...        or numerix.allclose(mesh._cellVertexIDs,
        ...                            [[1, 2, 2],
        ...                             [0, 1, 0]]))
        True
    """
    def __init__(self, dx = 1., nx = None, overlap=2):
        Grid1D.__init__(self, dx = dx, nx = nx, overlap=overlap,
                        BuilderClass=PeriodicGrid1DBuilder)
        self.nonPeriodicCellFaceIDs = numerix.array(super(Grid1D, self).cellFaceIDs)
        self._makePeriodic()

    def _makePeriodic(self):
        if self.occupiedNodes == 1:
            self._connectFaces(numerix.nonzero(self.facesLeft),
                               numerix.nonzero(self.facesRight))

    @getsetDeprecated
    def _getGlobalOverlappingCellIDs(self):
        return self._globalOverlappingCellIDs

    @property
    def _globalOverlappingCellIDs(self):
        return super(PeriodicGrid1D, self)._globalOverlappingCellIDs % self.args['nx']

    @property
    def cellCenters(self):
        """Defined outside of a geometry class since we need the `CellVariable`
        version of `cellCenters`; that is, the `cellCenters` defined in
        fipy.meshes.mesh and not in any geometry (since a `CellVariable` requires
        a reference to a mesh)."""
        return super(PeriodicGrid1D, self).cellCenters \
                % numerix.sum(self.globalNumberOfCells * self.args['dx'])

    def _translate(self, vector):
        """
        Test for ticket:298.

        >>> from fipy import *
        >>> m = PeriodicGrid1D(nx=2) + [[-1]]
        >>> print CellVariable(mesh=m, value=m.cellCenters[0])
        [-0.5  0.5]
        
        """
        newCoords = self.vertexCoords + vector
        newmesh = self.__class__(**self.args)
        from fipy.meshes.mesh1D import Mesh1D
        Mesh1D.__init__(newmesh, newCoords, numerix.array(self.faceVertexIDs), self.nonPeriodicCellFaceIDs)
        newmesh._makePeriodic()
        return newmesh
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
