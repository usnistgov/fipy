#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "periodicGrid3D.py"
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
3D periodic rectangular Mesh
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools import parallelComm
from fipy.meshes.nonUniformGrid3D import NonUniformGrid3D

__all__ = ["PeriodicGrid3D", "PeriodicGrid3DLeftRight", "PeriodicGrid3DTopBottom",
           "PeriodicGrid3DFrontBack", "PeriodicGrid3DLeftRightTopBottom",
           "PeriodicGrid3DLeftRightFrontBack", "PeriodicGrid3DTopBottomFrontBack"]

class _BasePeriodicGrid3D(NonUniformGrid3D):
    def __init__(self, dx=1., dy=1., dz=1., nx=None, ny=None, nz=None, overlap=2, communicator=parallelComm, *args, **kwargs):
        super(_BasePeriodicGrid3D, self).__init__(dx=dx, dy=dy, dz=dz, nx=nx, ny=ny, nz=nz, overlap=overlap, communicator=communicator, *args, **kwargs)
        self._nonPeriodicCellVertexIDs = super(_BasePeriodicGrid3D, self)._cellVertexIDs
        self._orderedCellVertexIDs_data = super(_BasePeriodicGrid3D, self)._orderedCellVertexIDs        
        self._nonPeriodicCellFaceIDs = numerix.array(super(_BasePeriodicGrid3D, self).cellFaceIDs)
        self._makePeriodic()

    @property
    def _cellVertexIDs(self):
        return self._nonPeriodicCellVertexIDs

    def _translate(self, vector):
        """
        Test for ticket:298.

        >>> from fipy import *
        >>> m = PeriodicGrid3DLeftRight(nx=2, ny=2, nz=1) + [[-1], [0], [0]]
        >>> orderedCellVertexIDs = [[13, 14, 16, 17],
        ...                         [12, 13, 15, 16],
        ...                         [10, 11, 13, 14],
        ...                         [9, 10, 12, 13],
        ...                         [4, 5, 7, 8],
        ...                         [3, 4, 6, 7],
        ...                         [1, 2, 4, 5],
        ...                         [0, 1, 3, 4]]
        >>> print numerix.allclose(m._orderedCellVertexIDs, orderedCellVertexIDs)  # doctest: +PROCESSOR_0
        True
        >>> print CellVariable(mesh=m, value=m.cellCenters[0])
        [-0.5  0.5 -0.5  0.5]
        """
        newCoords = self.vertexCoords + vector
        newmesh = self.__class__(**self.args)
        from fipy.meshes.mesh import Mesh
        Mesh.__init__(newmesh, newCoords, self.faceVertexIDs, self._nonPeriodicCellFaceIDs, communicator=self.communicator)
        newmesh._makePeriodic()
        return newmesh

class PeriodicGrid3D(_BasePeriodicGrid3D):
    """
    Creates a periodic3D grid mesh with horizontal faces numbered
    first and then vertical faces. Vertices and cells are numbered 
    in the usual way.

        >>> from fipy import numerix

        >>> mesh = PeriodicGrid3D(dx=1., dy=0.5, dz=2., nx=2, ny=2, nz=1)
        >>> print numerix.allclose(numerix.nonzero(mesh.exteriorFaces)[0],
        ...                        [4, 5, 6, 7, 12, 13, 16, 19])  # doctest: +PROCESSOR_0
        True

        >>> print numerix.allclose(mesh.faceCellIDs.filled(-1),
        ...                        [[0, 1, 2, 3, 0, 1, 2, 3, 2, 3, 
        ...                          0, 1, 2, 3, 1, 0, 1, 3, 2, 3],
        ...                         [0, 1, 2, 3, -1, -1, -1, -1, 0, 1, 
        ...                          2, 3, -1, -1, 0, 1, -1, 2, 3, -1]]) # doctest: +PROCESSOR_0
        True


        >>> print numerix.allclose(mesh._cellDistances,
        ...                        [2., 2., 2., 2., 1., 1., 1., 1., 0.5, 0.5,
        ...                         0.5, 0.5, 0.25, 0.25, 1., 1., 0.5, 1., 1., 0.5]) # doctest: +PROCESSOR_0
        True

        >>> print numerix.allclose(mesh.cellFaceIDs,
        ...                        [[14,15,17,18],
        ...                         [15, 14, 18,17],
        ...                         [8, 9, 10, 11],
        ...                         [10, 11, 8, 9],
        ...                         [0, 1, 2, 3],
        ...                         [0, 1, 2, 3]]) # doctest: +PROCESSOR_0
        True

        >>> print numerix.allclose(mesh._cellToCellDistances,
        ...                        [[1., 1., 1., 1.],
        ...                         [1., 1., 1., 1.],        
        ...                         [0.5, 0.5, 0.5, 0.5],
        ...                         [0.5, 0.5, 0.5, 0.5],
        ...                         [2., 2., 2., 2.],        
        ...                         [2., 2., 2., 2.]]) # doctest: +PROCESSOR_0
        True

        >>> normals = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
        ...            [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        ...            [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

        >>> print numerix.allclose(mesh.faceNormals, normals) # doctest: +PROCESSOR_0
        True

        >>> print numerix.allclose(mesh._cellVertexIDs,
        ...                        [[13, 14, 16, 17],
        ...                         [12, 13, 15, 16],
        ...                         [10, 11, 13, 14],
        ...                         [9, 10, 12, 13],
        ...                         [4, 5, 7, 8],
        ...                         [3, 4, 6, 7],
        ...                         [1, 2, 4, 5],
        ...                         [0, 1, 3, 4]]) # doctest: +PROCESSOR_0
        True
    """

    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft), 
                           numerix.nonzero(self.facesRight))
        self._connectFaces(numerix.nonzero(self.facesBottom), 
                           numerix.nonzero(self.facesTop))
        self._connectFaces(numerix.nonzero(self.facesFront), 
                           numerix.nonzero(self.facesBack))

    def _test(self):
        """
        Test to check that diffusion works correctly by checking that
        the elements one step away for element 0 on the diagonal are
        equal. They wouldn't be equal for a non-periodic grid.

        >>> import fipy as fp
        >>> m = fp.PeriodicGrid3D(nx=3, ny=3, nz=3)
        >>> v = fp.CellVariable(mesh=m)
        >>> v[0] = 1.
        >>> (fp.TransientTerm() == fp.DiffusionTerm()).solve(v, dt=1.)
        >>> assert numerix.allclose(v[13], v[26])
        
        """

        pass

class PeriodicGrid3DLeftRight(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft),
                           numerix.nonzero(self.facesRight))

class PeriodicGrid3DLeftRightTopBottom(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft),
                           numerix.nonzero(self.facesRight))
        self._connectFaces(numerix.nonzero(self.facesBottom),
                           numerix.nonzero(self.facesTop))

class PeriodicGrid3DLeftRightFrontBack(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesLeft),
                           numerix.nonzero(self.facesRight))
        self._connectFaces(numerix.nonzero(self.facesFront),
                           numerix.nonzero(self.facesBack))

class PeriodicGrid3DTopBottom(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesBottom),
                           numerix.nonzero(self.facesTop))

class PeriodicGrid3DTopBottomFrontBack(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesBottom),
                           numerix.nonzero(self.facesTop))
        self._connectFaces(numerix.nonzero(self.facesFront),
                           numerix.nonzero(self.facesBack))

class PeriodicGrid3DFrontBack(_BasePeriodicGrid3D):
    def _makePeriodic(self):
        self._connectFaces(numerix.nonzero(self.facesFront),
                           numerix.nonzero(self.facesBack))


        
def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
