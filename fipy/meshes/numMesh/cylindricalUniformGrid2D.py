#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "CylindricalUniformGrid2D.py"
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
2D cylindrical rectangular Mesh with constant spacing in x and constant spacing in y
"""
__docformat__ = 'restructuredtext'

from fipy import numerix
from fipy.meshes.numMesh.uniformGrid2D import UniformGrid2D
from fipy.tools import parallel

class CylindricalUniformGrid2D(UniformGrid2D):
    r"""
    Creates a 2D cylindrical grid in the radial and axial directions,
    appropriate for axial symmetry.
    """
    def __init__(self, dx=1., dy=1., nx=1, ny=1, origin=((0,),(0,)), overlap=2, communicator=parallel):
        UniformGrid2D.__init__(self, dx=dx, dy=dy, nx=nx, ny=ny, origin=origin, overlap=overlap, communicator=communicator)

    def _getAreaProjections(self):
        return self._getAreaProjectionsPy()
            
    def _getFaceAreas(self):
        faceAreas = numerix.zeros(self.numberOfFaces, 'd')
        faceAreas[:self.numberOfHorizontalFaces] = self.dx
        faceAreas[self.numberOfHorizontalFaces:] = self.dy
        return faceAreas * self.getFaceCenters()[0]
        
    def getCellVolumes(self):
        return numerix.ones(self.numberOfCells, 'd') * self.dx * self.dy * self.getCellCenters()[0]

    def _getCellAreas(self):
        areas = numerix.ones((4, self.numberOfCells), 'd')
        areas[0] = self.dx * self.getCellCenters()[0]
        areas[1] = self.dy * (self.getCellCenters()[0] + self.dx / 2)
        areas[2] = self.dx * self.getCellCenters()[0]
        areas[3] = self.dy * (self.getCellCenters()[0] - self.dx / 2)
        return areas

    def _translate(self, vector):
        return CylindricalUniformGrid2D(dx = self.args['dx'], nx = self.args['nx'], 
                                        dy = self.args['dy'], ny = self.args['ny'], 
                                        origin=numerix.array(self.args['origin']) + vector,
                                        overlap=self.args['overlap'])

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> from fipy.tools import parallel

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> mesh = CylindricalUniformGrid2D(nx = nx, ny = ny, dx = dx, dy = dy)     
            
            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 
            ...                            2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 
            ...                            1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array([[dx], [dy]])
            >>> print parallel.procID > 0 or numerix.allequal(vertices, mesh._createVertices())
            True
        
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 
            ...                         0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 
            ...                         4, 1, 2, 3, 8, 5, 6, 7)))
            >>> print parallel.procID > 0 or numerix.allequal(faces, mesh._createFaces())
            True

            >>> cells = numerix.array(((0,  1,  2,  3,  4,  5),
            ...                       (10, 11, 12, 14, 15, 16),
            ...                       ( 3,  4,  5,  6,  7,  8),
            ...                       ( 9, 10, 11, 13, 14, 15)))
            >>> print parallel.procID > 0 or numerix.allequal(cells, mesh._createCells())
            True

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> print parallel.procID > 0 or  numerix.allequal(externalFaces, 
            ...                                                numerix.nonzero(mesh.getExteriorFaces()))
            True

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print parallel.procID > 0 or numerix.allequal(internalFaces, 
            ...                                               numerix.nonzero(mesh.getInteriorFaces()))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  1,  2, 0,  1,  2,  3,  4, 
            ...                                   5,  0,  0, 1,  2,  3,  3,  4,  5),
            ...                                 (-1, -1, -1, 3,  4,  5, -1, -1, 
            ...                                  -1, -1,  1, 2, -1, -1,  4,  5, -1)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            True
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy))
            >>> if parallel.procID == 0: 
            ...     faceAreas = faceAreas * mesh.getFaceCenters()[0]
            >>> print parallel.procID > 0 or numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:]) / 2.
            >>> print parallel.procID > 0 or numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., 
            ...                               -1., 1., 1., 1., -1., 1., 1., 1.),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 
            ...                               1., 0, 0, 0, 0, 0, 0, 0, 0)))
            >>> print parallel.procID > 0 or numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToFaceOrientations = numerix.array(((1,  1,  1, -1, -1, -1), 
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1,  1,  1,  1,  1,  1),
            ...                                         (1, -1, -1,  1, -1, -1)))
            >>> print parallel.procID > 0 or numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            True
                                             
            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy))
            >>> if parallel.procID == 0:
            ...     cellVolumes = cellVolumes * mesh.getCellCenters()[0]
            >>> print numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2.,    dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2.,    dy/2.,    dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> print numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.,
            ...                                dx / 2., dx, dx,
            ...                                dx / 2.))
            >>> print parallel.procID > 0 or numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            True
            
            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2.),
            ...                                         (     -1,      -1,      -1, dy / 2., dy / 2., dy / 2.,      -1,      -1,      -1,      -1, dx / 2., dx / 2.,      -1,      -1, dx / 2., dx / 2.,      -1)), -1)
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print parallel.procID > 0 or numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., 
            ...                             -1., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0, 0, 0, 0, 0, 0, 0, 0, 0, -1., 1., 
            ...                             1., 1., -1., 1., 1., 1.)))
            >>> print parallel.procID > 0 or numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents2 = numerix.zeros((2, 17), 'd')
            >>> print parallel.procID > 0 or numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1,  0,  1,  2),
            ...                                   ( 1,  2, -1,  4,  5, -1),
            ...                                   ( 3,  4,  5, -1, -1, -1),
            ...                                   (-1,  0,  1, -1,  3,  4)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            True

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2.,      dy,      dy,      dy),
            ...                                         (     dx,      dx, dx / 2.,      dx,      dx, dx / 2.),
            ...                                         (     dy,      dy,      dy, dy / 2., dy / 2., dy / 2.),
            ...                                         (dx / 2.,      dx,      dx, dx / 2.,      dx,      dx)), -1)
            >>> print parallel.procID > 0 or numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellNormals = numerix.array(((( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               (-1, -1, -1, -1, -1, -1)),
            ...                              ((-1, -1, -1, -1, -1, -1),
            ...                               ( 0,  0,  0,  0,  0,  0),
            ...                               ( 1,  1,  1,  1,  1,  1),
            ...                               ( 0,  0,  0,  0,  0,  0))))
            >>> print parallel.procID > 0 or numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellAreaProjections = numerix.array((((0,) * 6, (dy,) * 6, (0,) * 6, (-dy,) * 6),
            ...                                      ((-dx,) * 6, (0,) * 6, (dx,) * 6, (0,) * 6)))
            
            >>> if parallel.procID == 0:
            ...     cellAreaProjections[:,0] = cellAreaProjections[:,0] * mesh.getCellCenters()[0]
            ...     cellAreaProjections[:,1] = cellAreaProjections[:,1] * (mesh.getCellCenters()[0] + mesh.dx / 2.)
            ...     cellAreaProjections[:,2] = cellAreaProjections[:,2] * mesh.getCellCenters()[0]
            ...     cellAreaProjections[:,3] = cellAreaProjections[:,3] * (mesh.getCellCenters()[0] - mesh.dx / 2.)
            >>> print parallel.procID > 0 or numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellVertexIDs = MA.masked_array(((5, 6, 7, 9, 10, 11),
            ...                                  (4, 5, 6, 8, 9, 10),
            ...                                  (1, 2, 3, 5, 6, 7),
            ...                                  (0, 1, 2, 4, 5, 6)), -1000)

            >>> print parallel.procID > 0 or numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            True

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allclose(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            True
            
            >>> faceVertexIDs = [[ 0, 1, 2, 4, 5, 6, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7],
            ...                  [ 1, 2, 3, 5, 6, 7, 9, 10, 11, 4, 5, 6, 7, 8, 9, 10, 11]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getFaceVertexIDs(), faceVertexIDs)
            True

            >>> mesh = CylindricalUniformGrid2D(nx=3)
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[0],
            ...                                               [0, 1, 2, 0, 1, 2, 0, 0, 1, 2])
            True
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[1],
            ...                                               [0, 1, 2, 0, 1, 2, 0, 1, 2, 2])
            True
            >>> faceCellIDs = [[0, 1, 2, 0, 1, 2, 0, 0, 1, 2],
            ...                [-1, -1, -1, -1, -1, -1, -1, 1, 2, -1]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh.getFaceCellIDs().filled(-1),
            ...                                               faceCellIDs)
            True

            >>> mesh = CylindricalUniformGrid2D(ny=3)
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[0],
            ...                                               [0, 0, 1, 2, 0, 0, 1, 1, 2, 2])
            True
            >>> print parallel.procID > 0 or numerix.allequal(mesh._getAdjacentCellIDs()[1],
            ...                                               [0, 1, 2, 2, 0, 0, 1, 1, 2, 2])
            True
            >>> faceCellIDs = [[0, 0, 1, 2, 0, 0, 1, 1, 2, 2],
            ...                [-1, 1, 2, -1, -1, -1, -1, -1, -1, -1]]
            >>> print parallel.procID > 0 or numerix.allequal(mesh.getFaceCellIDs().filled(-1),
            ...                                               faceCellIDs)
            True

        Following test added to change nx, ny argment to integer when its a float to prevent
        warnings from the solver.

            >>> from fipy import *
            >>> mesh = CylindricalUniformGrid2D(nx=3., ny=3., dx=1., dy=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
