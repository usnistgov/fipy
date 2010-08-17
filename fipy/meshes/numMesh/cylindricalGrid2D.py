#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "CylindricalGrid2D.py"
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
2D rectangular Mesh
"""
__docformat__ = 'restructuredtext'


from fipy.tools import numerix

from fipy.meshes.numMesh.grid2D import Grid2D
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import parallel

class CylindricalGrid2D(Grid2D):
    """
    Creates a 2D cylindrical grid mesh with horizontal faces numbered
    first and then vertical faces.
    """
    def __init__(self, dx=1., dy=1., nx=None, ny=None, origin=((0.,), (0.,)), overlap=2, communicator=parallel):
        scale = PhysicalField(value=1, unit=PhysicalField(value=dx).getUnit())
        self.origin = PhysicalField(value=origin)
        self.origin /= scale
        
        Grid2D.__init__(self, dx=dx, dy=dy, nx=nx, ny=ny, overlap=overlap, communicator=communicator)

        self.args['origin'] = self.origin

    def _getFaceAreas(self):
        return Grid2D._getFaceAreas(self) * self.getFaceCenters()[0]

    def getCellVolumes(self):
        return Grid2D.getCellVolumes(self) * self.getCellCenters()[0]

    def _translate(self, vector):
        return CylindricalGrid2D(dx=self.args['dx'], nx=self.args['nx'], 
                                 dy=self.args['dy'], ny=self.args['ny'], 
                                 origin=self.args['origin'] + vector,
                                 overlap=self.args['overlap'])

    def __mul__(self, factor):
        if numerix.shape(factor) is ():
            factor = numerix.resize(factor, (2,1))
        
        return CylindricalGrid2D(dx=self.args['dx'] * numerix.array(factor[0]), nx=self.args['nx'], 
                                 dy=self.args['dy'] * numerix.array(factor[1]), ny=self.args['ny'], 
                                 origin=self.args['origin'] * factor,
                                 overlap=self.args['overlap'])

    def getVertexCoords(self):
        return self.vertexCoords + self.origin

    def getCellCenters(self):
        return Grid2D.getCellCenters(self) + self.origin

    def getFaceCenters(self):
        return self.faceCenters + self.origin
    
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.

            >>> from fipy.tools import parallel

            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> mesh = CylindricalGrid2D(nx = nx, ny = ny, dx = dx, dy = dy)     
            
            >>> vertices = numerix.array(((0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> print parallel.procID > 0 or numerix.allequal(vertices, mesh._createVertices())
            True
        
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8, 9, 10, 0, 5, 6, 7, 4, 9, 10, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 4, 1, 2, 3, 8, 5, 6, 7)))
            >>> print parallel.procID > 0 or numerix.allequal(faces, mesh._createFaces())
            True

            >>> cells = numerix.array(((0, 1, 2, 3, 4, 5),
            ...                        (10, 11, 12, 14, 15, 16),
            ...                        (3, 4, 5, 6, 7, 8),
            ...                        (9, 10, 11, 13, 14, 15)))
            >>> print parallel.procID > 0 or numerix.allequal(cells, mesh._createCells())
            True

            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9 , 12, 13, 16))
            >>> print parallel.procID > 0 or numerix.allequal(externalFaces, 
            ...                                               numerix.nonzero(mesh.getExteriorFaces()))
            True

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 14, 15))
            >>> print parallel.procID > 0 or numerix.allequal(internalFaces, 
            ...                                               numerix.nonzero(mesh.getInteriorFaces()))
            True

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values(((0, 1, 2, 0, 1, 2, 3, 4, 5, 0, 0, 1, 2, 3, 3, 4, 5),
            ...                                 (-1, -1, -1, 3, 4, 5, -1, -1, -1, -1, 1, 2, -1, -1, 4, 5, -1)), -1)
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

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0., 0.)))
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

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2.)))
            >>> print numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            True

            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dy / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2., dx / 2.),
            ...                                         (-1, -1, -1, dy / 2., dy / 2., dy / 2., -1, -1, -1, -1, dx / 2., dx / 2., -1, -1, dx / 2., dx / 2., -1)), -1)
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
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
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print parallel.procID > 0 or numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            True

            >>> areaProjections = faceNormals * faceAreas
            >>> print parallel.procID > 0 or numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 1., 1., 1., -1., 1., 1., 1.)))
            >>> print parallel.procID > 0 or numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            True

            >>> tangents2 = numerix.zeros((2, 17), 'd')
            >>> print parallel.procID > 0 or numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1, 0, 1, 2),
            ...                                   (1, 2, -1, 4, 5, -1),
            ...                                   (3, 4, 5, -1, -1, -1),
            ...                                   (-1, 0, 1, -1, 3, 4)), -1)
            >>> print parallel.procID > 0 or numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            True

            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2.,      dy,      dy,      dy),
            ...                                         (     dx,      dx, dx / 2.,      dx,      dx, dx / 2.),
            ...                                         (     dy,      dy,      dy, dy / 2., dy / 2., dy / 2.),
            ...                                         (dx / 2.,      dx,      dx, dx / 2.,      dx,      dx)), -1)
            >>> print parallel.procID > 0 or numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            True

            >>> interiorCellIDs = numerix.array(())
            >>> print parallel.procID > 0 or numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            True

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5))
            >>> print parallel.procID > 0 or numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
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

            >>> cellAreaProjections = numerix.array((((  0,  0,  0,  0,  0,  0),
            ...                                       ( dy, dy, dy, dy, dy, dy),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       (-dy,-dy,-dy,-dy,-dy,-dy)),
            ...                                      ((-dx,-dx,-dx,-dx,-dx,-dx),
            ...                                       (  0,  0,  0,  0,  0,  0),
            ...                                       ( dx, dx, dx, dx, dx, dx),
            ...                                       (  0,  0,  0,  0,  0,  0))))
            >>> if parallel.procID == 0:
            ...     cellAreaProjections[:,0] = cellAreaProjections[:,0] * mesh.getCellCenters()[0]
            ...     cellAreaProjections[:,1] = cellAreaProjections[:,1] * (mesh.getCellCenters()[0] + mesh.dx / 2.)
            ...     cellAreaProjections[:,2] = cellAreaProjections[:,2] * mesh.getCellCenters()[0]
            ...     cellAreaProjections[:,3] = cellAreaProjections[:,3] * (mesh.getCellCenters()[0] - mesh.dx / 2.)
            >>> print parallel.procID > 0 or numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            True

            >>> cellVertexIDs = MA.masked_values(((5, 6, 7, 9, 10, 11),
            ...                                   (4, 5, 6, 8,  9, 10),
            ...                                   (1, 2, 3, 5,  6,  7),
            ...                                   (0, 1, 2, 4,  5,  6)), -1000)

            >>> print parallel.procID > 0 or numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            True

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')            
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allclose(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            True

            >>> mesh = CylindricalGrid2D(dx=(1., 2.), dy=(1.,)) + ((1.,),(0.,))
            >>> print mesh.getCellCenters()
            [[ 1.5  3. ]
             [ 0.5  0.5]]
            >>> print mesh.getCellVolumes()
            [ 1.5  6. ]

            
        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
