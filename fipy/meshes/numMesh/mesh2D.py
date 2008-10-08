#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh2D.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 10/8/08 {3:01:36 PM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

"""
Generic mesh class using numerix to do the calculations

Meshes contain cells, faces, and vertices.

This is built for a non-mixed element mesh.
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.meshes.numMesh.mesh import Mesh
from fipy.tools import numerix
from fipy.tools import vector


def _orderVertices(vertexCoords, vertices):
    coordinates = numerix.take(vertexCoords, vertices)
    centroid = numerix.add.reduce(coordinates) / coordinates.shape[0]
    coordinates = coordinates - centroid
    coordinates = numerix.where(coordinates == 0, 1.e-100, coordinates) ## to prevent division by zero
    angles = numerix.arctan(coordinates[:, 1] / coordinates[:, 0]) + numerix.where(coordinates[:, 0] < 0, numerix.pi, 0) ## angles go from -pi / 2 to 3*pi / 2
    sortorder = numerix.argsort(angles)
    return numerix.take(vertices, sortorder)

class Mesh2D(Mesh):
    def _calcFaceAreas(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs)
        tangent = faceVertexCoords[:,1] - faceVertexCoords[:,0]
        self.faceAreas = numerix.sqrtDot(tangent, tangent)

    def _calcFaceNormals(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        self.faceNormals = t1.copy()
        self.faceNormals[:,0] = -t1[:,1] / numerix.sqrt(t1[:,1]**2 + t1[:,0]**2)
        self.faceNormals[:,1] = t1[:,0] / numerix.sqrt(t1[:,1]**2 + t1[:,0]**2)
        
        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.cellDistanceVectors) < 0)
        self.faceNormals = self.faceNormals * orientation[..., numerix.NewAxis]


    def _calcFaceTangents(self):
        tmp = numerix.transpose(numerix.array((-self.faceNormals[:,1], self.faceNormals[:,0])))
        mag = numerix.sqrtDot(tmp, tmp)
        self.faceTangents1 = tmp / mag[:,numerix.NewAxis]
        self.faceTangents2 = numerix.zeros(self.faceTangents1.shape, 'd')

    def _calcHigherOrderScalings(self):
        self.scale['area'] = self.scale['length']
        self.scale['volume'] = self.scale['length']**2

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh2D(newCoords, self.faceVertexIDs, self.cellFaceIDs)
        return newmesh

    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh2D(newCoords, self.faceVertexIDs, self.cellFaceIDs)
        return newmesh

    def _concatenate(self, other, smallNumber):
        return Mesh2D(**self._getAddedMeshValues(other, smallNumber))

    def _getOrderedCellVertexIDs_old(self):
        cellFaceVertexIDs = numerix.take(self.faceVertexIDs, self.cellFaceIDs)
        cellVertexIDs = numerix.reshape(cellFaceVertexIDs, (len(self.getCellCenters()), numerix.size(cellFaceVertexIDs) / len(self.getCellCenters())))
        cellVertexIDs = numerix.sort(cellVertexIDs)
        cellVertexIDs = cellVertexIDs[:, ::2]
        orderedCellVertexIDs = numerix.zeros(cellVertexIDs.shape)
        for i in range(cellVertexIDs.shape[0]):
            orderedCellVertexIDs[i, :] = _orderVertices(self.getVertexCoords(), cellVertexIDs[i, :])
        return orderedCellVertexIDs
                                                     
    def _getOrderedCellVertexIDs(self):
        from fipy.tools.numerix import take
        NFac = self._getMaxFacesPerCell()
        # numpy 1.1's MA.take doesn't like FlatIter. Call ravel() instead.
        cellVertexIDs0 = take(self._getFaceVertexIDs()[:,0], self._getCellFaceIDs().ravel())
        cellVertexIDs1 = take(self._getFaceVertexIDs()[:,1], self._getCellFaceIDs().ravel())
        cellVertexIDs = MA.where(self.cellToFaceOrientations.ravel() > 0,
                             cellVertexIDs0, cellVertexIDs1)

        cellVertexIDs = MA.reshape(cellVertexIDs, (-1, NFac))
        return cellVertexIDs
    
    def _getNonOrthogonality(self):
        exteriorFaceArray = numerix.zeros((self.faceCellIDs.shape[0],))
        numerix.put(exteriorFaceArray, self.getExteriorFaces(), 1)
        unmaskedFaceCellIDs = MA.filled(self.faceCellIDs, 0) ## what we put in for the "fill" doesn't matter because only exterior faces have anything masked, and exterior faces have their displacement vectors set to zero.
        ## if it's an exterior face, make the "displacement vector" equal to zero so the cross product will be zero.
        faceDisplacementVectors = numerix.where(numerix.array(zip(exteriorFaceArray, exteriorFaceArray)), 0.0, numerix.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 1]) - numerix.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 0]))
        faceCrossProducts = (faceDisplacementVectors[:, 0] * self.faceNormals[:, 1]) - (faceDisplacementVectors[:, 1] * self.faceNormals[:, 0])
        faceDisplacementVectorLengths = numerix.maximum(((faceDisplacementVectors[:, 0] ** 2) + (faceDisplacementVectors[:, 1] ** 2)) ** 0.5, 1.e-100)
        faceWeightedNonOrthogonalities = abs(faceCrossProducts / faceDisplacementVectorLengths) * self.faceAreas
        cellFaceWeightedNonOrthogonalities = numerix.take(faceWeightedNonOrthogonalities, self.cellFaceIDs)
        cellFaceAreas = numerix.take(self.faceAreas, self.cellFaceIDs)
        cellTotalWeightedValues = numerix.add.reduce(cellFaceWeightedNonOrthogonalities, axis = 1)
        cellTotalFaceAreas = numerix.add.reduce(cellFaceAreas, axis = 1)
        return (cellTotalWeightedValues / cellTotalFaceAreas)

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> vertices = numerix.array(((0., 0.), (1., 0.), (2., 0.), (3., 0.),
            ...                           (0., 1.), (1., 1.), (2., 1.), (3., 1.),
            ...                           (0., 2.), (1., 2.), (2., 2.), (3., 2.),
            ...                           (4., 1.)))
            >>> vertices *= numerix.array((dx, dy))
            >>> faces = numerix.array(((1, 0), (2, 1), (3, 2),
            ...                        (4, 5), (5, 6), (6, 7),
            ...                        (8, 9), (9, 10), (10, 11),
            ...                        (0, 4), (5, 1), (6, 2), (7, 3),
            ...                        (4, 8), (9, 5), (10, 6), (11, 7),
            ...                        (12, 3), (7, 12), (11, 12)))
            >>> cells = MA.masked_values(((0, 10, 3, 9),
            ...                           (1 , 11, 4, 10),
            ...                           (2, 12, 5, 11),
            ...                           (3, 14, 6, 13),
            ...                           (4, 15, 7, 14),
            ...                           (5, 16, 8, 15),
            ...                           (17, 18, 12, -1),
            ...                           (18, 19 ,16, -1)), -1)
                
            >>> mesh = Mesh2D(vertexCoords = vertices, faceVertexIDs = faces, cellFaceIDs = cells)
            
            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9, 13, 17, 19))
            >>> numerix.allequal(externalFaces, mesh.getExteriorFaces())
            1

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 12, 14, 15, 16, 18))
            >>> numerix.allequal(internalFaces, mesh.getInteriorFaces())
            1

            >>> faceCellIds = MA.masked_values(((0, -1), (1, -1), (2, -1),
            ...                                 (0, 3), (1, 4), (2, 5),
            ...                                 (3, -1), (4, -1), (5, -1),
            ...                                 (0, -1), (0, 1), (1, 2), (2, 6),
            ...                                 (3, -1), (3, 4), (4, 5), (5, 7),
            ...                                 (6, -1), (6, 7), (7, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy,
            ...                            numerix.sqrt(dx**2 + dy**2), dx, numerix.sqrt(dx**2 + dy**2)))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, faces)
            >>> faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array(((0., -1.), (0., -1.), (0., -1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0),
            ...                              (dy / numerix.sqrt(dx**2 + dy**2), -dx / numerix.sqrt(dx**2 + dy**2)),
            ...                              (0., 1.),
            ...                              (dy / numerix.sqrt(dx**2 + dy**2), dx / numerix.sqrt(dx**2 + dy**2))))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, 1, 1, 1), (1, 1, 1, -1), (1, 1, 1, -1),
            ...                                            (-1, 1, 1, 1), (-1, 1, 1, -1), (-1, 1, 1, -1),
            ...                                            (1, 1, -1, 0), (-1, 1, -1, 0)), 0)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy / 2., dx*dy / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2.,dy/2.), (3.*dx/2.,dy/2.), (5.*dx/2.,dy/2.), 
            ...                              (dx/2.,3.*dy/2.), (3.*dx/2.,3.*dy/2.), (5.*dx/2.,3.*dy/2.),
            ...                              (3.*dx+dx/3.,2.*dy/3.), (3.*dx+dx/3.,4.*dy/3.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dy / 2., dy / 2.), (dy / 2., dy / 2.), (dy / 2., dy / 2.),
            ...                                         (dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
            ...                                         (numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1),
            ...                                         (numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2), numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)),
            ...                                         (numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = numerix.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                numerix.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
            ...                                dx / 2., dx, dx,
            ...                                numerix.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
            ...                                numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2),
            ...                                2. * dy / 3.,
            ...                                numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,numerix.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = numerix.array(((1., 0), (1., 0),(1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.),
            ...                            (dx / numerix.sqrt(dx**2 +dy**2), dy / numerix.sqrt(dx**2 +dy**2)),
            ...                            (-1, 0.),
            ...                            (-dx / numerix.sqrt(dx**2 +dy**2), dy / numerix.sqrt(dx**2 +dy**2))))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = numerix.array(((0., 0), (0., 0),(0., 0),
            ...                            (-0., 0), (-0., 0),(-0., 0),
            ...                            (-0., 0), (-0., 0),(-0., 0),
            ...                            (0., -0.), (0., 0.), (0., 0.), (0., 0.),
            ...                            (0., -0.), (0., 0.), (0., 0.), (0., 0.),
            ...                            (0., 0), (0., 0),(0., 0)))
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1, 1, 3, -1),
            ...                                   (-1, 2, 4, 0),
            ...                                   (-1, 6, 5, 1),
            ...                                   (0, 4, -1, -1),
            ...                                   (1, 5, -1, 3),
            ...                                   (2, 7, -1, 4),
            ...                                   (-1, 7, 2, -1),
            ...                                   (6, -1, 5, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> d1 = numerix.sqrt((5. * dx / 6.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> cellToCellDistances = MA.masked_values(((dy / 2., dx, dy, dx / 2.),
            ...                                         (dy / 2., dx, dy, dx),
            ...                                         (dy / 2., d1, dy, dx),
            ...                                         (dy, dx, dy / 2., dx / 2.),
            ...                                         (dy, dx, dy / 2., dx),
            ...                                         (dy, d1, dy / 2., dx),
            ...                                         (d2, 2. * dy / 3., d1, -1),
            ...                                         (2. * dy / 3., d2, d1, -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1, 2, 3, 4, 5, 6, 7))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> nx = dy / numerix.sqrt(dx**2 + dy**2)
            >>> ny = dx / numerix.sqrt(dx**2 + dy**2)
            >>> cellNormals = MA.masked_values(((( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 (( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 (( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 (( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 (( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 (( 0,  -1), ( 1,  0), ( 0, 1), (   -1,     0)),
            ...                                 ((nx, -ny), ( 0,  1), (-1, 0), (-1000, -1000)),
            ...                                 (( 0,  -1), (nx, ny), (-1, 0), (-1000, -1000))), -1000 )
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> vv = numerix.array(((0, -dx), (dy, 0), (0, dx), (-dy, 0)))
            >>> area = numerix.sqrt(dx**2 + dy**2)
            >>> cellAreaProjections = MA.masked_values((vv,vv,vv,vv,vv,vv,
            ...                                        ((nx * area, -ny * area), (0, dx), (-dy, 0), (-1000, -1000)),
            ...                                        ((0, -dx), (nx * area, ny * area), (-dy, 0), (-1000, -1000))), -1000)
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_array(((5, 4, 1, 0),
            ...                                  (6, 5, 2, 1),
            ...                                  (7, 6, 3, 2),
            ...                                  (9, 8, 5, 4),
            ...                                  (10, 9, 6, 5),
            ...                                  (11, 10, 7, 6),
            ...                                  (12, 7, 3, -1000),
            ...                                  (12, 11, 7, -1000)), -1000)

            >>> numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            1
            

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            1
        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
