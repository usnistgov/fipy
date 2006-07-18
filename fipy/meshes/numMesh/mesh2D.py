#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh2D.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 3/4/06 {3:55:06 PM} 
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
Generic mesh class using Numeric to do the calculations

Meshes contain cells, faces, and vertices.

This is built for a non-mixed element mesh.
"""
__docformat__ = 'restructuredtext'

import Numeric
import MA
from fipy.meshes.numMesh.mesh import Mesh
from fipy.tools import numerix
from fipy.tools import vector


def _orderVertices(vertexCoords, vertices):
    coordinates = Numeric.take(vertexCoords, vertices)
    centroid = Numeric.add.reduce(coordinates) / coordinates.shape[0]
    coordinates = coordinates - centroid
    coordinates = Numeric.where(coordinates == 0, 1.e-100, coordinates) ## to prevent division by zero
    angles = Numeric.arctan(coordinates[:, 1] / coordinates[:, 0]) + Numeric.where(coordinates[:, 0] < 0, Numeric.pi, 0) ## angles go from -pi / 2 to 3*pi / 2
    sortorder = Numeric.argsort(angles)
    return Numeric.take(vertices, sortorder)

class Mesh2D(Mesh):
    def _calcFaceAreas(self):
        faceVertexCoords = Numeric.take(self.vertexCoords, self.faceVertexIDs)
        tangent = faceVertexCoords[:,1] - faceVertexCoords[:,0]
        self.faceAreas = numerix.sqrtDot(tangent, tangent)
        
    def _calcFaceNormals(self):
        faceVertexCoords = Numeric.take(self.vertexCoords, self.faceVertexIDs)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        self.faceNormals = t1.copy()
        self.faceNormals[:,0] = -t1[:,1] / Numeric.sqrt(t1[:,1]**2 + t1[:,0]**2)
        self.faceNormals[:,1] = t1[:,0] / Numeric.sqrt(t1[:,1]**2 + t1[:,0]**2)

    def _calcFaceTangents(self):
        tmp = Numeric.transpose(Numeric.array((-self.faceNormals[:,1], self.faceNormals[:,0])))
        mag = numerix.sqrtDot(tmp, tmp)
        self.faceTangents1 = tmp / mag[:,Numeric.NewAxis]
        self.faceTangents2 = Numeric.zeros(self.faceTangents1.shape, 'd')

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
        cellFaceVertexIDs = Numeric.take(self.faceVertexIDs, self.cellFaceIDs)
        cellVertexIDs = Numeric.reshape(cellFaceVertexIDs, (len(self.getCellCenters()), Numeric.size(cellFaceVertexIDs) / len(self.getCellCenters())))
        cellVertexIDs = Numeric.sort(cellVertexIDs)
        cellVertexIDs = cellVertexIDs[:, ::2]
        orderedCellVertexIDs = Numeric.zeros(cellVertexIDs.shape)
        for i in range(cellVertexIDs.shape[0]):
            orderedCellVertexIDs[i, :] = _orderVertices(self.getVertexCoords(), cellVertexIDs[i, :])
        return orderedCellVertexIDs
                                                     
    def _getOrderedCellVertexIDs(self):
        from fipy.tools.numerix import MAtake
        NFac = self._getMaxFacesPerCell()
        cellVertexIDs0 = MAtake(self._getFaceVertexIDs()[:,0], self._getCellFaceIDs().flat)
        cellVertexIDs1 = MAtake(self._getFaceVertexIDs()[:,1], self._getCellFaceIDs().flat)
        cellVertexIDs = MA.where(self.cellToFaceOrientations.flat > 0,
                                 cellVertexIDs0,
                                 cellVertexIDs1)

        cellVertexIDs = MA.reshape(cellVertexIDs, (-1, NFac))
        return cellVertexIDs
    
    def _getNonOrthogonality(self):
        exteriorFaceArray = Numeric.zeros((self.faceCellIDs.shape[0],))
        Numeric.put(exteriorFaceArray, self.getExteriorFaces(), 1)
        unmaskedFaceCellIDs = MA.filled(self.faceCellIDs, 0) ## what we put in for the "fill" doesn't matter because only exterior faces have anything masked, and exterior faces have their displacement vectors set to zero.
        ## if it's an exterior face, make the "displacement vector" equal to zero so the cross product will be zero.
        faceDisplacementVectors = Numeric.where(Numeric.array(zip(exteriorFaceArray, exteriorFaceArray)), 0.0, Numeric.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 1]) - Numeric.take(self.getCellCenters(), unmaskedFaceCellIDs[:, 0]))
        faceCrossProducts = (faceDisplacementVectors[:, 0] * self.faceNormals[:, 1]) - (faceDisplacementVectors[:, 1] * self.faceNormals[:, 0])
        faceDisplacementVectorLengths = Numeric.maximum(((faceDisplacementVectors[:, 0] ** 2) + (faceDisplacementVectors[:, 1] ** 2)) ** 0.5, 1.e-100)
        faceWeightedNonOrthogonalities = abs(faceCrossProducts / faceDisplacementVectorLengths) * self.faceAreas
        cellFaceWeightedNonOrthogonalities = Numeric.take(faceWeightedNonOrthogonalities, self.cellFaceIDs)
        cellFaceAreas = Numeric.take(self.faceAreas, self.cellFaceIDs)
        cellTotalWeightedValues = Numeric.add.reduce(cellFaceWeightedNonOrthogonalities, axis = 1)
        cellTotalFaceAreas = Numeric.add.reduce(cellFaceAreas, axis = 1)
        return (cellTotalWeightedValues / cellTotalFaceAreas)

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> vertices = Numeric.array(((0., 0.), (1., 0.), (2., 0.), (3., 0.),
            ...                           (0., 1.), (1., 1.), (2., 1.), (3., 1.),
            ...                           (0., 2.), (1., 2.), (2., 2.), (3., 2.),
            ...                           (4., 1.)))
            >>> vertices *= Numeric.array((dx, dy))
            >>> faces = Numeric.array(((1, 0), (2, 1), (3, 2),
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
            
            >>> externalFaces = Numeric.array((0, 1, 2, 6, 7, 8, 9, 13, 17, 19))
            >>> numerix.allequal(externalFaces, mesh.getExteriorFaces())
            1

            >>> internalFaces = Numeric.array((3, 4, 5, 10, 11, 12, 14, 15, 16, 18))
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
            
            >>> faceAreas = Numeric.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy,
            ...                            Numeric.sqrt(dx**2 + dy**2), dx, Numeric.sqrt(dx**2 + dy**2)))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = Numeric.take(vertices, faces)
            >>> faceCenters = (faceCoords[:,0] + faceCoords[:,1]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = Numeric.array(((0., -1.), (0., -1.), (0., -1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (0., 1.), (0., 1.), (0., 1.),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0),
            ...                              (-1., 0), (1., 0), (1., 0), (1., 0),
            ...                              (dy / Numeric.sqrt(dx**2 + dy**2), -dx / Numeric.sqrt(dx**2 + dy**2)),
            ...                              (0., 1.),
            ...                              (dy / Numeric.sqrt(dx**2 + dy**2), dx / Numeric.sqrt(dx**2 + dy**2))))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, 1, 1, 1), (1, 1, 1, -1), (1, 1, 1, -1),
            ...                                            (-1, 1, 1, 1), (-1, 1, 1, -1), (-1, 1, 1, -1),
            ...                                            (1, 1, -1, 0), (-1, 1, -1, 0)), 0)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = Numeric.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy / 2., dx*dy / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = Numeric.array(((dx/2.,dy/2.), (3.*dx/2.,dy/2.), (5.*dx/2.,dy/2.), 
            ...                              (dx/2.,3.*dy/2.), (3.*dx/2.,3.*dy/2.), (5.*dx/2.,3.*dy/2.),
            ...                              (3.*dx+dx/3.,2.*dy/3.), (3.*dx+dx/3.,4.*dy/3.)))
            >>> numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dy / 2., dy / 2.), (dy / 2., dy / 2.), (dy / 2., dy / 2.),
            ...                                         (dy / 2., -1), (dy / 2., -1), (dy / 2., -1),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
            ...                                         (dx / 2., -1), (dx / 2., dx / 2.), (dx / 2., dx / 2.),
            ...                                         (dx / 2., Numeric.sqrt((dx / 3.)**2 + (dy / 6.)**2)),
            ...                                         (Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1),
            ...                                         (Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2), Numeric.sqrt((dx / 6.)**2 + (dy / 3.)**2)),
            ...                                         (Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2), -1)), -1)
            >>> numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
                                              
            >>> cellDistances = Numeric.array((dy / 2., dy / 2., dy / 2.,
            ...                                dy, dy, dy,
            ...                                dy / 2., dy / 2., dy / 2.,
            ...                                dx / 2., dx, dx,
            ...                                Numeric.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
            ...                                dx / 2., dx, dx,
            ...                                Numeric.sqrt((dx / 3. + dx / 2.)**2 + (dy / 6.)**2),
            ...                                Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2),
            ...                                2. * dy / 3.,
            ...                                Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2)))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[...,0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas[...,Numeric.NewAxis]
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = Numeric.array(((1., 0), (1., 0),(1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (-1., 0), (-1., 0),(-1., 0),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.),
            ...                            (0., -1.), (0., 1.), (0., 1.), (0., 1.),
            ...                            (dx / Numeric.sqrt(dx**2 +dy**2), dy / Numeric.sqrt(dx**2 +dy**2)),
            ...                            (-1, 0.),
            ...                            (-dx / Numeric.sqrt(dx**2 +dy**2), dy / Numeric.sqrt(dx**2 +dy**2))))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = Numeric.array(((0., 0), (0., 0),(0., 0),
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

            >>> d1 = Numeric.sqrt((5. * dx / 6.)**2 + (dy / 6.)**2)
            >>> d2 = Numeric.sqrt((dx / 6.)**2 + (dy / 6.)**2)
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

            >>> interiorCellIDs = Numeric.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = Numeric.array((0, 1, 2, 3, 4, 5, 6, 7))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> nx = dy / Numeric.sqrt(dx**2 + dy**2)
            >>> ny = dx / Numeric.sqrt(dx**2 + dy**2)
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

            >>> vv = Numeric.array(((0, -dx), (dy, 0), (0, dx), (-dy, 0)))
            >>> area = Numeric.sqrt(dx**2 + dy**2)
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
