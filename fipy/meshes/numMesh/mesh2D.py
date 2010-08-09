#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh2D.py"
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
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
        tangent = faceVertexCoords[:,1] - faceVertexCoords[:,0]
        self.faceAreas = numerix.sqrtDot(tangent, tangent)

    def _calcFaceNormals(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        self.faceNormals = t1.copy()
##         mag = numerix.sqrtDot(t1, t1)
        mag = numerix.sqrt(t1[1]**2 + t1[0]**2)
        self.faceNormals[0] = -t1[1] / mag
        self.faceNormals[1] = t1[0] / mag
        
        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.cellDistanceVectors) < 0)
        self.faceNormals = self.faceNormals * orientation


    def _calcFaceTangents(self):
        tmp = numerix.array((-self.faceNormals[1], self.faceNormals[0]))
        ## copy required to get internal memory ordering correct for inlining.
        tmp = tmp.copy()
        mag = numerix.sqrtDot(tmp, tmp)
        self.faceTangents1 = tmp / mag
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

    @property
    def _concatenatedClass(self):
        return Mesh2D
        
    def _getOrderedCellVertexIDs(self):
        from fipy.tools.numerix import take
        NFac = self._getMaxFacesPerCell()

        # numpy 1.1's MA.take doesn't like FlatIter. Call ravel() instead.
        cellVertexIDs0 = take(self._getFaceVertexIDs()[0], self._getCellFaceIDs().ravel())
        cellVertexIDs1 = take(self._getFaceVertexIDs()[1], self._getCellFaceIDs().ravel())
        cellVertexIDs = MA.where(self.cellToFaceOrientations.ravel() > 0,
                             cellVertexIDs0, cellVertexIDs1)

        cellVertexIDs = numerix.reshape(cellVertexIDs, (NFac, -1))
        return cellVertexIDs
    
    def _getNonOrthogonality(self):
        
        exteriorFaceArray = numerix.zeros((self.faceCellIDs.shape[1],))
        numerix.put(exteriorFaceArray, numerix.nonzero(self.getExteriorFaces()), 1)
        unmaskedFaceCellIDs = MA.filled(self.faceCellIDs, 0) ## what we put in for the "fill" doesn't matter because only exterior faces have anything masked, and exterior faces have their displacement vectors set to zero.
        ## if it's an exterior face, make the "displacement vector" equal to zero so the cross product will be zero.
    
        faceDisplacementVectors = numerix.where(numerix.array(zip(exteriorFaceArray, exteriorFaceArray)), 0.0, numerix.take(self._getCellCenters().swapaxes(0,1), unmaskedFaceCellIDs[1, :]) - numerix.take(self._getCellCenters().swapaxes(0,1), unmaskedFaceCellIDs[0, :])).swapaxes(0,1)
        faceCrossProducts = (faceDisplacementVectors[0, :] * self.faceNormals[1, :]) - (faceDisplacementVectors[1, :] * self.faceNormals[0, :])
        faceDisplacementVectorLengths = numerix.maximum(((faceDisplacementVectors[0, :] ** 2) + (faceDisplacementVectors[1, :] ** 2)) ** 0.5, 1.e-100)
        faceWeightedNonOrthogonalities = abs(faceCrossProducts / faceDisplacementVectorLengths) * self.faceAreas
        cellFaceWeightedNonOrthogonalities = numerix.take(faceWeightedNonOrthogonalities, self.cellFaceIDs)
        cellFaceAreas = numerix.take(self.faceAreas, self.cellFaceIDs)
        cellTotalWeightedValues = numerix.add.reduce(cellFaceWeightedNonOrthogonalities, axis = 0)  
        cellTotalFaceAreas = numerix.add.reduce(cellFaceAreas, axis = 0)
  
        return (cellTotalWeightedValues / cellTotalFaceAreas)

    def extrude(self, extrudeFunc=lambda x: x + numerix.array((0, 0, 1))[:,numerix.newaxis] , layers=1):
        """
        This function returns a new 3D mesh. The 2D mesh is extruded
        using the extrudeFunc and the number of layers.

        :Parameters:
          - `extrudeFunc`: function that takes the vertex coordinates and returns the displaced values
          - `layers`: the number of layers in the extruded mesh (number of times extrudeFunc will be called)

        >>> from fipy.meshes.grid2D import Grid2D
        >>> print Grid2D(nx=2,ny=2).extrude(layers=2).getCellCenters()
        [[ 0.5  1.5  0.5  1.5  0.5  1.5  0.5  1.5]
         [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]
         [ 0.5  0.5  0.5  0.5  1.5  1.5  1.5  1.5]]

        >>> from fipy.meshes.tri2D import Tri2D
        >>> print Tri2D().extrude(layers=2).getCellCenters()
        [[ 0.83333333  0.5         0.16666667  0.5         0.83333333  0.5
           0.16666667  0.5       ]
         [ 0.5         0.83333333  0.5         0.16666667  0.5         0.83333333
           0.5         0.16666667]
         [ 0.5         0.5         0.5         0.5         1.5         1.5         1.5
           1.5       ]]
        """

        return self._extrude(self, extrudeFunc, layers)

    def _extrude(self, mesh, extrudeFunc, layers):
        ## should extrude cnahe self rather than creating a new mesh?

        ## the following allows the 2D mesh to be in 3D space, this can be the case for a
        ## GmshImporter2DIn3DSpace which would then be extruded.
        oldVertices = mesh.getVertexCoords()
        if oldVertices.shape[0] == 2:
            oldVertices = numerix.resize(oldVertices, (3, len(oldVertices[0])))
            oldVertices[2] = 0

        NCells = mesh.getNumberOfCells()
        NFac = mesh._getNumberOfFaces()
        NFacPerCell =  mesh._getMaxFacesPerCell()

        ## set up the initial data arrays
        new_shape = (max(NFacPerCell, 4), (1 + layers)*NCells + layers*NFac)
        faces = numerix.MA.masked_values(-numerix.ones(new_shape), value = -1)
        orderedVertices = mesh._getOrderedCellVertexIDs()
        faces[:NFacPerCell, :NCells] = orderedVertices
        vertices = oldVertices
        vert0 = mesh._getFaceVertexIDs()
        faceCount = NCells
        
        for layer in range(layers):

            ## need this later
            initialFaceCount = faceCount

            ## build the vertices
            newVertices = extrudeFunc(oldVertices)
            vertices = numerix.concatenate((vertices, newVertices), axis=1)

            ## build the faces along the layers
            faces[:NFacPerCell, faceCount: faceCount + NCells] = orderedVertices + len(oldVertices[0]) * (layer + 1)
            try:
                # numpy 1.1 doesn't copy right side before assigning slice
                # See: http://www.mail-archive.com/numpy-discussion@scipy.org/msg09843.html
                faces[:NFacPerCell, faceCount: faceCount + NCells] = faces[:NFacPerCell, faceCount: faceCount + NCells][::-1,:].copy()
            except:
                faces[:NFacPerCell, faceCount: faceCount + NCells] = faces[:NFacPerCell, faceCount: faceCount + NCells][::-1,:]

            faceCount = faceCount + NCells

            vert1 = (vert0 + len(oldVertices[0]))[::-1,:]

            ## build the faces between the layers
            faces[:4, faceCount: faceCount + NFac] = numerix.concatenate((vert0, vert1), axis = 0)[::-1,:]

            vert0 = vert0 + len(oldVertices[0])

            NCells = mesh.getNumberOfCells()

            
            ## build the cells, the first layer has slightly different ordering
            if layer == 0:
                c0 =  numerix.reshape(numerix.arange(NCells), (1, NCells))
                cells = numerix.concatenate((c0, c0 + NCells, mesh._getCellFaceIDs() + 2 * NCells), axis = 0)
            else:
                newCells = numerix.concatenate((c0, c0 + initialFaceCount, mesh._getCellFaceIDs() + faceCount), axis=0)
                newCells[0] = cells[1,-NCells:]
                cells = numerix.concatenate((cells, newCells), axis=1)

            ## keep a count of things for the next layer
            faceCount = faceCount + NFac
            oldVertices = newVertices

        ## return a new mesh, extrude could just as easily act on self
        return Mesh(vertices, faces, cells, communicator=mesh.communicator)

    def _getVTKCellType(self):
        from enthought.tvtk.api import tvtk
        return tvtk.Polygon().cell_type
        
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 0.5
            >>> dy = 2.
            >>> nx = 3
            >>> ny = 2
            
            >>> vertices = numerix.array(((0., 1., 2., 3., 
            ...                            0., 1., 2., 3., 
            ...                            0., 1., 2., 3., 4.),
            ...                           (0., 0., 0., 0., 
            ...                            1., 1., 1., 1., 
            ...                            2., 2., 2., 2., 1.)))
            >>> vertices *= numerix.array(((dx,), (dy,)))
            >>> faces = numerix.array(((1, 2, 3, 4, 5, 6, 8,  9, 10, 
            ...                         0, 5, 6, 7, 4, 9, 10, 11, 12,  7, 11),
            ...                        (0, 1, 2, 5, 6, 7, 9, 10, 11, 
            ...                         4, 1, 2, 3, 8, 5,  6,  7,  3, 12, 12)))
            >>> cells = MA.masked_values((( 0,  1,  2,  3,  4,  5, 17, 18),
            ...                           (10, 11, 12, 14, 15, 16, 18, 19),
            ...                           ( 3,  4,  5,  6,  7,  8, 12, 16),
            ...                           ( 9, 10, 11, 13, 14, 15, -1, -1)), -1)
                
            >>> mesh = Mesh2D(vertexCoords = vertices, faceVertexIDs = faces, cellFaceIDs = cells)
            
            >>> externalFaces = numerix.array((0, 1, 2, 6, 7, 8, 9, 13, 17, 19))
            >>> print numerix.allequal(externalFaces, 
            ...                        numerix.nonzero(mesh.getExteriorFaces()))
            1

            >>> internalFaces = numerix.array((3, 4, 5, 10, 11, 12, 14, 15, 16, 18))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.getInteriorFaces()))
            1

            >>> faceCellIds = MA.masked_values((( 0,  1,  2, 0,  1, 2,  3,  4,  5,  
            ...                                   0,  0,  1, 2,  3, 3,  4,  5,  6, 6,  7),
            ...                                 (-1, -1, -1, 3,  4, 5, -1, -1, -1, 
            ...                                  -1,  1,  2, 6, -1, 4,  5,  7, -1, 7, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> faceAreas = numerix.array((dx, dx, dx, dx, dx, dx, dx, dx, dx,
            ...                            dy, dy, dy, dy, dy, dy, dy, dy,
            ...                            numerix.sqrt(dx**2 + dy**2), dx, numerix.sqrt(dx**2 + dy**2)))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, faces, axis=1)
            >>> faceCenters = (faceCoords[...,0,:] + faceCoords[...,1,:]) / 2.
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array(((0., 0., 0., 0., 0., 0., 0., 0., 0., 
            ...                               -1., 1., 1., 1., -1., 1., 1., 1., 
            ...                               dy / numerix.sqrt(dx**2 + dy**2), 
            ...                               0., 
            ...                               dy / numerix.sqrt(dx**2 + dy**2)),
            ...                              (-1., -1., -1., 1., 1., 1., 1., 1., 1., 
            ...                               0., 0., 0., 0., 0., 0., .0, 0., 
            ...                               -dx / numerix.sqrt(dx**2 + dy**2), 
            ...                               1., 
            ...                               dx / numerix.sqrt(dx**2 + dy**2))))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1,  1,  1, -1, -1, -1,  1, -1),
            ...                                            (1,  1,  1,  1,  1,  1,  1,  1),
            ...                                            (1,  1,  1,  1,  1,  1, -1, -1),
            ...                                            (1, -1, -1,  1, -1, -1,  0,  0)), 0)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy, dx*dy / 2., dx*dy / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., 3.*dx/2., 5.*dx/2., dx/2., 3.*dx/2., 5.*dx/2., 3.*dx+dx/3., 3.*dx+dx/3.),
            ...                              (dy/2., dy/2., dy/2., 3.*dy/2., 3.*dy/2., 3.*dy/2., 2.*dy/3., 4.*dy/3.)))
            >>> print numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> faceToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., 
            ...                                          dy / 2., dy / 2., dy / 2., 
            ...                                          dy / 2., dy / 2., dy / 2., 
            ...                                          dx / 2., dx / 2., dx / 2., dx / 2., 
            ...                                          dx / 2., dx / 2., dx / 2., dx / 2., 
            ...                                          numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2), 
            ...                                          numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2), 
            ...                                          numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)),
            ...                                         (-1, -1, -1, 
            ...                                          dy / 2., dy / 2., dy / 2., 
            ...                                          -1, -1, -1, 
            ...                                          -1, dx / 2., dx / 2., numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2), 
            ...                                          -1, dx / 2., dx / 2., numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2), 
            ...                                          -1, 
            ...                                          numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2), 
            ...                                          -1)), -1)
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
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents1 = numerix.array(((1., 1., 1., -1., -1., -1., -1., -1., -1., 
            ...                             0., 0., 0., 0., 0., 0., 0., 0., 
            ...                             dx / numerix.sqrt(dx**2 +dy**2), 
            ...                             -1., 
            ...                             -dx / numerix.sqrt(dx**2 +dy**2)),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., 
            ...                             -1., 1., 1., 1., -1., 1., 1., 1., 
            ...                             dy / numerix.sqrt(dx**2 +dy**2), 
            ...                             0., 
            ...                             dy / numerix.sqrt(dx**2 +dy**2))))
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tangents2 = numerix.array(((0., 0., 0., 0., -0., -0., -0., -0., -0., 
            ...                             -0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.),
            ...                            (0., 0., 0., 0., 0., 0., 0., 0., 0., 
            ...                             -0., 0., 0., 0., -0., 0., 0., 0., 0., 0., 0.)))
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1, -1, -1,  0,  1,  2, -1,  6),
            ...                                   ( 1,  2,  6,  4,  5,  7,  7, -1),
            ...                                   ( 3,  4,  5, -1, -1, -1,  2,  5),
            ...                                   (-1,  0,  1, -1,  3,  4, -1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> d1 = numerix.sqrt((5. * dx / 6.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> cellToCellDistances = MA.masked_values(((dy / 2., dy / 2., dy / 2., dy, dy, dy, d2, 2. * dy / 3.),
            ...                                         (dx, dx, d1, dx, dx, d1, 2. * dy / 3., d2),
            ...                                         (dy, dy, dy, dy / 2., dy / 2., dy / 2., d1, d1),
            ...                                         (dx / 2., dx, dx, dx / 2., dx, dx, -1, -1)), -1)
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
            >>> cellNormals = MA.masked_values(((( 0,  0,  0,  0,  0,  0,    nx,     0),
            ...                                  ( 1,  1,  1,  1,  1,  1,     0,    nx),
            ...                                  ( 0,  0,  0,  0,  0,  0,    -1,    -1),
            ...                                  (-1, -1, -1, -1, -1, -1, -1000, -1000)),
            ...                                 ((-1, -1, -1, -1, -1, -1,   -ny,    -1),
            ...                                  ( 0,  0,  0,  0,  0,  0,     1,    ny),
            ...                                  ( 1,  1,  1,  1,  1,  1,     0,     0),
            ...                                  ( 0,  0,  0,  0,  0,  0, -1000, -1000))), -1000)
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> area = numerix.sqrt(dx**2 + dy**2)
            >>> cellAreaProjections = MA.masked_values((((  0,  0,  0,  0,  0,  0,  nx * area,         0),
            ...                                          ( dy, dy, dy, dy, dy, dy,          0, nx * area),
            ...                                          (  0,  0,  0,  0,  0,  0,        -dy,       -dy),
            ...                                          (-dy,-dy,-dy,-dy,-dy,-dy,      -1000,     -1000)),
            ...                                         ((-dx,-dx,-dx,-dx,-dx,-dx, -ny * area,       -dx),
            ...                                          (  0,  0,  0,  0,  0,  0,         dx, ny * area),
            ...                                          ( dx, dx, dx, dx, dx, dx,          0,         0),
            ...                                          (  0,  0,  0,  0,  0,  0,      -1000,     -1000))), -1000)
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_array(((5, 6, 7, 9, 10, 11,    12,    12),
            ...                                  (4, 5, 6, 8,  9, 10,     7,    11),
            ...                                  (1, 2, 3, 5,  6,  7,     3,     7),
            ...                                  (0, 1, 2, 4,  5,  6, -1000, -1000)), -1000)

            >>> numerix.allclose(mesh._getCellVertexIDs(), cellVertexIDs)
            1
            

            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            True

            

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
