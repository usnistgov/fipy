#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: Alexander Mont <alexander.mont@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.meshes.cell import Cell
from fipy.meshes.topologies import _MeshTopology
from fipy.meshes.geometries import _MeshGeometry

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serial

class MeshAdditionError(Exception):
    pass
    
class Mesh(object):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serial):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""
         
        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.masked_values(faceVertexIDs, -1)
        self.cellFaceIDs = MA.masked_values(cellFaceIDs, -1)
        self.communicator = communicator

        self.dim = self.vertexCoords.shape[0]

        if not hasattr(self, "numberOfFaces"):
            self.numberOfFaces = self.faceVertexIDs.shape[-1]
        if not hasattr(self, "numberOfCells"):
            self.numberOfCells = self.cellFaceIDs.shape[-1]
        if not hasattr(self, "globalNumberOfCells"):
            self.globalNumberOfCells = self.numberOfCells
        if not hasattr(self, "globalNumberOfFaces"):
            self.globalNumberOfFaces = self.numberOfFaces

        self.faceCellIDs = self._calcFaceCellIDs() 

        self._setTopology()
        self._setGeometry(scaleLength = 1.)

    """Topology methods"""
    def _setTopology(self):
        self._topology = _MeshTopology(self.cellFaceIDs, 
                                       self.faceCellIDs, 
                                       self.numberOfCells,
                                       self._maxFacesPerCell,
                                       self) # `self` only for int/ext face calc

    def _setGeometry(self, scaleLength = 1.):
        self._geometry = _MeshGeometry(self.dim,
                                      self.faceVertexIDs,
                                      self.vertexCoords,
                                      self.faceCellIDs,
                                      self.cellFaceIDs,
                                      self.numberOfCells,
                                      self._maxFacesPerCell,
                                      self.cellToFaceOrientations,
                                      scaleLength)
                                      
    def setScale(self, scaleLength = 1.):
        """
        Sets scale of geometry.

        :Parameters:
          - `scaleLength`: The desired scale length.
        """
        self._geometry.scale = scaleLength

    @property
    def _concatenatedClass(self):
        return Mesh

    """Topology properties"""

    def _setExteriorFaces(self, newExtFaces):
        self._topology.exteriorFaces = newExtFaces

    interiorFaces          = property(lambda s: s._topology.interiorFaces)
    exteriorFaces          = property(lambda s: s._topology.exteriorFaces,
                                      _setExteriorFaces)
    interiorCellIDs        = property(lambda s: s._topology.interiorCellIDs)
    exteriorCellIDs        = property(lambda s: s._topology.exteriorCellIDs)
    cellToFaceOrientations = property(lambda s: s._topology.cellToFaceOrientations)
    adjacentCellIDs        = property(lambda s: s._topology.adjacentCellIDs)
    cellToCellIDs          = property(lambda s: s._topology.cellToCellIDs)
    cellToCellIDsFilled    = property(lambda s: s._topology.cellToCellIDsFilled)

    """geometry properties"""
    faceAreas                 = property(lambda s: s._geometry.scaledFaceAreas)
    faceCenters               = property(lambda s: s._geometry.faceCenters)

    def _setFaceToCellDistances(self, v):
        self._geometry.faceToCellDistances = v

    faceToCellDistances = property(lambda s: s._geometry.faceToCellDistances,
                                   _setFaceToCellDistances)

    def _setCellDistances(self, v):
        self._geometry.cellDistances = v

    cellDistances = property(lambda s: s._geometry.scaledCellDistances,
                             _setCellDistances)

    def _setFaceNormals(self, v):
        self._geometry.faceNormals = v

    faceNormals = property(lambda s: s._geometry.faceNormals,
                           _setFaceNormals)

    cellToFaceDistanceVectors = property(lambda s: s._geometry.cellToFaceDistanceVectors)
    cellDistanceVectors       = property(lambda s: s._geometry.cellDistanceVectors)
    orientedFaceNormals       = property(lambda s: s._geometry.orientedFaceNormals)
    cellVolumes               = property(lambda s: s._geometry.scaledCellVolumes)

    @property
    def cellCenters(self):
        from fipy.variables.cellVariable import CellVariable
        return CellVariable(mesh=self, value=self._geometry.scaledCellCenters,
                            rank=1)

    faceCellToCellNormals     = property(lambda s: s._geometry.faceCellToCellNormals)
    faceTangents1             = property(lambda s: s._geometry.faceTangents1)
    faceTangents2             = property(lambda s: s._geometry.faceTangents2)
    cellToCellDistances       = property(lambda s: s._geometry.scaledCellToCellDistances)
    cellAreas                 = property(lambda s: s._geometry.cellAreas)
    cellNormals               = property(lambda s: s._geometry.cellNormals)

    """scaled geometery properties
    
    These should not exist."""
    scale                     = property(lambda s: s._geometry.scale,
                                         setScale)
    scaledFaceAreas           = property(lambda s: s._geometry.scaledFaceAreas)
    scaledCellVolumes         = property(lambda s: s._geometry.scaledCellVolumes)
    scaledCellCenters         = property(lambda s: s._geometry.scaledCellCenters)
    scaledFaceToCellDistances = property(lambda s: \
                                         s._geometry.scaledFaceToCellDistances)
    scaledCellDistances       = property(lambda s: \
                                         s._geometry.scaledCellDistances)
    scaledCellToCellDistances = property(lambda s: \
                                         s._geometry.scaledCellToCellDistances)
    areaProjections           = property(lambda s: \
                                         s._geometry.areaProjections)
    orientedAreaProjections   = property(lambda s: \
                                         s._geometry.orientedAreaProjections)
    faceToCellDistanceRatio   = property(lambda s: \
                                         s._geometry.faceToCellDistanceRatio)
    faceAspectRatios          = property(lambda s: \
                                         s._geometry.faceAspectRatios)  
        
    def __add__(self, other):
        """
        Either translate a `Mesh` or concatenate two `Mesh` objects.
        
            >>> from fipy.meshes.grid2D import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5]
             [ 0.5  0.5  1.5  1.5]]
             
        If a vector is added to a `Mesh`, a translated `Mesh` is returned
        
            >>> translatedMesh = baseMesh + ((5,), (10,))
            >>> print translatedMesh.getCellCenters()
            [[  5.5   6.5   5.5   6.5]
             [ 10.5  10.5  11.5  11.5]]

             
        If a `Mesh` is added to a `Mesh`, a concatenation of the two 
        `Mesh` objects is returned
        
            >>> addedMesh = baseMesh + (baseMesh + ((2,), (0,)))
            >>> print addedMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5  2.5  3.5  2.5  3.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]
        
        The two `Mesh` objects need not be properly aligned in order to concatenate them
        but the resulting mesh may not have the intended connectivity
        
            >>> from fipy.meshes.mesh import MeshAdditionError
            >>> addedMesh = baseMesh + (baseMesh + ((3,), (0,))) 
            >>> print addedMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5  3.5  4.5  3.5  4.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]

            >>> addedMesh = baseMesh + (baseMesh + ((2,), (2,)))
            >>> print addedMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5  2.5  3.5  2.5  3.5]
             [ 0.5  0.5  1.5  1.5  2.5  2.5  3.5  3.5]]

        No provision is made to avoid or consolidate overlapping `Mesh` objects
        
            >>> addedMesh = baseMesh + (baseMesh + ((1,), (0,)))
            >>> print addedMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5  1.5  2.5  1.5  2.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]
            
        Different `Mesh` classes can be concatenated
         
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx = 1.0, dy = 1.0, nx = 2, ny = 1)
            >>> triMesh = triMesh + ((2,), (0,))
            >>> triAddedMesh = baseMesh + triMesh
            >>> cellCenters = [[0.5, 1.5, 0.5, 1.5, 2.83333333,  3.83333333,
            ...                 2.5, 3.5, 2.16666667, 3.16666667, 2.5, 3.5],
            ...                [0.5, 0.5, 1.5, 1.5, 0.5, 0.5, 0.83333333, 0.83333333, 
            ...                 0.5, 0.5, 0.16666667, 0.16666667]]
            >>> print numerix.allclose(triAddedMesh.getCellCenters(),
            ...                        cellCenters)
            True

        again, their faces need not align, but the mesh may not have 
        the desired connectivity
        
            >>> triMesh = Tri2D(dx = 1.0, dy = 2.0, nx = 2, ny = 1)
            >>> triMesh = triMesh + ((2,), (0,))
            >>> triAddedMesh = baseMesh + triMesh
            >>> cellCenters = [[ 0.5, 1.5, 0.5, 1.5, 2.83333333, 3.83333333,
            ...                  2.5, 3.5, 2.16666667, 3.16666667, 2.5, 3.5],
            ...                [ 0.5, 0.5, 1.5, 1.5, 1., 1.,
            ...                  1.66666667, 1.66666667, 1., 1., 0.33333333, 0.33333333]]
            >>> print numerix.allclose(triAddedMesh.getCellCenters(),
            ...                        cellCenters)
            True

        `Mesh` concatenation is not limited to 2D meshes
        
            >>> from fipy.meshes.grid3D import Grid3D
            >>> threeDBaseMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                         nx = 2, ny = 2, nz = 2)
            >>> threeDSecondMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                           nx = 1, ny = 1, nz = 1)
            >>> threeDAddedMesh = threeDBaseMesh + (threeDSecondMesh + ((2,), (0,), (0,)))
            >>> print threeDAddedMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5  0.5  1.5  0.5  1.5  2.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5  0.5]
             [ 0.5  0.5  0.5  0.5  1.5  1.5  1.5  1.5  0.5]]

        but the different `Mesh` objects must, of course, have the same 
        dimensionality.
        
            >>> InvalidMesh = threeDBaseMesh + baseMesh
            Traceback (most recent call last):
            ...
            MeshAdditionError: Dimensions do not match
        """  
        if(isinstance(other, Mesh)):
            return self._concatenatedClass(**self._getAddedMeshValues(other=other))
        else:
            return self._translate(other)

    __radd__ = __add__
    
    def __mul__(self, factor):
        """
        Dilate a `Mesh` by `factor`.
        
            >>> from fipy.meshes.grid2D import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.getCellCenters()
            [[ 0.5  1.5  0.5  1.5]
             [ 0.5  0.5  1.5  1.5]]

        The `factor` can be a scalar
        
            >>> dilatedMesh = baseMesh * 3
            >>> print dilatedMesh.getCellCenters()
            [[ 1.5  4.5  1.5  4.5]
             [ 1.5  1.5  4.5  4.5]]

        or a vector
        
            >>> dilatedMesh = baseMesh * ((3,), (2,))
            >>> print dilatedMesh.getCellCenters()
            [[ 1.5  4.5  1.5  4.5]
             [ 1.   1.   3.   3. ]]

        
        but the vector must have the same dimensionality as the `Mesh`
        
            >>> dilatedMesh = baseMesh * ((3,), (2,), (1,))
            Traceback (most recent call last):
            ...
            ValueError: shape mismatch: objects cannot be broadcast to a single shape
        """ 
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(vertexCoords=newCoords, 
                       faceVertexIDs=numerix.array(self.faceVertexIDs), 
                       cellFaceIDs=numerix.array(self.cellFaceIDs))
        return newmesh

    __rmul__ = __mul__

    def __repr__(self):
        return "%s()" % self.__class__.__name__
     
    def _connectFaces(self, faces0, faces1):
        """
        
        Merge faces on the same mesh. This is used to create periodic
        meshes. The first list of faces, `faces1`, will be the faces
        that are used to add to the matrix diagonals. The faces in
        `faces2` will not be used. They aren't deleted but their
        adjacent cells are made to point at `faces1`. The list
        `faces2` are not altered, they still remain as members of
        exterior faces.

           >>> from fipy.meshes.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)

           >>> from fipy.tools import parallel
           >>> print parallel.procID != 0 or (mesh.cellFaceIDs == [[0, 1, 2, 3],
           ...                                                           [7, 8, 10, 11],
           ...                                                           [2, 3, 4, 5],
           ...                                                           [6, 7, 9, 10]]).flatten().all()
           True

           >>> mesh._connectFaces(numerix.nonzero(mesh.getFacesLeft()), numerix.nonzero(mesh.getFacesRight()))

           >>> print parallel.procID != 0 or (mesh.cellFaceIDs == [[0, 1, 2, 3],
           ...                                                           [7, 6, 10, 9],
           ...                                                           [2, 3, 4, 5],
           ...                                                           [6, 7, 9, 10]]).flatten().all()
           True

        """

        ## check for errors

        ## check that faces are members of exterior faces
        from fipy.variables.faceVariable import FaceVariable
        faces = FaceVariable(mesh=self, value=False)
        faces[faces0] = True
        faces[faces1] = True
        assert (faces | self.exteriorFaces == self.exteriorFaces).all()

        ## following assert checks number of faces are equal, normals are opposite and areas are the same
        assert numerix.alltrue(numerix.take(self.areaProjections, faces0, axis=1) 
                               == numerix.take(-self.areaProjections, faces1, axis=1))

        ## extract the adjacent cells for both sets of faces
        faceCellIDs0 = self.faceCellIDs[0]
        faceCellIDs1 = self.faceCellIDs[1]
        ## set the new adjacent cells for `faces0`
        MA.put(faceCellIDs1, faces0, MA.take(faceCellIDs0, faces0))
        MA.put(faceCellIDs0, faces0, MA.take(faceCellIDs0, faces1))
        self.faceCellIDs[0] = faceCellIDs0
        self.faceCellIDs[1] = faceCellIDs1
        
        ## extract the face to cell distances for both sets of faces
        faceToCellDistances0 = self.faceToCellDistances[0]
        faceToCellDistances1 = self.faceToCellDistances[1]
        ## set the new faceToCellDistances for `faces0`
        MA.put(faceToCellDistances1, faces0, MA.take(faceToCellDistances0, faces0))
        MA.put(faceToCellDistances0, faces0, MA.take(faceToCellDistances0, faces1))

        """
        Abandon hope, all ye who enter.

        Some very hacky stuff going on here with property assignment. Temporary
        variables are often used because slice and index assignments DO NOT call
        the property-setters, instead they act on the underlying numpy reference
        directly.

        Does Guido know about this?
        """

        connectedFaceToCellDs = self.faceToCellDistances
        connectedFaceToCellDs[0] = faceToCellDistances0
        connectedFaceToCellDs[1] = faceToCellDistances1
        self.faceToCellDistances = connectedFaceToCellDs

        tempCellDist = self.cellDistances
        ## calculate new cell distances and add them to faces0
        numerix.put(tempCellDist, faces0, MA.take(faceToCellDistances0 + faceToCellDistances1, faces0))
        self.cellDistances = tempCellDist

        tempFaceNormals = self.faceNormals
        ## change the direction of the face normals for faces0
        for dim in range(self.getDim()):
            faceNormals = tempFaceNormals[dim].copy()
            numerix.put(faceNormals, faces0, MA.take(faceNormals, faces1))
            tempFaceNormals[dim] = faceNormals

        self.faceNormals = tempFaceNormals

        ## Cells that are adjacent to faces1 are changed to point at faces0
        ## get the cells adjacent to faces1
        faceCellIDs = MA.take(self.faceCellIDs[0], faces1)
        ## get all the adjacent faces for those particular cells
        cellFaceIDs = numerix.take(self.cellFaceIDs, faceCellIDs, axis=1)
        for i in range(cellFaceIDs.shape[0]):
            ## if the faces is a member of faces1 then change the face to point at
            ## faces0
            cellFaceIDs[i] = MA.where(cellFaceIDs[i] == faces1,
                                      faces0,
                                      cellFaceIDs[i])
            ## add those faces back to the main self.cellFaceIDs
            tmp = self.cellFaceIDs[i]
            numerix.put(tmp, faceCellIDs, cellFaceIDs[i])
            self.cellFaceIDs[i] = tmp

        ## calculate new topology
        self._setTopology()

        ## calculate new geometry
        self._geometry.handleFaceConnection()
        
        self.setScale(self.scale['length'])
        
    def _getConcatenableMesh(self):
        return self
        
    def _getAddedMeshValues(self, other, resolution=1e-2):
        """Calculate the parameters to define a concatenation of `other` with `self`
        
        :Parameters:
          - `other`: The :class:`~fipy.meshes.Mesh` to concatenate with `self`
          - `resolution`: How close vertices have to be (relative to the smallest 
            cell-to-cell distance in either mesh) to be considered the same

        :Returns:
          A `dict` with 3 elements: the new mesh vertexCoords, faceVertexIDs, and cellFaceIDs.
        """
        
        selfc = self._getConcatenableMesh()
        other = other._getConcatenableMesh()

        selfNumFaces = selfc.faceVertexIDs.shape[-1]
        selfNumVertices = selfc.vertexCoords.shape[-1]
        otherNumFaces = other.faceVertexIDs.shape[-1]
        otherNumVertices = other.vertexCoords.shape[-1]
        ## check dimensions
        if(selfc.vertexCoords.shape[0] != other.vertexCoords.shape[0]):
            raise MeshAdditionError, "Dimensions do not match"
            
        ## compute vertex correlates

        ## only try to match exterior (X) vertices
        self_Xvertices = numerix.unique(selfc.faceVertexIDs.filled()[...,
            self.exteriorFaces.getValue()].flatten())
        other_Xvertices = numerix.unique(other.faceVertexIDs.filled()[...,
            other.exteriorFaces.getValue()].flatten())

        self_XvertexCoords = selfc.vertexCoords[..., self_Xvertices]
        other_XvertexCoords = other.vertexCoords[..., other_Xvertices]
        
        # lifted from Mesh._getNearestCellID()
        other_vertexCoordMap = numerix.resize(other_XvertexCoords, 
                                              (self_XvertexCoords.shape[-1], 
                                               other_XvertexCoords.shape[0], 
                                               other_XvertexCoords.shape[-1])).swapaxes(0,1)
        tmp = self_XvertexCoords[..., numerix.newaxis] - other_vertexCoordMap
        closest = numerix.argmin(numerix.dot(tmp, tmp), axis=0)
        
        # just because they're closest, doesn't mean they're close
        tmp = self_XvertexCoords[..., closest] - other_XvertexCoords
        distance = numerix.sqrtDot(tmp, tmp)
        # only want vertex pairs that are 100x closer than the smallest 
        # cell-to-cell distance
        close = distance < resolution * min(selfc._getCellToCellDistances().min(), 
                                            other._getCellToCellDistances().min())
        vertexCorrelates = numerix.array((self_Xvertices[closest[close]],
                                          other_Xvertices[close]))
        
        # warn if meshes don't touch, but allow it
        if (selfc._getNumberOfVertices() > 0 
            and other._getNumberOfVertices() > 0 
            and vertexCorrelates.shape[-1] == 0):
            import warnings
            warnings.warn("Vertices are not aligned", UserWarning, stacklevel=4)

        ## compute face correlates

        # ensure that both sets of faceVertexIDs have the same maximum number of (masked) elements
        self_faceVertexIDs = selfc.faceVertexIDs
        other_faceVertexIDs = other.faceVertexIDs

        diff = self_faceVertexIDs.shape[0] - other_faceVertexIDs.shape[0]
        if diff > 0:
            other_faceVertexIDs = numerix.append(other_faceVertexIDs, 
                                                 -1 * numerix.ones((diff,) 
                                                                   + other_faceVertexIDs.shape[1:]),
                                                 axis=0)
            other_faceVertexIDs = MA.masked_values(other_faceVertexIDs, -1)
        elif diff < 0:
            self_faceVertexIDs = numerix.append(self_faceVertexIDs, 
                                                -1 * numerix.ones((-diff,) 
                                                                  + self_faceVertexIDs.shape[1:]),
                                                axis=0)
            self_faceVertexIDs = MA.masked_values(self_faceVertexIDs, -1)

        # want self's Faces for which all faceVertexIDs are in vertexCorrelates
        self_matchingFaces = numerix.in1d(self_faceVertexIDs, 
                                          vertexCorrelates[0]).reshape(self_faceVertexIDs.shape).all(axis=0).nonzero()[0]

        # want other's Faces for which all faceVertexIDs are in vertexCorrelates
        other_matchingFaces = numerix.in1d(other_faceVertexIDs, 
                                           vertexCorrelates[1]).reshape(other_faceVertexIDs.shape).all(axis=0).nonzero()[0]
                                           
        # map other's Vertex IDs to new Vertex IDs, 
        # accounting for overlaps with self's Vertex IDs
        vertex_map = numerix.empty(otherNumVertices, dtype=int)
        verticesToAdd = numerix.delete(numerix.arange(otherNumVertices), vertexCorrelates[1])
        vertex_map[verticesToAdd] = numerix.arange(otherNumVertices - len(vertexCorrelates[1])) + selfNumVertices
        vertex_map[vertexCorrelates[1]] = vertexCorrelates[0]

        # calculate hashes of faceVertexIDs for comparing Faces
        
        if self_matchingFaces.shape[-1] == 0:
            self_faceHash = numerix.empty(self_matchingFaces.shape[:-1] + (0,), dtype="str")
        else:
            # sort each of self's Face's vertexIDs for canonical comparison
            self_faceHash = numerix.sort(self_faceVertexIDs[..., self_matchingFaces], axis=0)
            # then hash the Faces for comparison (NumPy set operations are only for 1D arrays)
            self_faceHash = numerix.apply_along_axis(str, axis=0, arr=self_faceHash)
            
        face_sort = numerix.argsort(self_faceHash)
        self_faceHash = self_faceHash[face_sort]
        self_matchingFaces = self_matchingFaces[face_sort]

        if other_matchingFaces.shape[-1] == 0:
            other_faceHash = numerix.empty(other_matchingFaces.shape[:-1] + (0,), dtype="str")
        else:
            # convert each of other's Face's vertexIDs to new IDs
            other_faceHash = vertex_map[other_faceVertexIDs[..., other_matchingFaces]]
            # sort each of other's Face's vertexIDs for canonical comparison
            other_faceHash = numerix.sort(other_faceHash, axis=0)
            # then hash the Faces for comparison (NumPy set operations are only for 1D arrays)
            other_faceHash = numerix.apply_along_axis(str, axis=0, arr=other_faceHash)

        face_sort = numerix.argsort(other_faceHash)
        other_faceHash = other_faceHash[face_sort]
        other_matchingFaces = other_matchingFaces[face_sort]

        self_matchingFaces = self_matchingFaces[numerix.in1d(self_faceHash, 
                                                             other_faceHash)]
        other_matchingFaces = other_matchingFaces[numerix.in1d(other_faceHash, 
                                                               self_faceHash)]
        
        faceCorrelates = numerix.array((self_matchingFaces,
                                        other_matchingFaces))

        # warn if meshes don't touch, but allow it
        if (selfc._getNumberOfFaces() > 0 
            and other._getNumberOfFaces() > 0 
            and faceCorrelates.shape[-1] == 0):
            import warnings
            warnings.warn("Faces are not aligned", UserWarning, stacklevel=4)

        # map other's Face IDs to new Face IDs, 
        # accounting for overlaps with self's Face IDs
        face_map = numerix.empty(otherNumFaces, dtype=int)
        facesToAdd = numerix.delete(numerix.arange(otherNumFaces), faceCorrelates[1])
        face_map[facesToAdd] = numerix.arange(otherNumFaces - len(faceCorrelates[1])) + selfNumFaces
        face_map[faceCorrelates[1]] = faceCorrelates[0]
        
        other_faceVertexIDs = vertex_map[other.faceVertexIDs[..., facesToAdd]]
        
        # ensure that both sets of cellFaceIDs have the same maximum number of (masked) elements
        self_cellFaceIDs = selfc.cellFaceIDs
        other_cellFaceIDs = face_map[other.cellFaceIDs]
        diff = self_cellFaceIDs.shape[0] - other_cellFaceIDs.shape[0]
        if diff > 0:
            other_cellFaceIDs = numerix.append(other_cellFaceIDs, 
                                               -1 * numerix.ones((diff,) 
                                                                 + other_cellFaceIDs.shape[1:]),
                                               axis=0)
            other_cellFaceIDs = MA.masked_values(other_cellFaceIDs, -1)
        elif diff < 0:
            self_cellFaceIDs = numerix.append(self_cellFaceIDs, 
                                              -1 * numerix.ones((-diff,) 
                                                                + self_cellFaceIDs.shape[1:]),
                                              axis=0)
            self_cellFaceIDs = MA.masked_values(self_cellFaceIDs, -1)

        # concatenate everything and return
        return {
            'vertexCoords': numerix.concatenate((selfc.vertexCoords, 
                                                 other.vertexCoords[..., verticesToAdd]), axis=1), 
            'faceVertexIDs': numerix.concatenate((self_faceVertexIDs, 
                                                  other_faceVertexIDs), axis=1), 
            'cellFaceIDs': MA.concatenate((self_cellFaceIDs, 
                                           other_cellFaceIDs), axis=1)
            }

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    """calc Topology methods"""

    def _getNumberOfFacesPerCell(self):
        cellFaceIDs = self.cellFaceIDs
        if type(cellFaceIDs) is type(MA.array(0)):
            ## bug in count returns float values when there is no mask
            return numerix.array(cellFaceIDs.count(axis=0), 'l')
        else:
            return self._getMaxFacesPerCell() * numerix.ones(cellFaceIDs.shape[-1], 'l')
      
    """
    TODO: Does this really belong in mesh? I don't think so.
    """
    def _calcFaceCellIDs(self):
        array = MA.array(MA.indices(self.cellFaceIDs.shape, 'l')[1], 
                         mask=MA.getmask(self.cellFaceIDs))
        faceCellIDs = MA.zeros((2, self.numberOfFaces), 'l')

        ## Nasty bug: MA.put(arr, ids, values) fills its ids and
        ## values arguments when masked!  This was not the behavior
        ## that was assumed when used below.  It was only working
        ## because the old fill value was 0 and the first element of
        ## the array needed to be 0 since the cell's face was
        ## 0. numerix.put() has been changed to deal with this
        ## properly.

##         MA.put(firstRow, cellFaceIDsFlat[::-1], array[::-1])
##         MA.put(secondRow, cellFaceIDsFlat, array)
        firstRow = faceCellIDs[0]
        secondRow = faceCellIDs[1]
        numerix.put(firstRow, self.cellFaceIDs[::-1,::-1], array[::-1,::-1])
        numerix.put(secondRow, self.cellFaceIDs, array)
        
        mask = ((False,) * self.numberOfFaces, (firstRow == secondRow))
        return MA.sort(MA.array(faceCellIDs, mask = mask),
                                   axis=0)

    def _calcNumPts(self, d, n = None, axis = "x"):
        """
        Calculate the number of cells along the specified axis, based
        on either the specified number or on the number elements in the
        cell  `d` spacings.
        
        Used by the `Grid` meshes.

        This tests a bug that was occuring with PeriodicGrid1D when
        using a numpy float as the argument for the grid spacing.

           >>> from fipy.meshes.periodicGrid1D import PeriodicGrid1D
           >>> PeriodicGrid1D(nx=2, dx=numerix.float32(1.))
           PeriodicGrid1D(dx=1.0, nx=2)

        """

        if type(d) in [type(1), type(1.)] or not hasattr(d, '__len__'):
            n = int(n or 1)
        else:
            n = int(n or len(d))
            if n != len(d) and len(d) != 1:
                raise IndexError, "n%s != len(d%s)" % (axis, axis)
                
        return n

    def _calcVertexCoordinates(self, d, n):
        """
        Calculate the positions of the vertices along an axis, based on the 
        specified `Cell` `d` spacing or list of `d` spacings.
        
        Used by the `Grid` meshes.
        """
        x = numerix.zeros((n + 1), 'd')
        if n > 0:
            x[1:] = d
        return numerix.add.accumulate(x)

    """get Topology methods"""
    
    def getVertexCoords(self):
        """TODO: replace this with a warning."""
        if hasattr(self, 'vertexCoords'):
            return self.vertexCoords
        else:
            return self._createVertices()

    def getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        TODO: replace with a warning.
        """
        return self.exteriorFaces
            
    def getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        TODO: replace with a warning.
        """
        return self.interiorFaces

    def getFaceCellIDs(self):
        return self.faceCellIDs

    def _getMaxFacesPerCell(self):
        return self._maxFacesPerCell

    @property
    def _maxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    def _getExteriorCellIDs(self):
        """ Why do we have this?!? It's only used for testing against itself? """
        return self.exteriorCellIDs

    def _getInteriorCellIDs(self):
        """ Why do we have this?!? It's only used for testing against itself? """
        return self.interiorCellIDs

    def _getCellFaceOrientations(self):
        return self.cellToFaceOrientations

    def getNumberOfCells(self):
        return self.numberOfCells

    def _isOrthogonal(self):
        return False
    
    def _getNumberOfVertices(self):
        if hasattr(self, 'numberOfVertices'):
            return self.numberOfVertices
        else:
            return self.vertexCoords.shape[-1]
        
    def _getAdjacentCellIDs(self):
        return self.adjacentCellIDs

    def getDim(self):
        return self.dim

    def _getGlobalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfCells)

    def _getGlobalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6] for mesh A

            A        B
        ------------------
        | 4 | 5 || 6 | 7 |
        ------------------
        | 0 | 1 || 2 | 3 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfCells)

    def _getLocalNonOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3] for mesh A

            A        B
        ------------------
        | 3 | 4 || 4 | 5 |
        ------------------
        | 0 | 1 || 1 | 2 |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfCells)

    def _getLocalOverlappingCellIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5] for mesh A

            A        B
        ------------------
        | 3 | 4 || 5 |   |
        ------------------
        | 0 | 1 || 2 |   |
        ------------------
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfCells)

    def _getGlobalNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Does not include the IDs of boundary cells.

        E.g., would return [0, 1, 4, 5, 8, 9, 12, 13, 14, 17, 18, 19]
        for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfFaces)

    def _getGlobalOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in the context of the
        global parallel mesh. Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 
        14, 15, 17, 18, 19, 20] for mesh A

            A   ||   B
        --8---9---10--11--
       17   18  19  20   21
        --4---5----6---7--
       12   13  14  15   16
        --0---1----2---3--
                ||
                
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfFaces)

    def _getLocalNonOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Does not include the IDs of boundary cells.
        
        E.g., would return [0, 1, 3, 4, 6, 7, 9, 10, 11, 13, 14, 15]
        for mesh A

            A   ||   B
        --6---7-----7---8--
       13   14 15/14 15   16
        --3---4-----4---5--
        9   10 11/10 11   12
        --0---1-----1---2--
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfFaces)

    def _getLocalOverlappingFaceIDs(self):
        """
        Return the IDs of the local mesh in isolation. 
        Includes the IDs of boundary cells.
        
        E.g., would return [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 
        12, 13, 14, 15, 16] for mesh A

            A   ||   B
        --6---7----8------
       13   14  15  16   |
        --3---4----5------
        9   10  11  12   |
        --0---1----2------
                ||
        
        .. note:: Trivial except for parallel meshes
        """
        return numerix.arange(self.numberOfFaces)

    def getFacesLeft(self):
        """
        Return face on left boundary of Grid1D as list with the
        x-axis running from left to right.

            >>> from fipy import Grid2D, Grid3D
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((21, 25), 
            ...                              numerix.nonzero(mesh.getFacesLeft())[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> print parallel.procID > 0 or numerix.allequal((9, 13), 
            ...                              numerix.nonzero(mesh.getFacesLeft())[0])
            True

        """
        x = self.getFaceCenters()[0]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=x == _madmin(x))

    def getFacesRight(self):
        """
        Return list of faces on right boundary of Grid3D with the
        x-axis running from left to right. 

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((24, 28), 
            ...                              numerix.nonzero(mesh.getFacesRight())[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)    
            >>> print parallel.procID > 0 or numerix.allequal((12, 16), 
            ...                                               numerix.nonzero(mesh.getFacesRight())[0])
            True
            
        """
        x = self.getFaceCenters()[0]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=x == _madmax(x))

    def getFacesBottom(self):
        """
        Return list of faces on bottom boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((12, 13, 14), 
            ...                              numerix.nonzero(mesh.getFacesBottom())[0])
            1
            >>> x, y, z = mesh.getFaceCenters()
            >>> print parallel.procID > 0 or numerix.allequal((12, 13), 
            ...                              numerix.nonzero(mesh.getFacesBottom() & (x < 1))[0])
            1
            
        """
        y = self.getFaceCenters()[1]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=y == _madmin(y))

    getFacesDown = getFacesBottom

    def getFacesTop(self):
        """
        Return list of faces on top boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((18, 19, 20), 
            ...                              numerix.nonzero(mesh.getFacesTop())[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> print parallel.procID > 0 or numerix.allequal((6, 7, 8), 
            ...                              numerix.nonzero(mesh.getFacesTop())[0])
            True
            
        """
        y = self.getFaceCenters()[1]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=y == _madmax(y))

    getFacesUp = getFacesTop

    def getFacesBack(self):
        """
        Return list of faces on back boundary of Grid3D with the
        z-axis running from front to back. 

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((6, 7, 8, 9, 10, 11), 
            ...                              numerix.nonzero(mesh.getFacesBack())[0])
            True

        """
        z = self.getFaceCenters()[2] 
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=z == _madmax(z))

    def getFacesFront(self):
        """
        Return list of faces on front boundary of Grid3D with the
        z-axis running from front to back. 

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((0, 1, 2, 3, 4, 5), 
            ...                              numerix.nonzero(mesh.getFacesFront())[0])
            True

        """
        z = self.getFaceCenters()[2]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=z == _madmin(z))
    
    def _getNumberOfFaces(self):
        return self.numberOfFaces

    def _getCellToCellIDs(self):
        return self.cellToCellIDs

    def _getCellToCellIDsFilled(self):
        return self.cellToCellIDsFilled
     
    """get geometry methods"""

    def _getFaceAreas(self):
        return self.faceAreas

    def _getFaceNormals(self):
        return self.faceNormals

    def _getFaceCellToCellNormals(self):
        return self.faceCellToCellNormals
        
    def getCellVolumes(self):
        return self.cellVolumes

    """
    This shit has GOT TO GO.
    """
    def _getCellCenters(self):
        return self.scaledCellCenters
        
    def getCellCenters(self):
        return self.cellCenters

    def _getFaceToCellDistances(self):
        return self.faceToCellDistances

    def _getCellDistances(self):
        return self.cellDistances

    def _getFaceToCellDistanceRatio(self):
        return self.faceToCellDistanceRatio

    def _getOrientedAreaProjections(self):
        return self.orientedAreaProjections

    def _getAreaProjections(self):
        return self.areaProjections

    def _getOrientedFaceNormals(self):
        return self.orientedFaceNormals

    def _getFaceTangents1(self):
        return self.faceTangents1

    def _getFaceTangents2(self):
        return self.faceTangents2
        
    def _getFaceAspectRatios(self):
        return self.faceAspectRatios
    
    def _getCellToCellDistances(self):
        return self.cellToCellDistances

    def _getCellNormals(self):
        return self.cellNormals

    def _getCellAreas(self):
        return self.cellAreas

    def _getCellAreaProjections(self):
        return self.cellNormals * self.cellAreas
         
    def getFaceCenters(self):
        return self.faceCenters

    def _getOrderedCellVertexIDs(self):
        return self._getCellVertexIDs()

    def _getCellDistanceNormals(self):
        return self.getCellDistanceVectors() / self.getCellDistances()
        
    def _getCellVertexIDs(self):

        ## Get all the vertices from all the faces for each cell
        cellFaceVertices = numerix.take(self.faceVertexIDs, self.cellFaceIDs, axis=1)

        ## get a sorted list of vertices for each cell 
        cellVertexIDs = numerix.reshape(cellFaceVertices, (-1, self.numberOfCells))
        cellVertexIDs = MA.sort(cellVertexIDs, axis=0, fill_value=-1)

        cellVertexIDs = MA.sort(MA.concatenate((cellVertexIDs[-1, numerix.newaxis], 
                                                MA.masked_where(cellVertexIDs[:-1] 
                                                                == cellVertexIDs[1:], 
                                                                cellVertexIDs[:-1]))), 
                                axis=0, fill_value=-1)
        
        ## resize the array to remove extra masked values
        if cellVertexIDs.shape[-1] == 0:
            length = 0
        else:
            length = min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis=0))
        return cellVertexIDs[length:][::-1]

    """
    Below is an ordered version of _getCellVertexIDs()
    It works for the test case in this file (other than the ordering, obviously)
    I've left it in as it may be useful when we need ordered vertices for cells
  
    def _getOrderedCellVertexIDs(self):

        ## Get all the vertices from all the faces for each cell
        from fipy.tools.numerix import take
        cellFaceVertices = take(self.faceVertexIDs, self.cellFaceIDs)

        ## get a sorted list of vertices for each cell
        NCells = self.numberOfCells
        cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
        newmask = MA.getmaskarray(cellVertexIDs).copy()

        for i in range(len(cellVertexIDs[0]) - 1):
            for j in range(len(cellVertexIDs[0]))[i + 1:]:

                newmask[:,j] = MA.where(newmask[:,j],
                                        newmask[:,j],
                                        MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,j],
                                                 1,
                                                 newmask[:,j]))


        cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)

        for i in range(len(cellVertexIDs[0]) - 1):
            j = i + 1
            while j < len(cellVertexIDs[0]):
                tmp = cellVertexIDs[:]
                tmp[:, i] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
                                              cellVertexIDs[:, i],
                                              cellVertexIDs[:, j]),
                                     cellVertexIDs[:, i])
                                                      
                tmp[:, j] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
                                              cellVertexIDs[:,j],
                                              cellVertexIDs[:,i]),
                                     cellVertexIDs[:, j])

                cellVertexIDs = tmp[:]

                j += 1


        length = len(cellVertexIDs[0]) - min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
        return cellVertexIDs[:, :length]
    """


    """scaling"""
    def _calcScaleArea(self):
        raise NotImplementedError

    def _calcScaleVolume(self):
        raise NotImplementedError
     
    def _getPointToCellDistances(self, point):
        tmp = self.getCellCenters() - PhysicalField(point)
        from fipy.tools import numerix
        return numerix.sqrtDot(tmp, tmp)

    def getNearestCell(self, point):
        return self._getCellsByID([self._getNearestCellID(point)])[0]

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m0 = Grid2D(dx=(.1, 1., 10.), dy=(.1, 1., 10.))
           >>> m1 = Grid2D(nx=2, ny=2, dx=5., dy=5.)
           >>> print m0._getNearestCellID(m1.getCellCenters().getGlobalValue())
           [4 5 7 8]
           
        """
        if self.globalNumberOfCells == 0:
            return numerix.arange(0)
            
        points = numerix.resize(points, (self.globalNumberOfCells, len(points), len(points[0]))).swapaxes(0,1)

        centers = self.getCellCenters().getGlobalValue()[...,numerix.newaxis]
        try:
            tmp = centers - points
        except TypeError:
            tmp = centers - PhysicalField(points)

        return numerix.argmin(numerix.dot(tmp, tmp, axis = 0), axis=0)
     

    """pickling"""

    def __getstate__(self):
        dict = {
            'vertexCoords' : self.vertexCoords *  self.scale['length'],            
            'faceVertexIDs' : self.faceVertexIDs,
            'cellFaceIDs' : self.cellFaceIDs }
        return dict

    def __setstate__(self, dict):
        Mesh.__init__(self, **dict)
##        self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
     
    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected.
        
            >>> dx = 2.
            >>> dy = 1.23456
            >>> dz = 1.e-1
            
            >>> vertices = numerix.array(((0., 1., 1., 0., 0., 1., 1., 0., 2., 2.),
            ...                           (0., 0., 1., 1., 0., 0., 1., 1., 0., 0.),
            ...                           (0., 0., 0., 0., 1., 1., 1., 1., 0., 1.)))
            >>> vertices *= numerix.array([[dx], [dy], [dz]])
            >>> faces = MA.masked_values(((0, 7, 3, 5, 1, 3, 1, 9, 8, 8),
            ...                           (1, 6, 7, 6, 0, 2, 8, 5, 1, 9),
            ...                           (2, 5, 4, 2, 4, 6, 2, 6, 5, 6),
            ...                           (3, 4, 0, 1, 5, 7, -1, -1, 9, 2)), -1)
            >>> cells = MA.masked_values(((0, 3),
            ...                           (1, 6), 
            ...                           (2, 7),
            ...                           (3, 8),
            ...                           (4, 9),
            ...                           (5, -1)), -1)

            >>> mesh = Mesh(vertexCoords=vertices, faceVertexIDs=faces, cellFaceIDs=cells)

            >>> externalFaces = numerix.array((0, 1, 2, 4, 5, 6, 7, 8, 9))
            >>> print numerix.allequal(externalFaces, 
            ...                        numerix.nonzero(mesh.exteriorFaces))
            1

            >>> internalFaces = numerix.array((3,))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.interiorFaces))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  0,  0, 0,  0,  0,  1,  1,  1,  1),
            ...                                 (-1, -1, -1, 1, -1, -1, -1, -1, -1, -1)), -1)
            >>> numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> numerix.allclose(faceAreas, mesh._getFaceAreas(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> numerix.allclose(faceCenters, mesh.getFaceCenters(), atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> numerix.allclose(faceNormals, mesh._getFaceNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> numerix.allclose(cellVolumes, mesh.getCellVolumes(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., dx+dx/3.),
            ...                              (dy/2.,    dy/3.),
            ...                              (dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.getCellCenters(), atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> d1 = numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = numerix.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., dx / 2., dy / 2., dy / 2., dz / 2., dz / 2., d2, d3),
            ...                                         (     -1,      -1,      -1,      d1,      -1,      -1,      -1,      -1, -1, -1)), -1)
            >>> print numerix.allclose(faceToCellDistances, mesh._getFaceToCellDistances(), atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> numerix.allclose(cellDistances, mesh._getCellDistances(), atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._getFaceToCellDistanceRatio(), atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._getAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents1, mesh._getFaceTangents1(), atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents2, mesh._getFaceTangents2(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._getCellToCellDistances(), atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
            1

            >>> cellNormals = MA.masked_values((((0, -1),
            ...                                  (0, 0),
            ...                                  (-1, 0),
            ...                                  (1, 0),
            ...                                  (0, dy / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (0, -1000)),
            ...                                 ((0, 0),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1),
            ...                                  (-1, dx / numerix.sqrt(dx**2 + dy**2)),
            ...                                  (1, -1000)),
            ...                                 ((-1, 0),
            ...                                  (1, -1),
            ...                                  (0, 1),
            ...                                  (0, 0),
            ...                                  (0, 0),
            ...                                  (0, -1000))), -1000)
            >>> numerix.allclose(cellNormals, mesh._getCellNormals(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellAreaProjections = MA.masked_values((((0, -dy * dz),
            ...                                          (0, 0),
            ...                                          (-dy * dz, 0),
            ...                                          ( dy * dz, 0),
            ...                                          (0, dy * dz),
            ...                                          (0, -1000)),
            ...                                         ((0, 0),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -dx * dz),
            ...                                          (-dx * dz, dx * dz),
            ...                                          ( dx * dz, -1000)),
            ...                                         ((-dx * dy, 0),
            ...                                          ( dx * dy, -dx * dy / 2.),
            ...                                          (0, dx * dy / 2.),
            ...                                          (0, 0),
            ...                                          (0, 0),
            ...                                          (0, -1000))), -1000)
            >>> numerix.allclose(cellAreaProjections, mesh._getCellAreaProjections(), atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), (9, 8, 6, 5, 2, 1, -1000, -1000)), -1000)
            >>> cellVertexIDs = MA.masked_values(((7, 9),
            ...                                   (6, 8),
            ...                                   (5, 6),
            ...                                   (4, 5),
            ...                                   (3, 2),
            ...                                   (2, 1),
            ...                                   (1, -1000),
            ...                                   (0, -1000)), -1000)
            >>> numerix.allclose(cellVertexIDs, mesh._getCellVertexIDs())
            1


            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.getCellCenters(), unpickledMesh.getCellCenters())
            True

            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 10
            >>> ny = 2
            >>> from fipy.meshes.grid2D import Grid2D
            >>> gridMesh = Grid2D(dx, dy, nx, ny)
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx, dy, nx, 1) + [[dx*nx], [0]]
            >>> bigMesh = gridMesh + triMesh
            >>> x, y = bigMesh.getCellCenters()
            >>> from fipy.variables.cellVariable import CellVariable
            >>> volumes = CellVariable(mesh=bigMesh, value=1.)
            >>> volumes[x > dx * nx] = 0.25
            >>> print numerix.allclose(bigMesh.getCellVolumes(), volumes)
            True
            
            Following test was added due to a bug in adding UniformGrids.

            >>> from fipy.meshes.uniformGrid1D import UniformGrid1D
            >>> a = UniformGrid1D(nx=10) + (10,)
            >>> print a.getCellCenters()
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            >>> b = 10 + UniformGrid1D(nx=10)
            >>> print b.getCellCenters()
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            
            >>> from fipy.tools import parallel
            >>> if parallel.Nproc == 1:
            ...     c =  UniformGrid1D(nx=10) + (UniformGrid1D(nx=10) + 10)
            >>> print (parallel.Nproc > 1 
            ...        or numerix.allclose(c.getCellCenters()[0],
            ...                            [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            ...                            12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5]))
            True

        """

    def _getVTKCellType(self):
        from enthought.tvtk.api import tvtk
        return tvtk.ConvexPointSet().cell_type
                
    def getVTKCellDataSet(self):
        """Returns a TVTK `DataSet` representing the cells of this mesh
        """
        cvi = self._getOrderedCellVertexIDs().swapaxes(0,1)
        from fipy.tools import numerix
        if type(cvi) is numerix.ma.masked_array:
            counts = cvi.count(axis=1)[:,None]
            cells = numerix.ma.concatenate((counts,cvi),axis=1).compressed()
        else:
            counts = numerix.array([cvi.shape[1]]*cvi.shape[0])[:,None]
            cells = numerix.concatenate((counts,cvi),axis=1).flatten()
        
        from enthought.tvtk.api import tvtk
        num = counts.shape[0]

        cps_type = self._getVTKCellType()
        cell_types = numerix.array([cps_type]*num)
        cell_array = tvtk.CellArray()
        cell_array.set_cells(num, cells)

        points = self.vertexCoords
        points = self._toVTK3D(points)
        ug = tvtk.UnstructuredGrid(points=points)
        
        offset = numerix.cumsum(counts[:,0]+1)
        if len(offset) > 0:
            offset -= offset[0]
        ug.set_cells(cell_types, offset, cell_array)

        return ug

    def getVTKFaceDataSet(self):
        """Returns a TVTK `DataSet` representing the face centers of this mesh
        """
        from enthought.tvtk.api import tvtk
        
        points = self.getFaceCenters()
        points = self._toVTK3D(points)
        ug = tvtk.UnstructuredGrid(points=points)
        
        num = len(points)
        counts = numerix.array([1] * num)[..., numerix.newaxis]
        cells = numerix.arange(self._getNumberOfFaces())[..., numerix.newaxis]
        cells = numerix.concatenate((counts, cells), axis=1)
        cell_types = numerix.array([tvtk.Vertex().cell_type]*num)
        cell_array = tvtk.CellArray()
        cell_array.set_cells(num, cells)

        counts = numerix.array([1] * num)
        offset = numerix.cumsum(counts+1)
        if len(offset) > 0:
            offset -= offset[0]
        ug.set_cells(cell_types, offset, cell_array)

        return ug

    def _toVTK3D(self, arr, rank=1):
        if arr.dtype.name is 'bool':
            # VTK can't do bool, and the exception isn't properly
            # thrown back to the user
            arr = arr.astype('int')
        if rank == 0:
            return arr
        else:
            arr = numerix.concatenate((arr, 
                                       numerix.zeros((3 - self.dim,) 
                                                     + arr.shape[1:])))
            return arr.swapaxes(-2, -1)

def _madmin(x):
    if len(x) == 0:
        return 0
    else:
        return min(x)
        
def _madmax(x):
    if len(x) == 0:
        return 0
    else:
        return max(x)
 
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
