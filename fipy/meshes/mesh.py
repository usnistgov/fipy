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

from fipy.meshes.topologies import _MeshTopology
from fipy.meshes.geometries import _MeshGeometry

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serial
from fipy.tools.decorators import getsetDeprecated

class MeshAdditionError(Exception):
    pass
    
class Mesh(object):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serial):
        """faceVertexIds and cellFacesIds must be padded with minus ones."""

        from fipy.variables.variable import Variable
        from fipy.variables.vertexVariable import _VertexVariable
        from fipy.variables.faceVariable import FaceVariable
        from fipy.variables.cellVariable import CellVariable
        
        if isinstance(vertexCoords, Variable):
            vertexCoords = vertexCoords._copyValue()
        if isinstance(faceVertexIDs, Variable):
            faceVertexIDs = faceVertexIDs._copyValue()
        if isinstance(cellFaceIDs, Variable):
            cellFaceIDs = cellFaceIDs._copyValue()
            
        self.vertexCoords = _VertexVariable(mesh=self, value=vertexCoords, 
                                            _bootstrap=True)
        self.faceVertexIDs = FaceVariable(mesh=self, 
                                          elementshape=faceVertexIDs.shape[:-1],
                                          value=MA.masked_values(faceVertexIDs, -1), 
                                          _bootstrap=True)
        self.cellFaceIDs = CellVariable(mesh=self, 
                                        elementshape=cellFaceIDs.shape[:-1], 
                                        value=MA.masked_values(cellFaceIDs, -1), 
                                        _bootstrap=True)
        self.communicator = communicator

        self.dim = self.vertexCoords.shape[0]

        if not hasattr(self, "numberOfVertices"):
            self.numberOfVertices = self.vertexCoords.shape[-1]
        if not hasattr(self, "numberOfFaces"):
            self.numberOfFaces = self.faceVertexIDs.shape[-1]
        if not hasattr(self, "numberOfCells"):
            self.numberOfCells = self.cellFaceIDs.shape[-1]
        if not hasattr(self, "globalNumberOfVertices"):
            self.globalNumberOfVertices = self.numberOfVertices
        if not hasattr(self, "globalNumberOfFaces"):
            self.globalNumberOfFaces = self.numberOfFaces
        if not hasattr(self, "globalNumberOfCells"):
            self.globalNumberOfCells = self.numberOfCells

        self.faceCellIDs = self._calcFaceCellIDs() 

        self._setTopology()
        self._setGeometry(scaleLength = 1.)

    @getsetDeprecated
    def _getFaceVertexIDs(self):
        return self.faceVertexIDs

    def _getCellFaceIDs(self):
        return self._cellFaceIDs

    def _setCellFaceIDs(self, newVal):
        """We must notify the helper classes."""
        if hasattr(self, "_topology"): 
            self._topology.cellFaceIDs = newVal
        if hasattr(self, "_geometry"):
            self._geometry.cellFaceIDs = newVal

        self._cellFaceIDs = newVal

    """TODO; This is all to enable `_connectFaces` to work properly. Yet another
    reason why `_connectFaces` should be moved lower down the mesh hierarchy,
    since only periodic grids use them."""
    cellFaceIDs = property(_getCellFaceIDs, _setCellFaceIDs)


    """Topology methods"""
    def _setTopology(self):
        self._topology = _MeshTopology(self.cellFaceIDs, 
                                       self.faceCellIDs, 
                                       self.numberOfCells,
                                       self._maxFacesPerCell,
                                       self) # `self` only for int/ext face calc

    def _setGeometry(self, scaleLength = 1.):
        self._geometry = _MeshGeometry(self,
                                       self.dim,
                                       self.faceVertexIDs,
                                       self.vertexCoords,
                                       self.faceCellIDs,
                                       self.cellFaceIDs,
                                       self.numberOfCells,
                                       self._maxFacesPerCell,
                                       self._cellToFaceOrientations,
                                       scaleLength)
                                      
    @getsetDeprecated
    def setScale(self, scaleLength = 1.):
        return self._setScale(scaleLength)

    def _setScale(self, scaleLength = 1.):
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

    interiorFaces          = property(lambda s: s._topology.interiorFaces)

    def _setExteriorFaces(self, newExtFaces):
        self._topology.exteriorFaces = newExtFaces

    exteriorFaces          = property(lambda s: s._topology.exteriorFaces,
                                      _setExteriorFaces)
    _interiorCellIDs        = property(lambda s: s._topology.interiorCellIDs)
    _exteriorCellIDs        = property(lambda s: s._topology.exteriorCellIDs)
    _cellToFaceOrientations = property(lambda s: s._topology.cellToFaceOrientations)
    _adjacentCellIDs        = property(lambda s: s._topology.adjacentCellIDs)
    _cellToCellIDs          = property(lambda s: s._topology.cellToCellIDs)
    _cellToCellIDsFilled    = property(lambda s: s._topology.cellToCellIDsFilled)

    """geometry properties"""
    _faceAreas                = property(lambda s: s._geometry.scaledFaceAreas)
    faceCenters               = property(lambda s: s._geometry.faceCenters)

    def _setFaceToCellDistances(self, v):
        self._geometry.faceToCellDistances = v

    _faceToCellDistances = property(lambda s: s._geometry.faceToCellDistances,
                                   _setFaceToCellDistances)

    def _setCellDistances(self, v):
        self._geometry.cellDistances = v

    _cellDistances = property(lambda s: s._geometry.scaledCellDistances,
                             _setCellDistances)

    def _setFaceNormals(self, v):
        self._geometry.faceNormals = v

    _faceNormals = property(lambda s: s._geometry.faceNormals,
                           _setFaceNormals)

    cellToFaceDistanceVectors  = property(lambda s: s._geometry.cellToFaceDistanceVectors)
    cellDistanceVectors        = property(lambda s: s._geometry.cellDistanceVectors)
    _orientedFaceNormals       = property(lambda s: s._geometry.orientedFaceNormals)
    cellVolumes                = property(lambda s: s._geometry.scaledCellVolumes)

    @property
    def cellCenters(self):
        from fipy.variables.cellVariable import CellVariable
        return CellVariable(mesh=self, value=self._geometry.scaledCellCenters,
                            rank=1)

    _faceCellToCellNormals    = property(lambda s: s._geometry.faceCellToCellNormals)
    _faceTangents1            = property(lambda s: s._geometry.faceTangents1)
    _faceTangents2            = property(lambda s: s._geometry.faceTangents2)
    _cellToCellDistances      = property(lambda s: s._geometry.scaledCellToCellDistances)
    _cellAreas                 = property(lambda s: s._geometry.cellAreas)
    _cellNormals               = property(lambda s: s._geometry.cellNormals)

    """scaled geometery properties
    
    These should not exist."""
    scale                     = property(lambda s: s._geometry.scale,
                                         _setScale)
    scaledFaceAreas           = property(lambda s: s._geometry.scaledFaceAreas)
    scaledCellVolumes         = property(lambda s: s._geometry.scaledCellVolumes)
    _scaledCellCenters         = property(lambda s: s._geometry.scaledCellCenters)
    scaledFaceToCellDistances = property(lambda s: \
                                         s._geometry.scaledFaceToCellDistances)
    scaledCellDistances       = property(lambda s: \
                                         s._geometry.scaledCellDistances)
    scaledCellToCellDistances = property(lambda s: \
                                         s._geometry.scaledCellToCellDistances)
    _areaProjections          = property(lambda s: \
                                         s._geometry.areaProjections)
    _orientedAreaProjections  = property(lambda s: \
                                         s._geometry.orientedAreaProjections)
    _faceToCellDistanceRatio  = property(lambda s: \
                                         s._geometry.faceToCellDistanceRatio)
    _faceAspectRatios         = property(lambda s: \
                                         s._geometry.faceAspectRatios)  
        
    def __add__(self, other):
        """
        Either translate a `Mesh` or concatenate two `Mesh` objects.
        
            >>> from fipy.meshes import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5]
             [ 0.5  0.5  1.5  1.5]]
             
        If a vector is added to a `Mesh`, a translated `Mesh` is returned
        
            >>> translatedMesh = baseMesh + ((5,), (10,))
            >>> print translatedMesh.cellCenters
            [[  5.5   6.5   5.5   6.5]
             [ 10.5  10.5  11.5  11.5]]

             
        If a `Mesh` is added to a `Mesh`, a concatenation of the two 
        `Mesh` objects is returned
        
            >>> addedMesh = baseMesh + (baseMesh + ((2,), (0,)))
            >>> print addedMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5  2.5  3.5  2.5  3.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]
        
        The two `Mesh` objects need not be properly aligned in order to concatenate them
        but the resulting mesh may not have the intended connectivity
        
            >>> from fipy.meshes.mesh import MeshAdditionError
            >>> addedMesh = baseMesh + (baseMesh + ((3,), (0,))) 
            >>> print addedMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5  3.5  4.5  3.5  4.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]

            >>> addedMesh = baseMesh + (baseMesh + ((2,), (2,)))
            >>> print addedMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5  2.5  3.5  2.5  3.5]
             [ 0.5  0.5  1.5  1.5  2.5  2.5  3.5  3.5]]

        No provision is made to avoid or consolidate overlapping `Mesh` objects
        
            >>> addedMesh = baseMesh + (baseMesh + ((1,), (0,)))
            >>> print addedMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5  1.5  2.5  1.5  2.5]
             [ 0.5  0.5  1.5  1.5  0.5  0.5  1.5  1.5]]
            
        Different `Mesh` classes can be concatenated
         
            >>> from fipy.meshes import Tri2D
            >>> triMesh = Tri2D(dx = 1.0, dy = 1.0, nx = 2, ny = 1)
            >>> triMesh = triMesh + ((2,), (0,))
            >>> triAddedMesh = baseMesh + triMesh
            >>> cellCenters = [[0.5, 1.5, 0.5, 1.5, 2.83333333,  3.83333333,
            ...                 2.5, 3.5, 2.16666667, 3.16666667, 2.5, 3.5],
            ...                [0.5, 0.5, 1.5, 1.5, 0.5, 0.5, 0.83333333, 0.83333333, 
            ...                 0.5, 0.5, 0.16666667, 0.16666667]]
            >>> print numerix.allclose(triAddedMesh.cellCenters,
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
            >>> print numerix.allclose(triAddedMesh.cellCenters,
            ...                        cellCenters)
            True

        `Mesh` concatenation is not limited to 2D meshes
        
            >>> from fipy.meshes import Grid3D
            >>> threeDBaseMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                         nx = 2, ny = 2, nz = 2)
            >>> threeDSecondMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                           nx = 1, ny = 1, nz = 1)
            >>> threeDAddedMesh = threeDBaseMesh + (threeDSecondMesh + ((2,), (0,), (0,)))
            >>> print threeDAddedMesh.cellCenters
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
        
            >>> from fipy.meshes import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.cellCenters
            [[ 0.5  1.5  0.5  1.5]
             [ 0.5  0.5  1.5  1.5]]

        The `factor` can be a scalar
        
            >>> dilatedMesh = baseMesh * 3
            >>> print dilatedMesh.cellCenters
            [[ 1.5  4.5  1.5  4.5]
             [ 1.5  1.5  4.5  4.5]]

        or a vector
        
            >>> dilatedMesh = baseMesh * ((3,), (2,))
            >>> print dilatedMesh.cellCenters
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
     
    """
    TODO: Put this in a PeriodicGrid specific class.
    """
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

           >>> mesh._connectFaces(numerix.nonzero(mesh.facesLeft), numerix.nonzero(mesh.facesRight))

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
        assert numerix.allclose(numerix.take(self._areaProjections, faces0, axis=1),
                               numerix.take(-self._areaProjections, faces1, axis=1))

        ## extract the adjacent cells for both sets of faces
        self.faceCellIDs = self.faceCellIDs.copy()
        ## set the new adjacent cells for `faces0`
        newFaces0 = self.faceCellIDs[0].take(faces0)
        newFaces1 = self.faceCellIDs[0].take(faces1)
        
        self.faceCellIDs[1].put(faces0, newFaces0)
        self.faceCellIDs[0].put(faces0, newFaces1)
        
        ## extract the face to cell distances for both sets of faces
        self.faceToCellDistances = self.faceToCellDistances.copy()
        ## set the new faceToCellDistances for `faces0`
        newDistances0 = self.faceToCellDistances[0].take(faces0)
        newDistances1 = self.faceToCellDistances[0].take(faces1)
        
        self.faceToCellDistances[1].put(faces0, newDistances0)
        self.faceToCellDistances[0].put(faces0, newDistances1)

        ## change the direction of the face normals for faces0
        self.faceNormals = self.faceNormals.copy()
        for dim in range(self.getDim()):
            faceNormals = self.faceNormals[dim].copy()
            numerix.put(faceNormals, faces0, faceNormals.take(faces1))
            self.faceNormals[dim] = faceNormals

        ## Cells that are adjacent to faces1 are changed to point at faces0
        ## get the cells adjacent to faces1
        faceCellIDs = self.faceCellIDs[0].take(faces1)
        ## get all the adjacent faces for those particular cells
        self.cellFaceIDs = self.cellFaceIDs.copy()
        cellFaceIDs = self.cellFaceIDs.take(faceCellIDs, axis=1).copy()
        
        for i in range(cellFaceIDs.shape[0]):
            ## if the faces is a member of faces1 then change the face to point at
            ## faces0
            facesInFaces1 = (cellFaceIDs[i] == faces1)
            cellFaceIDs[i] = (facesInFaces1 * faces0
                              + ~facesInFaces1 * cellFaceIDs[i])
            ## add those faces back to the main self.cellFaceIDs
            self.cellFaceIDs[i].put(faceCellIDs, cellFaceIDs[i])

        ## calculate new topology
        self._setTopology()

        ## calculate new geometry
        self._geometry.handleFaceConnection()
        
        self.scale = self.scale['length']
        
    def _getConcatenableMesh(self):
        return self
        
    def _getAddedMeshValues(self, other, resolution=1e-2):
        """Calculate the parameters to define a concatenation of `other` with `self`
        
        :Parameters:
          - `other`: The :class:`~fipy.meshes.Mesh` to concatenate with `self`
          - `resolution`: How close vertices have to be (relative to the smallest 
            cell-to-cell distance in either mesh) to be considered the same

        :Returns:
          A `dict` with 3 elements: the new mesh `vertexCoords`, 
          `faceVertexIDs`, and `cellFaceIDs`.
        """
        
        selfc = self._getConcatenableMesh()
        other = other._getConcatenableMesh()

        ## check dimensions
        if self.dim != other.dim:
            raise MeshAdditionError, "Dimensions do not match"
            
        ## compute vertex correlates

        ## only try to match exterior (X) vertices
        self_Xvertices = numerix.unique(selfc._getFaceVertexIDs()
                                        .filled()[..., selfc.getExteriorFaces()]
                                        .flatten().value)
        other_Xvertices = numerix.unique(other._getFaceVertexIDs()
                                         .filled()[..., other.getExteriorFaces()]
                                         .flatten().value)

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
        close = distance < resolution * min(selfc._cellToCellDistances.min(), 
                                            other._cellToCellDistances.min()).value
        vertexCorrelates = numerix.array((self_Xvertices[closest[close]],
                                          other_Xvertices[close]))
        
        # warn if meshes don't touch, but allow it
        if (selfc.numberOfVertices > 0 
            and other.numberOfVertices > 0 
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
        self_matchingFaces = numerix.in1d(self_faceVertexIDs.value, 
                                          vertexCorrelates[0]).reshape(self_faceVertexIDs.shape).all(axis=0).nonzero()[0]

        # want other's Faces for which all faceVertexIDs are in vertexCorrelates
        other_matchingFaces = numerix.in1d(other_faceVertexIDs.value, 
                                           vertexCorrelates[1]).reshape(other_faceVertexIDs.shape).all(axis=0).nonzero()[0]
                                           
        # map other's Vertex IDs to new Vertex IDs, 
        # accounting for overlaps with self's Vertex IDs
        vertex_map = numerix.empty(other.numberOfVertices, dtype=int)
        verticesToAdd = numerix.delete(numerix.arange(other.numberOfVertices), vertexCorrelates[1])
        vertex_map[verticesToAdd] = numerix.arange(other.numberOfVertices - len(vertexCorrelates[1])) + selfc.numberOfVertices
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
            other_faceHash = vertex_map[other_faceVertexIDs[..., other_matchingFaces].value]
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
        if (selfc.numberOfFaces > 0 
            and other.numberOfFaces > 0 
            and faceCorrelates.shape[-1] == 0):
            import warnings
            warnings.warn("Faces are not aligned", UserWarning, stacklevel=4)

        # map other's Face IDs to new Face IDs, 
        # accounting for overlaps with self's Face IDs
        face_map = numerix.empty(other._numberOfFaces, dtype=int)
        facesToAdd = numerix.delete(numerix.arange(other._numberOfFaces), faceCorrelates[1])
        face_map[facesToAdd] = numerix.arange(other._numberOfFaces - len(faceCorrelates[1])) + selfc._numberOfFaces
        face_map[faceCorrelates[1]] = faceCorrelates[0]
        
        other_faceVertexIDs = vertex_map[other.faceVertexIDs[..., facesToAdd].value]
        
        # ensure that both sets of cellFaceIDs have the same maximum number of (masked) elements
        self_cellFaceIDs = selfc.cellFaceIDs
        other_cellFaceIDs = face_map[other.cellFaceIDs.value]
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
        newmesh = Mesh(vertexCoords=newCoords,
                       faceVertexIDs=self.faceVertexIDs,
                       cellFaceIDs=self.cellFaceIDs)
        return newmesh

    """calc Topology methods"""

    @getsetDeprecated
    def _getNumberOfFacesPerCell(self):
        return self._numberOfFacesPerCell

    @property
    def _numberOfFacesPerCell(self):
        cellFaceIDs = self.cellFaceIDs
        if type(cellFaceIDs) is type(MA.array(0)):
            ## bug in count returns float values when there is no mask
            return numerix.array(cellFaceIDs.count(axis=0), 'l')
        else:
            return self._maxFacesPerCell * numerix.ones(cellFaceIDs.shape[-1], 'l')
      
    """
    TODO: Does this really belong in mesh? I don't think so.
    """
    def _calcFaceCellIDs(self):
        from fipy.variables.faceVariable import FaceVariable
        class _FaceCellIDsVariable(FaceVariable):
            def __init__(self, mesh):
                FaceVariable.__init__(self, mesh=mesh,
                                      elementshape=(2,))
                self._requires(mesh.cellFaceIDs)
                self._requires(mesh.numberOfFaces)
                
            def _calcValue(self):
                array = MA.array(MA.indices(self.mesh.cellFaceIDs.shape, 'l')[1], 
                                 mask=self.mesh.cellFaceIDs.getMask())
                faceCellIDs = MA.zeros((2, self.mesh.numberOfFaces), 'l')

                numerix.put(faceCellIDs[0], self.mesh.cellFaceIDs[::-1,::-1], array[::-1,::-1])
                numerix.put(faceCellIDs[1], self.mesh.cellFaceIDs, array)
                return MA.sort(MA.array(faceCellIDs,
                                        mask = ((False,) * self.mesh.numberOfFaces, 
                                                (faceCellIDs[0] == faceCellIDs[1]))),
                               axis=0)
                               
        return _FaceCellIDsVariable(mesh=self)

    def _calcVertexFaceIDs(self):
        from fipy.variables.vertexVariable import _VertexVariable
        class _VertexFaceIDsVariable(_VertexVariable):
            def __init__(self, mesh):
                _VertexVariable.__init__(self, mesh=mesh,
                                         elementshape=(2,))
                self._requires(mesh.vertexFaceIDs)
##                 self._requires(mesh.numberOfFaces)
                
            def _calcValue(self):
                array = MA.array(MA.indices(self.mesh.faceVertexIDs.shape, 'l')[1], 
                                 mask=self.mesh.faceVertexIDs.getMask())
                vertexFaceIDs = MA.zeros((2, self.mesh.numberOfVertices), 'l')
                firstRow = self.faceCellIDs[0]
                secondRow = self.faceCellIDs[1]
                numerix.put(firstRow, self.cellFaceIDs[::-1,::-1], array[::-1,::-1])
                numerix.put(secondRow, self.cellFaceIDs, array)
                
                mask = ((False,) * self.numberOfFaces, (firstRow == secondRow))
                return MA.sort(MA.array(self.faceCellIDs, mask = mask),
                               axis=0)
                               
        return _VertexFaceIDsVariable(mesh=self)

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
    
    @getsetDeprecated
    def getVertexCoords(self):
        """TODO: replace this with a warning."""
        if hasattr(self, 'vertexCoords'):
            return self.vertexCoords
        else:
            return self._createVertices()

    @getsetDeprecated
    def getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        TODO: replace with a warning.
        """
        return self.exteriorFaces
            
    @getsetDeprecated
    def getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        TODO: replace with a warning.
        """
        return self.interiorFaces
    
    @getsetDeprecated
    def getInteriorFaceIDs(self):
        return self.interiorFaceIDs

    @property
    def interiorFaceIDs(self):
        if not hasattr(self, '_interiorFaceIDs'):
            self._interiorFaceIDs = numerix.nonzero(self.interiorFaces)[0]
        return self._interiorFaceIDs

    @getsetDeprecated
    def getInteriorFaceCellIDs(self):
        return self.interiorFaceCellIDs

    @property
    def interiorFaceCellIDs(self):
        if not hasattr(self, '_interiorFaceCellIDs'):
            ## Commented line is better, but doesn't work for zero length arrays
            ##  self.interiorFaceCellIDs = self.getFaceCellIDs()[..., self.getInteriorFaceIDs()]
            self._interiorFaceCellIDs = numerix.take(self.faceCellIDs,
                                                     self.interiorFaceIDs, axis=1)
        return self._interiorFaceCellIDs
     
    @getsetDeprecated
    def getFaceCellIDs(self):
        return self.faceCellIDs

    @getsetDeprecated
    def _getMaxFacesPerCell(self):
        return self._maxFacesPerCell

    @property
    def _maxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    @getsetDeprecated
    def _getExteriorCellIDs(self):
        """ Why do we have this?!? It's only used for testing against itself? """
        return self._exteriorCellIDs

    @getsetDeprecated
    def _getInteriorCellIDs(self):
        """ Why do we have this?!? It's only used for testing against itself? """
        return self._interiorCellIDs

    @getsetDeprecated
    def _getCellFaceOrientations(self):
        return self._cellToFaceOrientations

    @getsetDeprecated
    def getNumberOfCells(self):
        return self.numberOfCells

    def _isOrthogonal(self):
        return False
    
    @getsetDeprecated
    def _getNumberOfVertices(self):
        return self.numberOfVertices

    @getsetDeprecated
    def _getAdjacentCellIDs(self):
        return self._adjacentCellIDs

    @getsetDeprecated
    def getDim(self):
        return self.dim

    @getsetDeprecated
    def _getGlobalNonOverlappingCellIDs(self):
        return self._globalNonOverlappingCellIDs

    @property
    def _globalNonOverlappingCellIDs(self):
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

    @getsetDeprecated
    def _getGlobalOverlappingCellIDs(self):
        return self._globalOverlappingCellIDs

    @property
    def _globalOverlappingCellIDs(self):
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

    @getsetDeprecated
    def _getLocalNonOverlappingCellIDs(self):
        return self._localNonOverlappingCellIDs

    @property
    def _localNonOverlappingCellIDs(self):
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

    @getsetDeprecated
    def _getLocalOverlappingCellIDs(self):
        return self._localOverlappingCellIDs

    @property
    def _localOverlappingCellIDs(self):
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

    @getsetDeprecated
    def _getGlobalNonOverlappingFaceIDs(self):
        return self._globalNonOverlappingFaceIDs

    @property
    def _globalNonOverlappingFaceIDs(self):
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

    @getsetDeprecated
    def _getGlobalOverlappingFaceIDs(self):
        return self._globalOverlappingFaceIDs

    @property
    def _globalOverlappingFaceIDs(self):
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

    @getsetDeprecated
    def _getLocalNonOverlappingFaceIDs(self):
        return self._localNonOverlappingFaceIDs

    @property
    def _localNonOverlappingFaceIDs(self):
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

    @getsetDeprecated
    def _getLocalOverlappingFaceIDs(self):
        return self._localOverlappingFaceIDs

    @property
    def _localOverlappingFaceIDs(self):
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

    @getsetDeprecated
    def getFacesLeft(self):
        return self.facesLeft

    @property
    def facesLeft(self):
        """
        Return face on left boundary of Grid1D as list with the
        x-axis running from left to right.

            >>> from fipy import Grid2D, Grid3D
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((21, 25), 
            ...                              numerix.nonzero(mesh.facesLeft)[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> print parallel.procID > 0 or numerix.allequal((9, 13), 
            ...                              numerix.nonzero(mesh.facesLeft)[0])
            True

        """
        x = self.faceCenters[0]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=x == _madmin(x))

    @getsetDeprecated
    def getFacesRight(self):
        return self.facesRight

    @property
    def facesRight(self):
        """
        Return list of faces on right boundary of Grid3D with the
        x-axis running from left to right. 

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((24, 28), 
            ...                              numerix.nonzero(mesh.facesRight)[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)    
            >>> print parallel.procID > 0 or numerix.allequal((12, 16), 
            ...                                               numerix.nonzero(mesh.facesRight)[0])
            True
            
        """
        x = self.faceCenters[0]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=x == _madmax(x))

    @getsetDeprecated
    def getFacesBottom(self):
        return self.facesBottom

    @property
    def facesBottom(self):
        """
        Return list of faces on bottom boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((12, 13, 14), 
            ...                              numerix.nonzero(mesh.facesBottom)[0])
            1
            >>> x, y, z = mesh.faceCenters
            >>> print parallel.procID > 0 or numerix.allequal((12, 13), 
            ...                              numerix.nonzero(mesh.facesBottom & (x < 1))[0])
            1
            
        """
        y = self.faceCenters[1]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=y == _madmin(y))


    getFacesDown = getFacesBottom
    facesDown = facesBottom

    @getsetDeprecated
    def getFacesTop(self):
        return self.facesTop

    @property
    def facesTop(self):
        """
        Return list of faces on top boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((18, 19, 20), 
            ...                              numerix.nonzero(mesh.facesTop)[0])
            True
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)        
            >>> print parallel.procID > 0 or numerix.allequal((6, 7, 8), 
            ...                              numerix.nonzero(mesh.facesTop)[0])
            True
            
        """
        y = self.faceCenters[1]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=y == _madmax(y))

    getFacesUp = getFacesTop
    facesUp = facesTop

    @getsetDeprecated
    def getFacesBack(self):
        return self.facesBack

    @property
    def facesBack(self):
        """
        Return list of faces on back boundary of Grid3D with the
        z-axis running from front to back. 

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((6, 7, 8, 9, 10, 11), 
            ...                              numerix.nonzero(mesh.facesBack)[0])
            True

        """
        z = self.faceCenters[2] 
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=z == _madmax(z))

    @getsetDeprecated
    def getFacesFront(self):
        return self.facesFront

    @property
    def facesFront(self):
        """
        Return list of faces on front boundary of Grid3D with the
        z-axis running from front to back. 

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> from fipy.tools import parallel
            >>> print parallel.procID > 0 or numerix.allequal((0, 1, 2, 3, 4, 5), 
            ...                              numerix.nonzero(mesh.facesFront)[0])
            True

        """
        z = self.faceCenters[2]
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=z == _madmin(z))
    
    @getsetDeprecated
    def _getNumberOfFaces(self):
        return self.numberOfFaces

    @getsetDeprecated
    def _getCellToCellIDs(self):
        return self._cellToCellIDs

    @getsetDeprecated
    def _getCellToCellIDsFilled(self):
        return self._cellToCellIDsFilled
     
    """get geometry methods"""

    @getsetDeprecated
    def _getFaceAreas(self):
        return self._faceAreas

    @getsetDeprecated
    def _getFaceNormals(self):
        return self._faceNormals

    @getsetDeprecated
    def _getFaceCellToCellNormals(self):
        return self._faceCellToCellNormals
        
    @getsetDeprecated
    def getCellVolumes(self):
        return self.cellVolumes

    """
    This shit has GOT TO GO.

    This is yet another repercussion of UniformGrid inheriting from mesh.
    """
    def _getCellCenters(self):
        return self._scaledCellCenters
        
    @getsetDeprecated
    def getCellCenters(self):
        return self.cellCenters

    @getsetDeprecated
    def _getFaceToCellDistances(self):
        return self._faceToCellDistances

    @getsetDeprecated
    def _getCellDistances(self):
        return self._cellDistances

    @getsetDeprecated
    def _getFaceToCellDistanceRatio(self):
        return self._faceToCellDistanceRatio

    @getsetDeprecated
    def _getOrientedAreaProjections(self):
        return self._orientedAreaProjections

    @getsetDeprecated
    def _getAreaProjections(self):
        return self._areaProjections

    @getsetDeprecated
    def _getOrientedFaceNormals(self):
        return self._orientedFaceNormals

    @getsetDeprecated
    def _getFaceTangents1(self):
        return self._faceTangents1

    @getsetDeprecated
    def _getFaceTangents2(self):
        return self._faceTangents2
        
    @getsetDeprecated
    def _getFaceAspectRatios(self):
        return self._faceAspectRatios
    
    @getsetDeprecated
    def _getCellToCellDistances(self):
        return self._cellToCellDistances

    @getsetDeprecated
    def _getCellNormals(self):
        return self._cellNormals

    @getsetDeprecated
    def _getCellAreas(self):
        return self._cellAreas

    @property
    def _cellAreaProjections(self):
        return self._cellNormals * self._cellAreas

    @getsetDeprecated
    def _getCellAreaProjections(self):
        return self._cellAreaProjections
         
    @getsetDeprecated
    def getFaceCenters(self):
        return self.faceCenters

    @getsetDeprecated
    def _getOrderedCellVertexIDs(self):
        return self._orderedCellVertexIDs

    @property
    def _orderedCellVertexIDs(self):
        return self._cellVertexIDs

    @getsetDeprecated
    def _getCellDistanceNormals(self):
        return self._cellDistanceNormals

    @property
    def _cellDistanceNormals(self):
        return self._cellDistanceNormals/ self._cellDistances
        
    @getsetDeprecated
    def _getCellVertexIDs(self):
        return self._cellVertexIDs

    @property
    def _cellVertexIDs(self):
        ## Get all the vertices from all the faces for each cell
        cellFaceVertices = numerix.take(self.faceVertexIDs, self.cellFaceIDs, axis=1)

        ## get a sorted list of vertices for each cell 
        cellVertexIDs = numerix.reshape(cellFaceVertices, (-1, self.numberOfCells))
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)

        tmp = cellVertexIDs[:-1].masked(cellVertexIDs[:-1] == cellVertexIDs[1:])
        cellVertexIDs = tmp.append(cellVertexIDs[-1, numerix.newaxis], axis=0)
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)
        
        ## resize the array to remove extra masked values
        if cellVertexIDs.shape[-1] == 0:
            length = 0
        else:
            length = cellVertexIDs.getMask().sum(axis=0).min()
        return cellVertexIDs[length:][::-1]

# #     Below is an ordered version of _getCellVertexIDs()
# #     It works for the test case in this file (other than the ordering, obviously)
# #     I've left it in as it may be useful when we need ordered vertices for cells
# #   
# #     def _getOrderedCellVertexIDs(self):
# # 
# #         ## Get all the vertices from all the faces for each cell
# #         from fipy.tools.numerix import take
# #         cellFaceVertices = take(self.faceVertexIDs, self.cellFaceIDs)
# # 
# #         ## get a sorted list of vertices for each cell
# #         NCells = self.numberOfCells
# #         cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
# #         newmask = MA.getmaskarray(cellVertexIDs).copy()
# # 
# #         for i in range(len(cellVertexIDs[0]) - 1):
# #             for j in range(len(cellVertexIDs[0]))[i + 1:]:
# # 
# #                 newmask[:,j] = MA.where(newmask[:,j],
# #                                         newmask[:,j],
# #                                         MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,j],
# #                                                  1,
# #                                                  newmask[:,j]))
# # 
# # 
# #         cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)
# # 
# #         for i in range(len(cellVertexIDs[0]) - 1):
# #             j = i + 1
# #             while j < len(cellVertexIDs[0]):
# #                 tmp = cellVertexIDs[:]
# #                 tmp[:, i] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
# #                                      MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
# #                                               cellVertexIDs[:, i],
# #                                               cellVertexIDs[:, j]),
# #                                      cellVertexIDs[:, i])
# #                                                       
# #                 tmp[:, j] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
# #                                      MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
# #                                               cellVertexIDs[:,j],
# #                                               cellVertexIDs[:,i]),
# #                                      cellVertexIDs[:, j])
# # 
# #                 cellVertexIDs = tmp[:]
# # 
# #                 j += 1
# # 
# # 
# #         length = len(cellVertexIDs[0]) - min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
# #         return cellVertexIDs[:, :length]


    """scaling"""

## ##     def setScale(self, value = 1.):
## ##         self.scale = 1.
    
##    def setScale(self):
##	self.scale = PhysicalField(value = 1.)
##        self.faceAreas = self.faceAreas * self.scale**2
##        self.faceCenters = self.faceCenters * self.scale
##        self.cellVolumes = self.cellVolumes * self.scale**3
##        self.cellCenters = self.cellCenters * self.scale
##        self.faceToCellDistances = self.faceToCellDistances * self.scale
##        self.cellDistances = self.cellDistances * self.scale
##        self.areaProjections = self.areaProjections * self.scale**2

    def getNearestCell(self, point):
        return self._getCellsByID([self._getNearestCellID(point)])[0]

    def _getNearestCellID(self, points):
        """
        Test cases

           >>> from fipy import *
           >>> m0 = Grid2D(dx=(.1, 1., 10.), dy=(.1, 1., 10.))
           >>> m1 = Grid2D(nx=2, ny=2, dx=5., dy=5.)
           >>> print m0._getNearestCellID(m1.cellCenters.globalValue)
           [4 5 7 8]
           
        """
        return numerix.nearest(data=self.cellCenters.globalValue, points=points)
        
    def _subscribe(self, var):
        if not hasattr(self, 'subscribedVariables'):
            self.subscribedVariables = []

        # we retain a weak reference to avoid a memory leak
        # due to circular references between the subscriber
        # and the subscribee
        import weakref
        self.subscribedVariables.append(weakref.ref(var))

    def getSubscribedVariables(self):
        if not hasattr(self, 'subscribedVariables'):
            self.subscribedVariables = []
           
        self.subscribedVariables = [sub for sub in self.subscribedVariables if sub() is not None]
       
        return self.subscribedVariables


    """pickling"""

    def __getstate__(self):
        dict = {
            'vertexCoords' : self.vertexCoords.value *  self.scale['length'],            
            'faceVertexIDs' : numerix.ma.filled(self.faceVertexIDs.value),
            'cellFaceIDs' : numerix.ma.filled(self.cellFaceIDs.value) }
        return dict

    def __setstate__(self, dict):
        self._concatenatedClass.__init__(self, **dict)
     
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
            >>> print numerix.allequal(faceCellIds, mesh.faceCellIDs)
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> print numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> print numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> print numerix.allclose(faceNormals, mesh._faceNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> print numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> print numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellCenters = numerix.array(((dx/2., dx+dx/3.),
            ...                              (dy/2.,    dy/3.),
            ...                              (dz/2.,    dz/2.)))
            >>> print numerix.allclose(cellCenters, mesh.cellCenters, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> d1 = numerix.sqrt((dx / 3.)**2 + (dy / 6.)**2)
            >>> d2 = numerix.sqrt((dx / 6.)**2 + (dy / 3.)**2)
            >>> d3 = numerix.sqrt((dx / 6.)**2 + (dy / 6.)**2)
            >>> d4 = numerix.sqrt((5 * dx / 6.)**2 + (dy / 6.)**2)
            >>> faceToCellDistances = MA.masked_values(((dz / 2., dz / 2., dx / 2., dx / 2., dy / 2., dy / 2., dz / 2., dz / 2., d2, d3),
            ...                                         (     -1,      -1,      -1,      d1,      -1,      -1,      -1,      -1, -1, -1)), -1)
            >>> print numerix.allclose(faceToCellDistances, mesh._faceToCellDistances, atol = 1e-10, rtol = 1e-10)
            True
                                              
            >>> cellDistances = numerix.array((dz / 2., dz / 2., dx / 2.,
            ...                                d4,
            ...                                dy / 2., dy / 2., dz / 2., dz / 2.,
            ...                                d2,
            ...                                d3))
            >>> print numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> print numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> print numerix.allequal(cellToCellIDs, mesh._cellToCellIDs)
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> print numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> print numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> print numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)
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
            >>> print numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)
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
            >>> print numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)
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
            >>> print numerix.allclose(cellVertexIDs, mesh._cellVertexIDs)
            1


            >>> from fipy.tools import dump            
            >>> (f, filename) = dump.write(mesh, extension = '.gz')
            >>> unpickledMesh = dump.read(filename, f)

            >>> print numerix.allequal(mesh.cellCenters, unpickledMesh.cellCenters)
            True

            >>> dx = 1.
            >>> dy = 1.
            >>> nx = 10
            >>> ny = 2
            >>> from fipy.meshes import Grid2D
            >>> gridMesh = Grid2D(dx, dy, nx, ny)
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx, dy, nx, 1) + [[dx*nx], [0]]
            >>> bigMesh = gridMesh + triMesh
            >>> x, y = bigMesh.cellCenters
            >>> from fipy.variables.cellVariable import CellVariable
            >>> volumes = CellVariable(mesh=bigMesh, value=1.)
            >>> volumes[x > dx * nx] = 0.25
            >>> print numerix.allclose(bigMesh.cellVolumes, volumes)
            True
            
            Following test was added due to a bug in adding UniformGrids.

            >>> from fipy.meshes.uniformGrid1D import UniformGrid1D
            >>> a = UniformGrid1D(nx=10) + (10,)
            >>> print a.cellCenters
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            >>> b = 10 + UniformGrid1D(nx=10)
            >>> print b.cellCenters
            [[ 10.5  11.5  12.5  13.5  14.5  15.5  16.5  17.5  18.5  19.5]]
            
            >>> from fipy.tools import parallel
            >>> if parallel.Nproc == 1:
            ...     c =  UniformGrid1D(nx=10) + (UniformGrid1D(nx=10) + 10)
            >>> print (parallel.Nproc > 1 
            ...        or numerix.allclose(c.cellCenters[0],
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
        cvi = self._orderedCellVertexIDs.swapaxes(0,1)
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
        
        points = self.faceCenters
        points = self._toVTK3D(points)
        ug = tvtk.UnstructuredGrid(points=points)
        
        num = len(points)
        counts = numerix.array([1] * num)[..., numerix.newaxis]
        cells = numerix.arange(self.numberOfFaces)[..., numerix.newaxis]
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
        from fipy.variables import Variable
        if isinstance(arr, Variable):
            arr = arr.getValue()
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
