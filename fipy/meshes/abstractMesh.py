#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "abstractMesh.py"
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

__docformat__ = 'restructuredtext'

__all__ = ["AbstractMesh"]

from fipy.tools import serialComm
from fipy.tools import numerix
from fipy.tools.decorators import deprecate
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.meshes.representations.abstractRepresentation import _AbstractRepresentation
from fipy.meshes.topologies.abstractTopology import _AbstractTopology

class MeshAdditionError(Exception):
    pass

class AbstractMesh(object):
    """
    A class encapsulating all commonalities among meshes in FiPy.
    """

    def __init__(self, communicator, _RepresentationClass=_AbstractRepresentation, _TopologyClass=_AbstractTopology):
        self.communicator = communicator
        self.representation = _RepresentationClass(mesh=self)
        self.topology = _TopologyClass(mesh=self)

    def _setTopology(self):
        raise NotImplementedError

    def _setGeometry(self):
        raise NotImplementedError

    """
    Scale business
    """
    def _setScale(self, scaleLength = 1.):
        """
        Sets scale of geometry.

        :Parameters:
          - `scaleLength`: The desired scale length.
        """
        self._scale['length'] = scaleLength

    scale = property(lambda s: s._scale, _setScale)

    def _calcScaleArea(self):
        raise NotImplementedError

    def _calcScaleVolume(self):
        raise NotImplementedError

    def _getPointToCellDistances(self, point):
        tmp = self.cellCenters - PhysicalField(point)
        return numerix.sqrtDot(tmp, tmp)

    def getNearestCell(self, point):
        return self._getCellsByID([self._getNearestCellID(point)])[0]

    def _getCellFaceIDsInternal(self):
        return self._cellFaceIDs

    def _setCellFaceIDsInternal(self, newVal):
        self._cellFaceIDs = newVal

    """This is to enable `_connectFaces` to work properly."""
    cellFaceIDs = property(_getCellFaceIDsInternal, _setCellFaceIDsInternal)

    """Topology properties"""

    interiorFaces = property(lambda s: s._interiorFaces)

    def _setExteriorFaces(self, newExtFaces):
        self._exteriorFaces = newExtFaces

    exteriorFaces           = property(lambda s: s._exteriorFaces,
                                      _setExteriorFaces)

    @property
    def _isOrthogonal(self):
        return self.topology._isOrthogonal

    """Geometry properties"""

    @property
    def faceCenters(self):
        from fipy.variables.faceVariable import FaceVariable
        return FaceVariable(mesh=self, value=self._faceCenters,
                            rank=1)

    cellToFaceDistanceVectors  = property(lambda s: s._cellToFaceDistanceVectors)
    cellDistanceVectors        = property(lambda s: s._cellDistanceVectors)
    cellVolumes                = property(lambda s: s._scaledCellVolumes)

    @property
    @deprecate(new_name="faceNormals", version=3.1)
    def _faceNormals(self):
        return self.faceNormals

    @property
    def cellCenters(self):
        from fipy.variables.cellVariable import CellVariable
        return CellVariable(mesh=self, value=self._scaledCellCenters,
                            rank=1)

    @property
    def x(self):
        """
        Equivalent to using :attr:`cellCenters`\ ``[0]``.

        >>> from fipy import *
        >>> print Grid1D(nx=2).x
        [ 0.5  1.5]
        """
        return self.cellCenters[0]

    @property
    def y(self):
        """
        Equivalent to using :attr:`cellCenters`\ ``[1]``.

        >>> from fipy import *
        >>> print Grid2D(nx=2, ny=2).y
        [ 0.5  0.5  1.5  1.5]
        >>> print Grid1D(nx=2).y
        Traceback (most recent call last):
          ...
        AttributeError: 1D meshes do not have a "y" attribute.
        """
        if self.dim > 1:
            return self.cellCenters[1]
        else:
            raise AttributeError, '1D meshes do not have a "y" attribute.'

    @property
    def z(self):
        """
        Equivalent to using :attr:`cellCenters`\ ``[2]``.

        >>> from fipy import *
        >>> print Grid3D(nx=2, ny=2, nz=2).z
        [ 0.5  0.5  0.5  0.5  1.5  1.5  1.5  1.5]
        >>> print Grid2D(nx=2, ny=2).z
        Traceback (most recent call last):
          ...
        AttributeError: 1D and 2D meshes do not have a "z" attribute.
        """
        if self.dim > 2:
            return self.cellCenters[2]
        else:
            raise AttributeError, '1D and 2D meshes do not have a "z" attribute.'

    @property
    def extents(self):
        ext = dict(min=[], max=[])

        for d in range(self.dim):
            X = numerix.take(self.vertexCoords[d], self._orderedCellVertexIDs)
            ext['min'].append(X.min())
            ext['max'].append(X.max())

        return ext

    """scaled geometery properties

    These should not exist."""
    scaledFaceAreas           = property(lambda s: s._scaledFaceAreas)
    scaledCellVolumes         = property(lambda s: s._scaledCellVolumes)
    scaledFaceToCellDistances = property(lambda s: s._scaledFaceToCellDistances)
    scaledCellDistances       = property(lambda s: s._scaledCellDistances)
    scaledCellToCellDistances = property(lambda s: s._scaledCellToCellDistances)

    def _connectFaces(self, faces0, faces1):
        """

        Merge faces on the same mesh. This is used to create periodic
        meshes. The first list of faces, `faces1`, will be the faces
        that are used to add to the matrix diagonals. The faces in
        `faces2` will not be used. They aren't deleted but their
        adjacent cells are made to point at `faces1`. The list
        `faces2` are not altered, they still remain as members of
        exterior faces.

           >>> from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D
           >>> mesh = NonUniformGrid2D(nx = 2, ny = 2, dx = 1., dy = 1.)

           >>> print (mesh.cellFaceIDs == [[0, 1, 2, 3],
           ...                             [7, 8, 10, 11],
           ...                             [2, 3, 4, 5],
           ...                             [6, 7, 9, 10]]).flatten().all() # doctest: +PROCESSOR_0
           True

           >>> mesh._connectFaces(numerix.nonzero(mesh.facesLeft), numerix.nonzero(mesh.facesRight))

           >>> print (mesh.cellFaceIDs == [[0, 1, 2, 3],
           ...                             [7, 6, 10, 9],
           ...                             [2, 3, 4, 5],
           ...                             [6, 7, 9, 10]]).flatten().all() # doctest: +PROCESSOR_0
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
        faceCellIDs0 = self.faceCellIDs[0]
        faceCellIDs1 = self.faceCellIDs[1]
        ## set the new adjacent cells for `faces0`
        MA.put(faceCellIDs1, faces0, MA.take(faceCellIDs0, faces0))
        MA.put(faceCellIDs0, faces0, MA.take(faceCellIDs0, faces1))
        self.faceCellIDs[0] = faceCellIDs0
        self.faceCellIDs[1] = faceCellIDs1

        ## extract the face to cell distances for both sets of faces
        faceToCellDistances0 = self._faceToCellDistances[0]
        faceToCellDistances1 = self._faceToCellDistances[1]
        ## set the new faceToCellDistances for `faces0`
        MA.put(faceToCellDistances1, faces0, MA.take(faceToCellDistances0, faces0))
        MA.put(faceToCellDistances0, faces0, MA.take(faceToCellDistances0, faces1))

        self._faceToCellDistances[0] = faceToCellDistances0
        self._faceToCellDistances[1] = faceToCellDistances1

        ## calculate new cell distances and add them to faces0
        numerix.put(self._cellDistances, faces0, MA.take(faceToCellDistances0 + faceToCellDistances1, faces0))

        ## change the direction of the face normals for faces0
        for dim in range(self.dim):
            faceNormals = self.faceNormals[dim].copy()
            numerix.put(faceNormals, faces0, MA.take(faceNormals, faces1))
            self.faceNormals[dim] = faceNormals

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
            numerix.put(self.cellFaceIDs[i], faceCellIDs, cellFaceIDs[i])

        ## calculate new topology
        self._setTopology()

        ## calculate new geometry
        self._handleFaceConnection()

        self.scale = self.scale['length']

    @property
    def _concatenableMesh(self):
        raise NotImplementedError

    def _translate(self, vector):
        raise NotImplementedError

    def _getAddedMeshValues(self, other, resolution=1e-2):
        """Calculate the parameters to define a concatenation of `other` with `self`

        :Parameters:
          - `other`: The :class:`~fipy.meshes.Mesh` to concatenate with `self`
          - `resolution`: How close vertices have to be (relative to the smallest
            cell-to-cell distance in either mesh) to be considered the same

        :Returns:
          A `dict` with 3 elements: the new mesh vertexCoords, faceVertexIDs, and cellFaceIDs.
        """

        selfc = self._concatenableMesh
        otherc = other._concatenableMesh

        selfNumFaces = selfc.faceVertexIDs.shape[-1]
        selfNumVertices = selfc.vertexCoords.shape[-1]
        otherNumFaces = otherc.faceVertexIDs.shape[-1]
        otherNumVertices = otherc.vertexCoords.shape[-1]
        ## check dimensions
        if(selfc.vertexCoords.shape[0] != otherc.vertexCoords.shape[0]):
            raise MeshAdditionError, "Dimensions do not match"

        ## compute vertex correlates

        """
        from fipy.tools.debug import PRINT
        PRINT("selfNumFaces", selfNumFaces)
        PRINT("otherNumFaces", otherNumVertices)
        PRINT("selfNumVertices", selfNumVertices)
        PRINT("otherNumVertices", otherNumVertices)

        from fipy.tools.debug import PRINT
        from fipy.tools.debug import PRINT
        PRINT("otherExt", otherc.exteriorFaces.value)
        raw_input()
        PRINT("selfExt", selfc.exteriorFaces.value)

        PRINT("self filled", selfc.faceVertexIDs.filled())
        PRINT("othe filled", otherc.faceVertexIDs.filled())
        raw_input()

        PRINT("selfc.faceVertexIDs.filled()\n",selfc.faceVertexIDs.filled())
        PRINT("flat\n",selfc.faceVertexIDs.filled()[...,
            selfc.exteriorFaces.value].flatten())
        PRINT("selfc.exteriorFaces.value\n",selfc.exteriorFaces.value)
        PRINT("extfaces type", type(selfc.exteriorFaces))
        PRINT("extfaces mesh", selfc.exteriorFaces.mesh)
        """

        ## only try to match along the operation manifold
        if hasattr(self, "opManifold"):
            self_faces = self.opManifold(selfc)
        else:
            self_faces = selfc.exteriorFaces.value
        if hasattr(other, "opManifold"):
            other_faces = other.opManifold(otherc)
        else:
            other_faces = otherc.exteriorFaces.value

        ## only try to match exterior (X) vertices
        self_Xvertices = numerix.unique(selfc.faceVertexIDs.filled()[...,
            self_faces].flatten())
        other_Xvertices = numerix.unique(otherc.faceVertexIDs.filled()[...,
            other_faces].flatten())

        self_XvertexCoords = selfc.vertexCoords[..., self_Xvertices]
        other_XvertexCoords = otherc.vertexCoords[..., other_Xvertices]

        closest = numerix.nearest(self_XvertexCoords, other_XvertexCoords)

        # just because they're closest, doesn't mean they're close
        tmp = self_XvertexCoords[..., closest] - other_XvertexCoords
        distance = numerix.sqrtDot(tmp, tmp)
        # only want vertex pairs that are 100x closer than the smallest
        # cell-to-cell distance
        close = distance < resolution * min(selfc._cellToCellDistances.min(),
                                            otherc._cellToCellDistances.min())
        vertexCorrelates = numerix.array((self_Xvertices[closest[close]],
                                          other_Xvertices[close]))

        # warn if meshes don't touch, but allow it
        if (selfc._numberOfVertices > 0
            and otherc._numberOfVertices > 0
            and vertexCorrelates.shape[-1] == 0):
            import warnings
            warnings.warn("Vertices are not aligned", UserWarning, stacklevel=4)

        ## compute face correlates

        # ensure that both sets of faceVertexIDs have the same maximum number of (masked) elements
        self_faceVertexIDs = selfc.faceVertexIDs
        other_faceVertexIDs = otherc.faceVertexIDs

        diff = self_faceVertexIDs.shape[0] - other_faceVertexIDs.shape[0]
        if diff > 0:
            other_faceVertexIDs = numerix.append(other_faceVertexIDs,
                                                 -1 * numerix.ones((diff,)
                                                                   + other_faceVertexIDs.shape[1:], 'l'),
                                                 axis=0)
            other_faceVertexIDs = MA.masked_values(other_faceVertexIDs, -1)
        elif diff < 0:
            self_faceVertexIDs = numerix.append(self_faceVertexIDs,
                                                -1 * numerix.ones((-diff,)
                                                                  + self_faceVertexIDs.shape[1:], 'l'),
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
        vertex_map = numerix.empty(otherNumVertices, dtype=numerix.INT_DTYPE)
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
        if (selfc.numberOfFaces > 0
            and otherc.numberOfFaces > 0
            and faceCorrelates.shape[-1] == 0):
            import warnings
            warnings.warn("Faces are not aligned", UserWarning, stacklevel=4)

        # map other's Face IDs to new Face IDs,
        # accounting for overlaps with self's Face IDs
        face_map = numerix.empty(otherNumFaces, dtype=numerix.INT_DTYPE)
        facesToAdd = numerix.delete(numerix.arange(otherNumFaces), faceCorrelates[1])
        face_map[facesToAdd] = numerix.arange(otherNumFaces - len(faceCorrelates[1])) + selfNumFaces
        face_map[faceCorrelates[1]] = faceCorrelates[0]

        other_faceVertexIDs = vertex_map[otherc.faceVertexIDs[..., facesToAdd]]

        # ensure that both sets of cellFaceIDs have the same maximum number of (masked) elements
        self_cellFaceIDs = selfc.cellFaceIDs
        other_cellFaceIDs = face_map[otherc.cellFaceIDs]
        diff = self_cellFaceIDs.shape[0] - other_cellFaceIDs.shape[0]
        if diff > 0:
            other_cellFaceIDs = numerix.append(other_cellFaceIDs,
                                               -1 * numerix.ones((diff,)
                                                                 + other_cellFaceIDs.shape[1:], 'l'),
                                               axis=0)
            other_cellFaceIDs = MA.masked_values(other_cellFaceIDs, -1)
        elif diff < 0:
            self_cellFaceIDs = numerix.append(self_cellFaceIDs,
                                              -1 * numerix.ones((-diff,)
                                                                + self_cellFaceIDs.shape[1:], 'l'),
                                              axis=0)
            self_cellFaceIDs = MA.masked_values(self_cellFaceIDs, -1)

        # concatenate everything and return
        return {
            'vertexCoords': numerix.concatenate((selfc.vertexCoords,
                                                 otherc.vertexCoords[..., verticesToAdd]), axis=1),
            'faceVertexIDs': numerix.concatenate((self_faceVertexIDs,
                                                  other_faceVertexIDs), axis=1),
            'cellFaceIDs': MA.concatenate((self_cellFaceIDs,
                                           other_cellFaceIDs), axis=1)
            }

    """
    Topology -- maybe should be elsewhere?
    """

    @property
    def interiorFaceIDs(self):
        if not hasattr(self, '_interiorFaceIDs'):
            self._interiorFaceIDs = numerix.nonzero(self.interiorFaces)[0]
        return self._interiorFaceIDs

    @property
    def interiorFaceCellIDs(self):
        if not hasattr(self, '_interiorFaceCellIDs'):
            ## Commented line is better, but doesn't work for zero length arrays
            ##  self.interiorFaceCellIDs = self.getFaceCellIDs()[..., self.getInteriorFaceIDs()]
            self._interiorFaceCellIDs = numerix.take(self.faceCellIDs,
                                                     self.interiorFaceIDs, axis=1)
        return self._interiorFaceCellIDs

    @property
    def _numberOfFacesPerCell(self):
        cellFaceIDs = self.cellFaceIDs
        if type(cellFaceIDs) is type(MA.array(0)):
            ## bug in count returns float values when there is no mask
            return numerix.array(cellFaceIDs.count(axis=0), 'l')
        else:
            return self._maxFacesPerCell * numerix.ones(cellFaceIDs.shape[-1], 'l')

    @property
    def _maxFacesPerCell(self):
        raise NotImplementedError

    @property
    def _numberOfVertices(self):
        if hasattr(self, 'numberOfVertices'):
            return self.numberOfVertices
        else:
            return self.vertexCoords.shape[-1]

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
        return self.topology._globalNonOverlappingCellIDs

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
        return self.topology._globalOverlappingCellIDs

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
        return self.topology._localNonOverlappingCellIDs

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
        return self.topology._localOverlappingCellIDs

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
        return self.topology._globalNonOverlappingFaceIDs

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
        return self.topology._globalOverlappingFaceIDs

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
        return self.topology._localNonOverlappingFaceIDs

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
        return self.topology._localOverlappingFaceIDs

    @property
    def facesLeft(self):
        """
        Return face on left boundary of Grid1D as list with the
        x-axis running from left to right.

            >>> from fipy import Grid2D, Grid3D
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((21, 25),
            ...                        numerix.nonzero(mesh.facesLeft)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesLeft.value # doctest: +PROCESSOR_NOT_0
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)
            >>> print numerix.allequal((9, 13),
            ...                        numerix.nonzero(mesh.facesLeft)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesLeft.value # doctest: +PROCESSOR_NOT_0

        """
        x = self.faceCenters[0]
        return x == _madmin(x)

    @property
    def facesRight(self):
        """
        Return list of faces on right boundary of Grid3D with the
        x-axis running from left to right.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((24, 28),
            ...                        numerix.nonzero(mesh.facesRight)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesRight.value # doctest: +PROCESSOR_NOT_0
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)
            >>> print numerix.allequal((12, 16),
            ...                        numerix.nonzero(mesh.facesRight)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesRight.value # doctest: +PROCESSOR_NOT_0

        """
        x = self.faceCenters[0]
        return x == _madmax(x)

    @property
    def facesBottom(self):
        """
        Return list of faces on bottom boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((12, 13, 14),
            ...                        numerix.nonzero(mesh.facesBottom)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesBottom.value # doctest: +PROCESSOR_NOT_0
            >>> x, y, z = mesh.faceCenters
            >>> print numerix.allequal((12, 13),
            ...                        numerix.nonzero(mesh.facesBottom & (x < 1))[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesBottom.value # doctest: +PROCESSOR_NOT_0

        """
        y = self.faceCenters[1]
        return y == _madmin(y)

    facesDown = facesBottom

    @property
    def facesTop(self):
        """
        Return list of faces on top boundary of Grid3D with the
        y-axis running from bottom to top.

            >>> from fipy import Grid2D, Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((18, 19, 20),
            ...                        numerix.nonzero(mesh.facesTop)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesTop.value # doctest: +PROCESSOR_NOT_0
            >>> mesh = Grid2D(nx = 3, ny = 2, dx = 0.5, dy = 2.)
            >>> print numerix.allequal((6, 7, 8),
            ...                        numerix.nonzero(mesh.facesTop)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesTop.value # doctest: +PROCESSOR_NOT_0

        """
        y = self.faceCenters[1]
        return y == _madmax(y)

    facesUp = facesTop

    @property
    def facesBack(self):
        """
        Return list of faces on back boundary of Grid3D with the
        z-axis running from front to back.

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((6, 7, 8, 9, 10, 11),
            ...                        numerix.nonzero(mesh.facesBack)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesBack.value # doctest: +PROCESSOR_NOT_0

        """
        z = self.faceCenters[2]
        return z == _madmax(z)

    @property
    def facesFront(self):
        """
        Return list of faces on front boundary of Grid3D with the
        z-axis running from front to back.

            >>> from fipy import Grid3D, numerix
            >>> mesh = Grid3D(nx = 3, ny = 2, nz = 1, dx = 0.5, dy = 2., dz = 4.)
            >>> print numerix.allequal((0, 1, 2, 3, 4, 5),
            ...                        numerix.nonzero(mesh.facesFront)[0]) # doctest: +PROCESSOR_0
            True
            >>> ignore = mesh.facesFront.value # doctest: +PROCESSOR_NOT_0

        """
        z = self.faceCenters[2]
        return z == _madmin(z)

    @property
    def _cellVertexIDs(self):
        raise NotImplementedError

    def _calcOrderedCellVertexIDs(self):
        return self._cellVertexIDs

    @property
    def _orderedCellVertexIDs(self):
        if hasattr(self, "_orderedCellVertexIDs_data"):
            return self._orderedCellVertexIDs_data
        else:
            return self._calcOrderedCellVertexIDs()

    @property
    def _cellDistanceNormals(self):
        return self._cellDistanceNormals/ self._cellDistances

    @property
    def _cellAreaProjections(self):
        return self._cellNormals * self._cellAreas

    """
    Special methods
    """

    @property
    def _concatenatedClass(self):
        return self.topology._concatenatedClass

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

        >>> InvalidMesh = threeDBaseMesh + baseMesh # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
        ...
        MeshAdditionError: Dimensions do not match
        """
        if isinstance(other, AbstractMesh):
            return self._concatenatedClass(**self._getAddedMeshValues(other=other))
        else:
            return self._translate(other)

    __radd__ = __add__

    def __mul__(self, other):
        raise NotImplementedError

    __rmul__ = __mul__

    def __sub__(self, other):
        """
        Tests.
        >>> from fipy import *
        >>> m = Grid1D()
        >>> print (m - ((1,))).cellCenters
        [[-0.5]]
        >>> ((1,)) - m
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand type(s) for -: 'tuple' and 'UniformGrid1D'

        """
        if isinstance(other, AbstractMesh):
            raise TypeError, "'-' is unsupported for meshes, use '+'"
        else:
            return self._translate(-numerix.array(other))

    def __truediv__(self, other):
        """
        Tests.
        >>> from fipy import *
        >>> print (Grid1D(nx=1) / 2.).cellCenters
        [[ 0.25]]
        >>> AbstractMesh(communicator=None) / 2.
        Traceback (most recent call last):
        ...
        NotImplementedError

        """
        return self.__mul__(1 / other)

    __div__ = __truediv__

    def __getstate__(self):
        return self.representation.getstate()

    def __setstate__(self, state):
        return state["_RepresentationClass"].setstate(self, state)

    def __repr__(self):
        return self.representation.repr()

    @property
    def _VTKCellType(self):
        try:
            from tvtk.api import tvtk
        except ImportError, e:
            from enthought.tvtk.api import tvtk
        return tvtk.ConvexPointSet().cell_type

    @property
    def VTKCellDataSet(self):
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

        try:
            from tvtk.api import tvtk
        except ImportError, e:
            from enthought.tvtk.api import tvtk
        num = counts.shape[0]

        cps_type = self._VTKCellType
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

    @property
    def VTKFaceDataSet(self):
        """Returns a TVTK `DataSet` representing the face centers of this mesh
        """
        try:
            from tvtk.api import tvtk
        except ImportError, e:
            from enthought.tvtk.api import tvtk

        points = self.faceCenters
        points = self._toVTK3D(numerix.array(points))
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
        if arr.dtype.name is 'bool':
            # VTK can't do bool, and the exception isn't properly
            # thrown back to the user
            arr = arr.astype('int')
        if rank == 0:
            return arr
        else:
            arr = numerix.concatenate((arr,
                                       numerix.zeros((3 - self.dim,)
                                                     + arr.shape[1:], 'l')))
            return arr.swapaxes(-2, -1)

    def _makePeriodic(self):
        raise NotImplementedError

    @property
    def _maxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    @property
    def _facesPerCell(self):
        if numerix.MA.is_masked(self.cellFaceIDs):
            facesPerCell = (~numerix.MA.getmask(self.cellFaceIDs)).sum(axis=0)
        else:
            facesPerCell = numerix.empty((self.numberOfCells,), dtype=numerix.INT_DTYPE)
            facesPerCell[:] = self._maxFacesPerCell
        return facesPerCell

    @property
    def _cellFaceVertices(self):
        return numerix.take(self.faceVertexIDs, self.cellFaceIDs, axis=1)

    @property
    def _unsortedNodesPerFace(self):
        cellFaceVertices = self._cellFaceVertices
        if numerix.MA.is_masked(cellFaceVertices):
            nodesPerFace = (~cellFaceVertices.mask).sum(axis=0)
        else:
            nodesPerFace = numerix.empty(cellFaceVertices.shape[1:], dtype=numerix.INT_DTYPE)
            nodesPerFace[:] = self.faceVertexIDs.shape[0]

        return nodesPerFace

    @property
    def _nodesPerFace(self):
        return numerix.sort(self._unsortedNodesPerFace, axis=0)[::-1]

    @property
    def _faceOrder(self):
        # note: argsort does not work as documented:
        #   Array of indices that sort a along the specified axis. In other words, a[index_array] yields a sorted a.
        # No, it's not.
        # We need both the sorted values and the sort order.
        return numerix.argsort(self._unsortedNodesPerFace, axis=0)[::-1]

    @property
    def _cellTopology(self):
        """return a map of the topology of each cell"""
        return self.topology._cellTopology

    def _makePeriodic(self):
        raise NotImplementedError

    @property
    def aspect2D(self):
        """The physical y:x aspect ratio of a 2D mesh
        """
        if self.dim != 2:
            raise NotImplementedError
        else:
            vertexIDs = self._orderedCellVertexIDs

            xCoords = numerix.take(self.vertexCoords[0], vertexIDs)
            yCoords = numerix.take(self.vertexCoords[1], vertexIDs)

            return float((yCoords.max() - yCoords.min()) / (xCoords.max() - xCoords.min()))


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
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
