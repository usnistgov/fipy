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

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.meshes.common.mesh import Mesh as _CommonMesh

from fipy.meshes.numMesh.cell import Cell

from fipy.tools.dimensions.physicalField import PhysicalField

class MeshAdditionError(Exception):
    pass
    
class Mesh(_CommonMesh):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs):
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

        _CommonMesh.__init__(self)
        
    """Topology methods"""

    @property
    def _concatenatedClass(self):
        return Mesh
        
    def __add__(self, other):
        if(isinstance(other, Mesh)):
            return self._concatenatedClass(**self._getAddedMeshValues(other=other))
        else:
            return self._translate(other)

    __radd__ = __add__
    
    def __mul__(self, factor):
        newCoords = self.vertexCoords * factor
        newmesh = Mesh(vertexCoords=newCoords, 
                       faceVertexIDs=numerix.array(self.faceVertexIDs), 
                       cellFaceIDs=numerix.array(self.cellFaceIDs))
        return newmesh

    __rmul__ = __mul__

    def _connectFaces(self, faces0, faces1):
        """
        
        Merge faces on the same mesh. This is used to create periodic
        meshes. The first list of faces, `faces1`, will be the faces
        that are used to add to the matrix diagonals. The faces in
        `faces2` will not be used. They aren't deleted but their
        adjacent cells are made to point at `faces1`. The list
        `faces2` are not altered, they still remain as members of
        exterior faces.

           >>> from fipy.meshes.numMesh.grid2D import Grid2D
           >>> mesh = Grid2D(nx = 2, ny = 2, dx = 1., dy = 1.)

           >>> from fipy.tools import parallel
           >>> print parallel.procID != 0 or (mesh._getCellFaceIDs() == [[0, 1, 2, 3],
           ...                                                           [7, 8, 10, 11],
           ...                                                           [2, 3, 4, 5],
           ...                                                           [6, 7, 9, 10]]).flatten().all()
           True

           >>> mesh._connectFaces(numerix.nonzero(mesh.getFacesLeft()), numerix.nonzero(mesh.getFacesRight()))

           >>> print parallel.procID != 0 or (mesh._getCellFaceIDs() == [[0, 1, 2, 3],
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
        assert (faces | self.getExteriorFaces() == self.getExteriorFaces()).all()

        ## following assert checks number of faces are equal, normals are opposite and areas are the same
        assert (numerix.take(self.areaProjections, faces0, axis=1) 
                == numerix.take(-self.areaProjections, faces1, axis=1)).all()

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

        ## calculate new cell distances and add them to faces0
        self.cellDistances = self.cellDistances.copy()
        self.cellDistances.put(faces0, (self.faceToCellDistances[0] 
                                        + self.faceToCellDistances[1]).take(faces0))

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
        _CommonMesh._calcTopology(self)

        ## calculate new geometry
        self._calcOrientedFaceNormals()
        self._calcFaceToCellDistanceRatio()
        self._calcCellToCellDistances()
        self._calcScaledGeometry()
        self._calcFaceAspectRatios()
        
    def _getConcatenableMesh(self):
        return self
        
    def _getAddedMeshValues(self, other, resolution=1e-2):
        """Calculate the parameters to define a concatenation of `other` with `self`
        
        :Parameters:
          - `other`: The :class:`~fipy.meshes.numMesh.Mesh` to concatenate with `self`
          - `resolution`: How close vertices have to be (relative to the smallest 
            cell-to-cell distance in either mesh) to be considered the same

        :Returns:
          A `dict` with 3 elements: the new mesh `vertexCoords`, 
          `faceVertexIDs`, and `cellFaceIDs`.
        """
        
        selfc = self._getConcatenableMesh()
        other = other._getConcatenableMesh()

        ## check dimensions
        if self.getDim() != other.getDim():
            raise MeshAdditionError, "Dimensions do not match"
            
        ## compute vertex correlates

        ## only try to match exterior (X) vertices
        self_Xvertices = numerix.unique(selfc._getFaceVertexIDs().filled()[..., selfc.getExteriorFaces()].flatten().getValue())
        other_Xvertices = numerix.unique(other._getFaceVertexIDs().filled()[..., other.getExteriorFaces()].flatten().getValue())

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
        close = (distance < resolution * min(selfc._getCellToCellDistances().min(), 
                                             other._getCellToCellDistances().min())).getValue()
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
        self_matchingFaces = numerix.in1d(self_faceVertexIDs.getValue(), 
                                          vertexCorrelates[0]).reshape(self_faceVertexIDs.shape).all(axis=0).nonzero()[0]

        # want other's Faces for which all faceVertexIDs are in vertexCorrelates
        other_matchingFaces = numerix.in1d(other_faceVertexIDs.getValue(), 
                                           vertexCorrelates[1]).reshape(other_faceVertexIDs.shape).all(axis=0).nonzero()[0]
                                           
        # map other's Vertex IDs to new Vertex IDs, 
        # accounting for overlaps with self's Vertex IDs
        vertex_map = numerix.empty(other._getNumberOfVertices(), dtype=int)
        verticesToAdd = numerix.delete(numerix.arange(other._getNumberOfVertices()), vertexCorrelates[1])
        vertex_map[verticesToAdd] = numerix.arange(other._getNumberOfVertices() - len(vertexCorrelates[1])) + selfc._getNumberOfVertices()
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
            other_faceHash = vertex_map[other_faceVertexIDs[..., other_matchingFaces].getValue()]
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
        face_map = numerix.empty(other._getNumberOfFaces(), dtype=int)
        facesToAdd = numerix.delete(numerix.arange(other._getNumberOfFaces()), faceCorrelates[1])
        face_map[facesToAdd] = numerix.arange(other._getNumberOfFaces() - len(faceCorrelates[1])) + selfc._getNumberOfFaces()
        face_map[faceCorrelates[1]] = faceCorrelates[0]
        
        other_faceVertexIDs = vertex_map[other.faceVertexIDs[..., facesToAdd].getValue()]
        
        # ensure that both sets of cellFaceIDs have the same maximum number of (masked) elements
        self_cellFaceIDs = selfc.cellFaceIDs
        other_cellFaceIDs = face_map[other.cellFaceIDs.getValue()]
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

    def _calcTopology(self):
        self.dim = self.vertexCoords.shape[0]
        if not hasattr(self, "numberOfFaces"):
            self.numberOfFaces = self.faceVertexIDs.shape[-1]
        if not hasattr(self, "numberOfCells"):
            self.numberOfCells = self.cellFaceIDs.shape[-1]
        if not hasattr(self, "globalNumberOfCells"):
            self.globalNumberOfCells = self.numberOfCells
        if not hasattr(self, "globalNumberOfFaces"):
            self.globalNumberOfFaces = self.numberOfFaces
        self._calcFaceCellIDs()
        
        _CommonMesh._calcTopology(self)


    """calc Topology methods"""

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
                               
        self.faceCellIDs = _FaceCellIDsVariable(mesh=self)

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
                               
        self.vertexFaceIDs = _VertexFaceIDsVariable(mesh=self)


    def _calcInteriorAndExteriorFaceIDs(self):
        self.exteriorFaces = self.faceCellIDs[1].getMaskArray()
        self.interiorFaces = ~self.exteriorFaces

    def _calcInteriorAndExteriorCellIDs(self):
        ids = numerix.take(self.faceCellIDs[0], self.getExteriorFaces(), axis=-1).filled().sorted()
        extras = numerix.array([True] * (len(ids) - len(ids[:-1])), dtype=bool)
        self.exteriorCellIDs = ids[(ids[:-1] != ids[1:]).append(extras)]
        
        from fipy.variables.cellVariable import CellVariable
        self.interiorCellIDs = CellVariable(mesh=self, value=numerix.arange(self.numberOfCells)).delete(self.exteriorCellIDs)
    
#         try:
#             import sets
#             self.exteriorCellIDs = sets.Set(self.faceCellIDs[0, self.getExteriorFaces()])
#             self.interiorCellIDs = list(sets.Set(range(self.numberOfCells)) - self.exteriorCellIDs)
#             self.exteriorCellIDs = list(self.exteriorCellIDs)
#         except:
#             self.exteriorCellIDs = self.faceCellIDs[0, self.getExteriorFaces()]
#             tmp = numerix.zeros(self.numberOfCells)
#             numerix.put(tmp, self.exteriorCellIDs, numerix.ones(len(self.exteriorCellIDs)))
#             self.exteriorCellIDs = numerix.nonzero(tmp)            
#             self.interiorCellIDs = numerix.nonzero(numerix.logical_not(tmp))
            
    def _calcCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs, axis=-1)
        self.cellToFaceOrientations = (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        mask = self.faceCellIDs[1].getMask()
        self.adjacentCellIDs = (self.faceCellIDs[0].filled(),
                                (mask * self.faceCellIDs[0].filled(0) 
                                + ~mask * self.faceCellIDs[1].filled(0)))


    def _calcCellToCellIDs(self):    
        self.cellToCellIDs = numerix.take(self.faceCellIDs, self.cellFaceIDs, axis=1)
        self.cellToCellIDs = ((self.cellToFaceOrientations == 1) * self.cellToCellIDs[1] 
                              + (self.cellToFaceOrientations != 1) * self.cellToCellIDs[0])
        
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
        if hasattr(self, 'vertexCoords'):
            return self.vertexCoords
        else:
            return self._createVertices()

    def getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        return self.exteriorFaces
            
    def getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        return self.interiorFaces
        
    def getFaceCellIDs(self):
        return self.faceCellIDs

    def _getMaxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    """Geometry methods"""

    def _calcGeometry(self):
        self._calcFaceCenters()
        _CommonMesh._calcGeometry(self)
        self._calcCellNormals()
        
    """calc geometry methods"""

    def _calcFaceAreas(self):
        faceVertexIDs = self.faceVertexIDs.filled(-1)
        substitute = numerix.repeat(faceVertexIDs[numerix.newaxis, 0], 
                                    faceVertexIDs.shape[0], axis=0)
        mask = self.faceVertexIDs.getMaskArray()
        faceVertexIDs = mask * substitute + ~mask * faceVertexIDs
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        faceOrigins = numerix.repeat(faceVertexCoords[:,0], faceVertexIDs.shape[0], axis=0)
        faceOrigins = numerix.reshape(faceOrigins, MA.shape(faceVertexCoords))
        faceVertexCoords = faceVertexCoords - faceOrigins
        left = range(faceVertexIDs.shape[0])
        right = left[1:] + [left[0]]
        cross = numerix.sum(numerix.cross(faceVertexCoords, 
                                          numerix.take(faceVertexCoords, right, axis=1),
                                          axis=0), 1)
        self.faceAreas = numerix.sqrtDot(cross, cross) / 2.

    def _calcFaceCenters(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
## 
##         self.faceCenters = MA.filled(MA.average(faceVertexCoords, axis=1))
        
        self.faceCenters = faceVertexCoords.mean(axis=1).filled()

        

    def _calcFaceNormals(self):
        faceVertexIDs = self.faceVertexIDs.filled(0)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.cross(t1, t2, axis=0)
        ## reordering norm's internal memory for inlining
        norm = norm.copy()
        norm = norm / numerix.sqrtDot(norm, norm)
        
        self.faceNormals = -norm
        
        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.cellDistanceVectors) < 0)
        self.faceNormals *= orientation

    def _calcFaceCellToCellNormals(self):
        faceCellCentersUp = numerix.take(self.cellCenters, self.getFaceCellIDs()[1], axis=1)
        faceCellCentersDown = numerix.take(self.cellCenters, self.getFaceCellIDs()[0], axis=1)
        mask = faceCellCentersUp.getMaskArray()
        faceCellCentersUp = mask * self.getFaceCenters() + ~mask * faceCellCentersUp.filled()

        diff = faceCellCentersDown - faceCellCentersUp
        mag = numerix.sqrtDot(diff, diff)
        self.faceCellToCellNormals = diff / mag

        orientation = 1 - 2 * (numerix.dot(self.faceNormals, self.faceCellToCellNormals) < 0)
        self.faceCellToCellNormals *= orientation

    def _calcOrientedFaceNormals(self):
        self.orientedFaceNormals = self.faceNormals
        
    def _calcCellVolumes(self):
        tmp = self.faceCenters[0] * self.faceAreas * self.faceNormals[0]
        tmp = numerix.take(tmp, self.cellFaceIDs, axis=-1) * self.cellToFaceOrientations
        self.cellVolumes = tmp.sum(0).filled()
        self.cellVolumes.name = self.__class__.__name__ + ".cellVolumes"

    def _calcCellCenters(self):
        tmp = numerix.take(self.faceCenters, self.cellFaceIDs, axis=1)
        self.cellCenters = tmp.mean(axis=1).filled()
        self.cellCenters.name = self.__class__.__name__ + ".cellCenters"
        
    def _calcFaceToCellDistances(self):
        tmp = self.faceCenters[...,numerix.newaxis,:]
        # array -= masked_array screws up masking for on numpy 1.1

        tmp = tmp - numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        self.cellToFaceDistanceVectors = tmp
        self.faceToCellDistances = (tmp * tmp).sum(axis=0).sqrt()
        self.faceToCellDistances.name = self.__class__.__name__ + ".faceToCellDistances"

    def _calcCellDistances(self):
        tmp = numerix.take(self.cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[...,1,:] - tmp[...,0,:]
        self.cellDistanceVectors = (tmp.getMask() * self.cellToFaceDistanceVectors[:,0].filled() 
                                    + ~tmp.getMask() * tmp.filled())
        self.cellDistances = self.cellDistanceVectors.getMag()
        self.cellDistances.name = self.__class__.__name__ + ".cellDistances"

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._getCellDistances()
        dFP = self._getFaceToCellDistances()[0]
        
        self.faceToCellDistanceRatio = (dFP / dAP).filled()
        self.faceToCellDistanceRatio.name = self.__class__.__name__ + ".faceToCellDistanceRatio"

    def _calcAreaProjections(self):
        self.areaProjections = self._getFaceNormals() * self._getFaceAreas()
        self.areaProjections.name = self.__class__.__name__ + ".areaProjections"
        
    def _calcOrientedAreaProjections(self):
        self.orientedAreaProjections = self.areaProjections

    def _calcFaceTangents(self):
        faceVertexCoord = numerix.take(self.vertexCoords, 
                                       self.faceVertexIDs[0], 
                                       axis=1)
        tmp = (self.faceCenters - faceVertexCoord).filled()
        self.faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)
        tmp = numerix.cross(self.faceTangents1, self.faceNormals, axis=0)
        self.faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)
        self.faceTangents1.name = self.__class__.__name__ + ".faceTangents1"
        self.faceTangents2.name = self.__class__.__name__ + ".faceTangents2"
        
    def _calcCellToCellDistances(self):
        self.cellToCellDistances = numerix.take(self.cellDistances, self._getCellFaceIDs(), axis=-1)
        self.cellToCellDistances.name = self.__class__.__name__ + ".cellToCellDistances"

    def _calcCellNormals(self):
        cellNormals = numerix.take(self._getFaceNormals(), self._getCellFaceIDs(), axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs, axis=-1)
        cellIDs = numerix.repeat(numerix.arange(self.getNumberOfCells())[numerix.newaxis,...], 
                                 self._getMaxFacesPerCell(), 
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        if self._getMaxFacesPerCell() > 0:
            self.cellNormals =  direction[numerix.newaxis, ...] * cellNormals
        else:
            self.cellNormals = cellNormals
        self.cellNormals.name = self.__class__.__name__ + ".cellNormals"
                         
    """get geometry methods"""

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
        cellVertexIDs = numerix.reshape(cellFaceVertices, (-1, self.getNumberOfCells()))
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

## Below is an ordered version of _getCellVertexIDs()
## It works for the test case in this file (other than the ordering, obviously)
## I've left it in as it may be useful when we need ordered vertices for cells
    
##    def _getOrderedCellVertexIDs(self):

##        ## Get all the vertices from all the faces for each cell
##        from fipy.tools.numerix import take
##        cellFaceVertices = take(self.faceVertexIDs, self.cellFaceIDs)

##        ## get a sorted list of vertices for each cell
##        NCells = self.getNumberOfCells()
##        cellVertexIDs = MA.reshape(cellFaceVertices.flat, (NCells, -1))
##        newmask = MA.getmaskarray(cellVertexIDs).copy()

##        for i in range(len(cellVertexIDs[0]) - 1):
##            for j in range(len(cellVertexIDs[0]))[i + 1:]:

##                newmask[:,j] = MA.where(newmask[:,j],
##                                        newmask[:,j],
##                                        MA.where(MA.filled(cellVertexIDs)[:,i] == MA.filled(cellVertexIDs)[:,j],
##                                                 1,
##                                                 newmask[:,j]))


##        cellVertexIDs = MA.masked_array(cellVertexIDs, newmask)

##        for i in range(len(cellVertexIDs[0]) - 1):
##            j = i + 1
##            while j < len(cellVertexIDs[0]):
##                tmp = cellVertexIDs[:]
##                tmp[:, i] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
##                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
##                                              cellVertexIDs[:, i],
##                                              cellVertexIDs[:, j]),
##                                     cellVertexIDs[:, i])
                                                        
##                tmp[:, j] = MA.where(MA.getmaskarray(cellVertexIDs[:,i]),
##                                     MA.where(MA.getmaskarray(cellVertexIDs[:, j]),
##                                              cellVertexIDs[:,j],
##                                              cellVertexIDs[:,i]),
##                                     cellVertexIDs[:, j])

##                cellVertexIDs = tmp[:]

##                j += 1


##        length = len(cellVertexIDs[0]) - min(numerix.sum(MA.getmaskarray(cellVertexIDs), axis = 1))
##        return cellVertexIDs[:, :length]


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

## pickling

##    self.__getinitargs__(self):
##        return (self.vertexCoords, self.faceVertexIDs, self.cellFaceIDs)
    

    def __getstate__(self):
        dict = {
            'vertexCoords' : self.vertexCoords.getValue(),            
            'faceVertexIDs' : self.faceVertexIDs.getValue(),
            'cellFaceIDs' : self.cellFaceIDs.getValue() }
        return dict

    def __setstate__(self, dict):
        Mesh.__init__(self, **dict)
     
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
            ...                        numerix.nonzero(mesh.getExteriorFaces()))
            1

            >>> internalFaces = numerix.array((3,))
            >>> print numerix.allequal(internalFaces, 
            ...                        numerix.nonzero(mesh.getInteriorFaces()))
            1

            >>> from fipy.tools.numerix import MA
            >>> faceCellIds = MA.masked_values((( 0,  0,  0, 0,  0,  0,  1,  1,  1,  1),
            ...                                 (-1, -1, -1, 1, -1, -1, -1, -1, -1, -1)), -1)
            >>> print numerix.allequal(faceCellIds, mesh.getFaceCellIDs())
            1
            
            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> print numerix.allclose(faceAreas, mesh._getFaceAreas(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> print numerix.allclose(faceCenters, mesh.getFaceCenters(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> print numerix.allclose(faceNormals, mesh._getFaceNormals(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> print numerix.allequal(cellToFaceOrientations, mesh._getCellFaceOrientations())
            1
                                             
            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> print numerix.allclose(cellVolumes, mesh.getCellVolumes(), 
            ...                        atol = 1e-10, rtol = 1e-10)
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
            >>> print numerix.allclose(cellDistances, mesh._getCellDistances(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1
            
            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> print numerix.allclose(faceToCellDistanceRatios, 
            ...                        mesh._getFaceToCellDistanceRatio(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> print numerix.allclose(areaProjections, mesh._getAreaProjections(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents1, mesh._getFaceTangents1(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> print numerix.allclose(tangents2, mesh._getFaceTangents2(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1), 
            ...                                   ( 1, -1), 
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> print numerix.allequal(cellToCellIDs, mesh._getCellToCellIDs())
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> print numerix.allclose(cellToCellDistances, 
            ...                        mesh._getCellToCellDistances(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> print numerix.allequal(interiorCellIDs, mesh._getInteriorCellIDs())
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> print numerix.allequal(exteriorCellIDs, mesh._getExteriorCellIDs())
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
            >>> print numerix.allclose(cellNormals, mesh._getCellNormals(), 
            ...                        atol = 1e-10, rtol = 1e-10)
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
            >>> print numerix.allclose(cellAreaProjections, 
            ...                        mesh._getCellAreaProjections(), 
            ...                        atol = 1e-10, rtol = 1e-10)
            1

            >>> cellVertexIDs = MA.masked_values(((7, 6, 5, 4, 3, 2, 1, 0), 
            ...                                   (9, 8, 6, 5, 2, 1, -1000, -1000)), 
            ...                                  -1000)
            >>> cellVertexIDs = MA.masked_values(((7, 9),
            ...                                   (6, 8),
            ...                                   (5, 6),
            ...                                   (4, 5),
            ...                                   (3, 2),
            ...                                   (2, 1),
            ...                                   (1, -1000),
            ...                                   (0, -1000)), -1000)
            >>> print numerix.allclose(cellVertexIDs, mesh._getCellVertexIDs())
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

            >>> from fipy.meshes.numMesh.uniformGrid1D import UniformGrid1D
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
        cvi = self._getOrderedCellVertexIDs().getValue().swapaxes(0,1)
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

        points = self.getVertexCoords()
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
        points = self._toVTK3D(points.getValue())
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


                      
### test test test    

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
