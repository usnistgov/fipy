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

from fipy.meshes.abstractMesh import AbstractMesh
from fipy.meshes.representations.meshRepresentation import _MeshRepresentation
from fipy.meshes.topologies.meshTopology import _MeshTopology

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serialComm

__all__ = ["MeshAdditionError", "Mesh"]

class MeshAdditionError(Exception):
    pass

class Mesh(AbstractMesh):
    """Generic mesh class using numerix to do the calculations

        Meshes contain cells, faces, and vertices.

        This is built for a non-mixed element mesh.
    """

    def __init__(self, vertexCoords, faceVertexIDs, cellFaceIDs, communicator=serialComm, _RepresentationClass=_MeshRepresentation, _TopologyClass=_MeshTopology):
        super(Mesh, self).__init__(communicator=communicator,
                                   _RepresentationClass=_RepresentationClass,
                                   _TopologyClass=_TopologyClass)

        """faceVertexIds and cellFacesIds must be padded with minus ones."""

        self.vertexCoords = vertexCoords
        self.faceVertexIDs = MA.masked_values(faceVertexIDs, -1)
        self.cellFaceIDs = MA.masked_values(cellFaceIDs, -1)

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

    """
    Topology set and calc
    """

    def _setTopology(self):
        (self._interiorFaces,
         self._exteriorFaces) = self._calcInteriorAndExteriorFaceIDs()
        (self._interiorCellIDs,
         self._exteriorCellIDs) = self._calcInteriorAndExteriorCellIDs()
        self._cellToFaceOrientations = self._calcCellToFaceOrientations()
        self._adjacentCellIDs = self._calcAdjacentCellIDs()
        self._cellToCellIDs = self._calcCellToCellIDs()
        self._cellToCellIDsFilled = self._calcCellToCellIDsFilled()

    def _calcInteriorAndExteriorFaceIDs(self):
        from fipy.variables.faceVariable import FaceVariable
        mask = MA.getmask(self.faceCellIDs[1])
        exteriorFaces = FaceVariable(mesh=self,
                                     value=mask)
        interiorFaces = FaceVariable(mesh=self,
                                     value=numerix.logical_not(mask))
        return interiorFaces, exteriorFaces

    def _calcInteriorAndExteriorCellIDs(self):
        try:
            import sys
            if sys.version_info < (2, 6):
                from sets import Set as set
            exteriorCellIDs = set(self.faceCellIDs[0, self._exteriorFaces.value])
            interiorCellIDs = list(set(range(self.numberOfCells)) - self._exteriorCellIDs)
            exteriorCellIDs = list(self._exteriorCellIDs)
        except:
            exteriorCellIDs = self.faceCellIDs[0, self._exteriorFaces.value]
            tmp = numerix.zeros(self.numberOfCells, 'l')
            numerix.put(tmp, exteriorCellIDs, numerix.ones(len(exteriorCellIDs), 'l'))
            exteriorCellIDs = numerix.nonzero(tmp)
            interiorCellIDs = numerix.nonzero(numerix.logical_not(tmp))
        return interiorCellIDs, exteriorCellIDs

    def _calcCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        return (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        return (MA.filled(self.faceCellIDs[0]),
                          MA.filled(MA.where(MA.getmaskarray(self.faceCellIDs[1]),
                              self.faceCellIDs[0],
                                             self.faceCellIDs[1])))

    def _calcCellToCellIDs(self):
        cellToCellIDs = numerix.take(self.faceCellIDs, self.cellFaceIDs, axis=1)
        cellToCellIDs = MA.where(self._cellToFaceOrientations == 1,
                                 cellToCellIDs[1], cellToCellIDs[0])
        return cellToCellIDs

    def _calcCellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self._maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        return MA.where(MA.getmaskarray(self._cellToCellIDs), cellIDs,
                        self._cellToCellIDs)

    """
    Geometry set and calc
    """

    def _setGeometry(self, scaleLength = 1.):
        self._faceCenters = self._calcFaceCenters()
        self._faceAreas = self._calcFaceAreas()
        self._cellCenters = self._calcCellCenters()
        (self._internalFaceToCellDistances,
         self._cellToFaceDistanceVectors) = self._calcFaceToCellDistAndVec()
        (self._internalCellDistances,
         self._cellDistanceVectors) = self._calcCellDistAndVec()
        self.faceNormals = self._calcFaceNormals()
        self._orientedFaceNormals = self._calcOrientedFaceNormals()
        self._cellVolumes = self._calcCellVolumes()
        self._faceCellToCellNormals = self._calcFaceCellToCellNormals()
        (self._faceTangents1,
         self._faceTangents2) = self._calcFaceTangents()
        self._cellToCellDistances = self._calcCellToCellDist()

        self._setScaledGeometry(self.scale['length'])

        self._cellAreas = self._calcCellAreas()
        self._cellNormals = self._calcCellNormals()

    def _calcFaceAreas(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, -1)
        substitute = numerix.repeat(faceVertexIDs[numerix.newaxis, 0],
                                    faceVertexIDs.shape[0], axis=0)
        faceVertexIDs = numerix.where(MA.getmaskarray(self.faceVertexIDs),
                                      substitute, faceVertexIDs)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        faceOrigins = numerix.repeat(faceVertexCoords[:,0], faceVertexIDs.shape[0], axis=0)
        faceOrigins = numerix.reshape(faceOrigins, MA.shape(faceVertexCoords))
        faceVertexCoords = faceVertexCoords - faceOrigins
        left = range(faceVertexIDs.shape[0])
        right = left[1:] + [left[0]]
        cross = numerix.sum(numerix.cross(faceVertexCoords,
                                          numerix.take(faceVertexCoords, right, 1),
                                          axis=0),
                            1)
        return numerix.sqrtDot(cross, cross) / 2.

    def _calcFaceCenters(self):
        maskedFaceVertexIDs = MA.filled(self.faceVertexIDs, 0)

        faceVertexCoords = numerix.take(self.vertexCoords, maskedFaceVertexIDs, axis=1)

        if MA.getmask(self.faceVertexIDs) is False:
            faceVertexCoordsMask = numerix.zeros(numerix.shape(faceVertexCoords), 'l')
        else:
            faceVertexCoordsMask = \
              numerix.repeat(MA.getmaskarray(self.faceVertexIDs)[numerix.newaxis,...],
                             self.dim, axis=0)

        faceVertexCoords = MA.array(data=faceVertexCoords, mask=faceVertexCoordsMask)

        return MA.filled(MA.average(faceVertexCoords, axis=1))

    @property
    def _rightHandOrientation(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.cross(t1, t2, axis=0)
        ## reordering norm's internal memory for inlining
        norm = norm.copy()
        norm = norm / numerix.sqrtDot(norm, norm)

        faceNormals = -norm

        return 1 - 2 * (numerix.dot(faceNormals, self.cellDistanceVectors) < 0)

    def _calcFaceNormals(self):
        faceVertexIDs = MA.filled(self.faceVertexIDs, 0)
        faceVertexCoords = numerix.take(self.vertexCoords, faceVertexIDs, axis=1)
        t1 = faceVertexCoords[:,1,:] - faceVertexCoords[:,0,:]
        t2 = faceVertexCoords[:,2,:] - faceVertexCoords[:,1,:]
        norm = numerix.cross(t1, t2, axis=0)
        ## reordering norm's internal memory for inlining
        norm = norm.copy()
        norm = norm / numerix.sqrtDot(norm, norm)

        faceNormals = -norm

        orientation = 1 - 2 * (numerix.dot(faceNormals, self.cellDistanceVectors) < 0)
        return faceNormals * orientation

    def _calcFaceCellToCellNormals(self):
        faceCellCentersUp = numerix.take(self._cellCenters, self.faceCellIDs[1], axis=1)
        faceCellCentersDown = numerix.take(self._cellCenters, self.faceCellIDs[0], axis=1)
        faceCellCentersUp = numerix.where(MA.getmaskarray(faceCellCentersUp),
                                          self._faceCenters,
                                          faceCellCentersUp)

        diff = faceCellCentersDown - faceCellCentersUp
        mag = numerix.sqrt(numerix.sum(diff**2))
        faceCellToCellNormals = diff / numerix.resize(mag, (self.dim, len(mag)))

        orientation = 1 - 2 * (numerix.dot(self.faceNormals, faceCellToCellNormals) < 0)
        return faceCellToCellNormals * orientation

    def _calcOrientedFaceNormals(self):
        return self.faceNormals

    def _calcCellVolumes(self):
        tmp = self._faceCenters[0] * self._faceAreas * self.faceNormals[0]
        tmp = numerix.take(tmp, self.cellFaceIDs) * self._cellToFaceOrientations
        return MA.filled(MA.sum(tmp, 0))

    def _calcCellCenters(self):
        tmp = numerix.take(self._faceCenters, self.cellFaceIDs, axis=1)
        return MA.filled(MA.average(tmp, 1))

    def _calcFaceToCellDistAndVec(self):
        tmp = MA.repeat(self._faceCenters[...,numerix.NewAxis,:], 2, 1)
        # array -= masked_array screws up masking for on numpy 1.1

        tmp = tmp - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        cellToFaceDistanceVectors = tmp
        faceToCellDistances = MA.sqrt(MA.sum(tmp * tmp,0))
        return faceToCellDistances, cellToFaceDistanceVectors

    def _calcCellDistAndVec(self):
        tmp = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[...,1,:] - tmp[...,0,:]
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp), self._cellToFaceDistanceVectors[:,0], tmp))
        cellDistanceVectors = tmp
        cellDistances = MA.filled(MA.sqrt(MA.sum(tmp * tmp, 0)))
        return cellDistances, cellDistanceVectors

    def _calcFaceTangents(self):
        faceVertexCoord = numerix.array(numerix.take(self.vertexCoords,
                                                     self.faceVertexIDs[0],
                                                     axis=1))
        tmp = self._faceCenters - faceVertexCoord
        faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)
        tmp = numerix.cross(faceTangents1, self.faceNormals, axis=0)
        faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)
        return faceTangents1, faceTangents2

    def _calcCellToCellDist(self):
        return numerix.take(self._cellDistances, self.cellFaceIDs)

    def _calcCellAreas(self):
        from fipy.tools.numerix import take
        return take(self._faceAreas, self.cellFaceIDs)

    def _calcCellNormals(self):
        cellNormals = numerix.take(self.faceNormals, self.cellFaceIDs, axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        cellIDs = numerix.repeat(numerix.arange(self.numberOfCells)[numerix.newaxis,...],
                                 self._maxFacesPerCell,
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        if self._maxFacesPerCell > 0:
            return direction[numerix.newaxis, ...] * cellNormals
        else:
            return cellNormals

    """settable geometry properties"""
    def _getFaceToCellDistances(self):
        return self._internalFaceToCellDistances

    def _setFaceToCellDistances(self, v):
        self._internalFaceToCellDistances = v
        self._setScaledValues()

    _faceToCellDistances = property(_getFaceToCellDistances,
                                    _setFaceToCellDistances)

    def _getCellDistances(self):
        return self._internalCellDistances

    def _setCellDistances(self, v):
        self._internalCellDistances = v
        self._setScaledValues()

    _cellDistances = property(_getCellDistances, _setCellDistances)

    """
    Scaled geometry set and calc
    """

    _scale = {
        'length': 1.,
        'area': 1.,
        'volume': 1.,
    }

    def _setScaledGeometry(self, val):
        """
        Set the scale by length.

        :Parameters:
          - `val`: The new scale length.
        """

        self._scale['length'] = PhysicalField(value=val)

        if self._scale['length'].unit.isDimensionless():
            self._scale['length'] = 1

        self._scale['area'] = self._calcAreaScale()
        self._scale['volume'] = self._calcVolumeScale()
        self._setScaledValues()

    def _setScaledValues(self):
        self._scaledFaceAreas = self._scale['area'] * self._faceAreas
        self._scaledCellVolumes = self._scale['volume'] * self._cellVolumes
        self._scaledCellCenters = self._scale['length'] * self._cellCenters
        self._scaledFaceToCellDistances = self._scale['length'] * self._faceToCellDistances
        self._scaledCellDistances = self._scale['length'] * self._cellDistances
        self._setFaceDependentScaledValues()

    def _setFaceDependentScaledValues(self):
        self._scaledCellToCellDistances = self._scale['length'] * self._cellToCellDistances
        self._areaProjections = self._calcAreaProjections()
        self._orientedAreaProjections = self._calcOrientedAreaProjections()
        self._faceToCellDistanceRatio = self._calcFaceToCellDistanceRatio()
        self._faceAspectRatios = self._calcFaceAspectRatios()

    def _calcAreaScale(self):
        return self.scale['length']**2

    def _calcVolumeScale(self):
        return self.scale['length']**3

    def _calcAreaProjections(self):
        return self.faceNormals * self._faceAreas

    def _calcOrientedAreaProjections(self):
        return self._areaProjections

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._cellDistances
        dFP = self._faceToCellDistances[0]

        return MA.filled(dFP / dAP)

    def _calcFaceAspectRatios(self):
        return self._scaledFaceAreas / self._cellDistances

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

            >>> dilatedMesh = baseMesh * ((3,), (2,), (1,)) #doctest: +IGNORE_EXCEPTION_DETAIL
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

    @property
    def _concatenableMesh(self):
        return self

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(newCoords, numerix.array(self.faceVertexIDs), numerix.array(self.cellFaceIDs))
        return newmesh

    def _handleFaceConnection(self):
        """
        The _faceCellToCellNormals were added to ensure faceNormals == _faceCellToCellNormals for periodic grids.

        >>> from fipy import *
        >>> m = PeriodicGrid2DLeftRight(nx=2, ny=2)
        >>> (m.faceNormals == m._faceCellToCellNormals).all()
        True

        """
        self._cellToCellDistances = self._calcCellToCellDist()
        self._faceCellToCellNormals = self._calcFaceCellToCellNormals()
        self._setFaceDependentScaledValues()

    """calc Topology methods"""

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

    """get Topology methods"""

    @property
    def _maxFacesPerCell(self):
        return self.cellFaceIDs.shape[0]

    @property
    def _cellVertexIDs(self):
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
            >>> numerix.allequal(faceCellIds, mesh.faceCellIDs)
            1

            >>> dxdy = dx * dy
            >>> dxdz = dx * dz
            >>> dydz = dy * dz
            >>> faceAreas = numerix.array((dxdy, dxdy, dydz, dydz, dxdz, dxdz,
            ...                            dxdy/2., dxdy/2., dxdz, numerix.sqrt(dx**2 + dy**2) * dz))
            >>> numerix.allclose(faceAreas, mesh._faceAreas, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceCoords = numerix.take(vertices, MA.filled(faces, 0), axis=1)
            >>> faceCenters = faceCoords[...,0,:] + faceCoords[...,1,:] + faceCoords[...,2,:] + faceCoords[...,3,:]
            >>> numVex = numerix.array((4., 4., 4., 4., 4., 4., 3., 3., 4., 4.))
            >>> faceCenters /= numVex
            >>> print numerix.allclose(faceCenters, mesh.faceCenters, atol = 1e-10, rtol = 1e-10)
            True

            >>> faceNormals = numerix.array((( 0., 0., -1., 1.,  0., 0.,  0., 0.,  0., dy / numerix.sqrt(dy**2 + dx**2)),
            ...                              ( 0., 0.,  0., 0., -1., 1.,  0., 0., -1., dx / numerix.sqrt(dy**2 + dx**2)),
            ...                              (-1., 1.,  0., 0.,  0., 0., -1., 1.,  0., 0.)))
            >>> numerix.allclose(faceNormals, mesh.faceNormals, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToFaceOrientations = MA.masked_values(((1, -1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, 1),
            ...                                            (1, -2)), -2)
            >>> numerix.allequal(cellToFaceOrientations, mesh._cellToFaceOrientations)
            1

            >>> cellVolumes = numerix.array((dx*dy*dz, dx*dy*dz / 2.))
            >>> numerix.allclose(cellVolumes, mesh.cellVolumes, atol = 1e-10, rtol = 1e-10)
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
            >>> numerix.allclose(cellDistances, mesh._cellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> faceToCellDistanceRatios = faceToCellDistances[0] / cellDistances
            >>> numerix.allclose(faceToCellDistanceRatios, mesh._faceToCellDistanceRatio, atol = 1e-10, rtol = 1e-10)
            1

            >>> areaProjections = faceNormals * faceAreas
            >>> numerix.allclose(areaProjections, mesh._areaProjections, atol = 1e-10, rtol = 1e-10)
            1

            >>> v1 = numerix.take(vertices, numerix.array(faces[0]), axis=1)
            >>> tmp = faceCenters - v1
            >>> tangents1 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents1, mesh._faceTangents1, atol = 1e-10, rtol = 1e-10)
            1

            >>> tmp = numerix.cross(tangents1, faceNormals, axis=0)
            >>> tangents2 = tmp / numerix.sqrtDot(tmp, tmp)
            >>> numerix.allclose(tangents2, mesh._faceTangents2, atol = 1e-10, rtol = 1e-10)
            1

            >>> cellToCellIDs = MA.masked_values(((-1,  0),
            ...                                   (-1, -1),
            ...                                   (-1, -1),
            ...                                   ( 1, -1),
            ...                                   (-1, -1),
            ...                                   (-1, -1)), -1)
            >>> numerix.allequal(cellToCellIDs, mesh._cellToCellIDs)
            1

            >>> cellToCellDistances = MA.masked_values(((dz / 2., d4),
            ...                                         (dz / 2., -1),
            ...                                         (dx / 2., -1),
            ...                                         (     d4, -1),
            ...                                         (dy / 2., -1),
            ...                                         (dy / 2., -1)), -1)
            >>> numerix.allclose(cellToCellDistances, mesh._cellToCellDistances, atol = 1e-10, rtol = 1e-10)
            1

            >>> interiorCellIDs = numerix.array(())
            >>> numerix.allequal(interiorCellIDs, mesh._interiorCellIDs)
            1

            >>> exteriorCellIDs = numerix.array((0, 1))
            >>> numerix.allequal(exteriorCellIDs, mesh._exteriorCellIDs)
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
            >>> numerix.allclose(cellNormals, mesh._cellNormals, atol = 1e-10, rtol = 1e-10)
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
            >>> numerix.allclose(cellAreaProjections, mesh._cellAreaProjections, atol = 1e-10, rtol = 1e-10)
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
            >>> numerix.allclose(cellVertexIDs, mesh._cellVertexIDs)
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

            >>> c = UniformGrid1D(nx=10) + (UniformGrid1D(nx=10) + 10) # doctest: +SERIAL
            >>> print numerix.allclose(c.cellCenters[0],
            ...                        [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5,
            ...                        12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5, 19.5])
            ... # doctest: +SERIAL
            True

        """

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
