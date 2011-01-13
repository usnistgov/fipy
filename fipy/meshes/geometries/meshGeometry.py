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

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField

from abstractGeometries import AbstractScaledMeshGeometry
from abstractGeometries import AbstractMeshGeometry

class ScaledMeshGeometry(AbstractScaledMeshGeometry):

    """
    Only internal attributes of `geom` are accessed here because the external
    counterparts are actually properties which, in some cases, refer back to the
    scaled attribute being defined (causing a circular dependency).

    TODO: Fix that. Accept attributes as params instead of coupling this class
    with geometry.
    """

    def __init__(self, scaleLength, meshGeom):
        self._scale = {
            'length': 1.,
            'area': 1.,
            'volume':1.
        }
        self.geom = meshGeom
        self._setScaleAndRecalculate(scaleLength)

    def _setScaledValues(self):
        self._scaledFaceAreas = self._scale['area'] * self.geom.faceAreas
        self._scaledCellVolumes = self._scale['volume'] * self.geom.cellVolumes
        self._scaledCellCenters = self._scale['length'] * self.geom.cellCenters
        self._scaledFaceToCellDistances = self._scale['length'] * self.geom.faceToCellDistances
        self._scaledCellDistances = self._scale['length'] * self.geom.cellDistances
        self._scaledCellToCellDistances = self._scale['length'] * self.geom.cellToCellDistances
        self._areaProjections = self._calcAreaProjections()
        self._orientedAreaProjections = self._calcOrientedAreaProjections()
        self._faceToCellDistanceRatio = self._calcFaceToCellDistanceRatio()
        self._faceAspectRatios = self._calcFaceAspectRatios()
         
    def _setScaleAndRecalculate(self, val):
        """
        Set the scale by length.

        :Parameters:
          - `val`: The new scale length.
        """
        self._scale['length'] = PhysicalField(value=val)
        
        if self._scale['length'].getUnit().isDimensionless():
            self._scale['length'] = 1    

        self._scale['area'] = self._calcAreaScale()
        self._scale['volume'] = self._calcVolumeScale()
        self._setScaledValues()

    def _calcAreaScale(self):
        return self.scale['length']**2

    def _calcVolumeScale(self):
        return self.scale['length']**3  
         

    """TODO: Don't like the scale ambiguity"""
    scale                     = property(lambda self: self._scale,
                                         lambda v: _setScaleAndRecalculate(v))
    scaledFaceAreas           = property(lambda self: self._scaledFaceAreas)
    scaledCellVolumes         = property(lambda self: self._scaledCellVolumes)
    scaledCellCenters         = property(lambda self: self._scaledCellCenters)
    scaledFaceToCellDistances = property(lambda self: self._scaledFaceToCellDistances)
    scaledCellDistances       = property(lambda self: self._scaledCellDistances)
    scaledCellToCellDistances = property(lambda self: self._scaledCellToCellDistances)
    areaProjections           = property(lambda self: self._areaProjections)
    orientedAreaProjections   = property(lambda self: self._orientedAreaProjections)
    faceToCellDistanceRatio   = property(lambda self: self._faceToCellDistanceRatio)
    faceAspectRatios          = property(lambda self: self._faceAspectRatios)
    
    def _calcAreaProjections(self):
        return self.geom.faceNormals * self.geom.faceAreas
        
    def _calcOrientedAreaProjections(self):
        return self.areaProjections

    def _calcFaceToCellDistanceRatio(self):
        dAP = self.geom._cellDistances
        dFP = self.geom._faceToCellDistances[0]
        
        return MA.filled(dFP / dAP)
       
    def _calcFaceAspectRatios(self):
        return self.scaledFaceAreas / self.geom._cellDistances
    
class MeshGeometry(AbstractMeshGeometry):

    def __init__(self, dim, 
                       faceVertexIDs,
                       vertexCoords,
                       faceCellIDs,
                       cellFaceIDs,
                       numberOfCells,
                       maxFacesPerCell,
                       cellToFaceOrientations,
                       scaleLength, ScaledGeom=ScaledMeshGeometry):

        self.dim = dim
        self.faceVertexIDs = faceVertexIDs
        self.vertexCoords = vertexCoords
        self.faceCellIDs = faceCellIDs
        self.cellFaceIDs = cellFaceIDs
        self.numberOfCells = numberOfCells
        self.maxFacesPerCell = maxFacesPerCell
        self.cellToFaceOrientations = cellToFaceOrientations

        self._faceAreas = self._calcFaceAreas()
        self._faceCenters = self._calcFaceCenters()
        self._cellCenters = self._calcCellCenters()
        (self._faceToCellDistances,
        self._cellToFaceDistanceVectors) = self._calcFaceToCellDistAndVec()
        (self._cellDistances,
        self._cellDistanceVectors) = self._calcCellDistAndVec()
        self._faceNormals = self._calcFaceNormals()
        self._orientedFaceNormals = self._calcOrientedFaceNormals()
        self._cellVolumes = self._calcCellVolumes()
        self._faceCellToCellNormals = self._calcFaceCellToCellNormals()
        (self._faceTangents1,
        self._faceTangents2) = self._calcFaceTangents()
        self._cellToCellDistances = self._calcCellToCellDist()

        self._scaledGeometry = ScaledGeom(scaleLength, self)

        self._cellAreas = self._calcCellAreas()
        self._cellNormals = self._calcCellNormals()

    def _getCellFaceIDs(self):
        return self._cellFaceIDs

    def _setCellFaceIDs(self, newVal):
        self._cellFaceIDs = newVal

    cellFaceIDs = property(_getCellFaceIDs, _setCellFaceIDs)
        
    def handleFaceConnection(self):
        self._cellToCellDistances = self._calcCellToCellDist()
     
    """internal `_calc*` methods"""

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
            faceVertexCoordsMask = numerix.zeros(numerix.shape(faceVertexCoords))
        else:
            faceVertexCoordsMask = \
              numerix.repeat(MA.getmaskarray(self.faceVertexIDs)[numerix.newaxis,...], 
                             self.dim, axis=0)
            
        faceVertexCoords = MA.array(data=faceVertexCoords, mask=faceVertexCoordsMask)

        return MA.filled(MA.average(faceVertexCoords, axis=1))
     
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
        tmp = numerix.take(tmp, self.cellFaceIDs) * self.cellToFaceOrientations
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
        tmp = MA.filled(MA.where(MA.getmaskarray(tmp), self.cellToFaceDistanceVectors[:,0], tmp))
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
        return take(self.faceAreas, self.cellFaceIDs)
    
    def _calcCellNormals(self):
        cellNormals = numerix.take(self.faceNormals, self.cellFaceIDs, axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        cellIDs = numerix.repeat(numerix.arange(self.numberOfCells)[numerix.newaxis,...], 
                                 self.maxFacesPerCell,
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        if self.maxFacesPerCell > 0:
            return direction[numerix.newaxis, ...] * cellNormals
        else:
            return cellNormals   
      
    """geometry properties"""
    @property
    def faceAreas(self):
        return self._faceAreas

    def _getFaceToCellDistances(self):
        return self._faceToCellDistances

    def _setFaceToCellDistances(self, v):
        self._faceToCellDistances = v
        self._scaledGeometry._setScaledValues()

    faceToCellDistances = property(_getFaceToCellDistances,
                                   _setFaceToCellDistances)

    def _getCellDistances(self):
        return self._cellDistances

    def _setCellDistances(self, v):
        self._cellDistances = v
        self._scaledGeometry._setScaledValues()

    cellDistances = property(_getCellDistances, _setCellDistances)

    def _getFaceNormals(self):
        return self._faceNormals

    def _setFaceNormals(self, v):
        self._faceNormals = v

    faceNormals = property(_getFaceNormals, _setFaceNormals)

    faceCenters               = property(lambda self: self._faceCenters)
    cellToFaceDistanceVectors = property(lambda self: self._cellToFaceDistanceVectors)
    cellDistanceVectors       = property(lambda self: self._cellDistanceVectors)
    orientedFaceNormals       = property(lambda self: self._orientedFaceNormals)


    @property
    def cellVolumes(self):
        return self._cellVolumes

    @property
    def cellCenters(self):
        return self._cellCenters

    faceCellToCellNormals     = property(lambda self: self._faceCellToCellNormals)
    faceTangents1             = property(lambda self: self._faceTangents1)
    faceTangents2             = property(lambda self: self._faceTangents2)
    cellToCellDistances       = property(lambda self: self._cellToCellDistances)
    cellAreas                 = property(lambda self: self._cellAreas)
    cellNormals               = property(lambda self: self._cellNormals)

    """scaled geometery properties"""
    scale                     = property(lambda self: self._scaledGeometry.scale,
                                         lambda s,v: s._scaledGeometry._setScaleAndRecalculate(v))
    scaledFaceAreas           = property(lambda self: self._scaledGeometry._scaledFaceAreas)
    scaledCellVolumes         = property(lambda self: self._scaledGeometry._scaledCellVolumes)

    @property
    def scaledCellCenters(self):
        return self._scaledGeometry.scaledCellCenters

    scaledFaceToCellDistances = property(lambda self: \
                                         self._scaledGeometry.scaledFaceToCellDistances)
    scaledCellDistances       = property(lambda self: \
                                         self._scaledGeometry.scaledCellDistances)
    scaledCellToCellDistances = property(lambda self: \
                                         self._scaledGeometry.scaledCellToCellDistances)
    areaProjections           = property(lambda self: \
                                         self._scaledGeometry.areaProjections)
    orientedAreaProjections   = property(lambda self: \
                                         self._scaledGeometry.orientedAreaProjections)
    faceToCellDistanceRatio   = property(lambda self: \
                                         self._scaledGeometry.faceToCellDistanceRatio)
    faceAspectRatios          = property(lambda self: \
                                          self._scaledGeometry.faceAspectRatios)      

