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

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import serial

class MeshAdditionError(Exception):
    pass
    
class Mesh(AbstractMesh):
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
        exteriorFaces = self.faceCellIDs[1].maskArray 
        interiorFaces = ~exteriorFaces
        
        return interiorFaces, exteriorFaces
           
    def _calcInteriorAndExteriorCellIDs(self):
        ids = numerix.take(self.faceCellIDs[0], self.exteriorFaces, axis=-1).filled().sorted() 
        extras = numerix.array([True] * (len(ids) - len(ids[:-1])), dtype=bool) 
        exteriorCellIDs = ids[(ids[:-1] != ids[1:]).append(extras)] 
        from fipy.variables.cellVariable import CellVariable
        interiorCellIDs = CellVariable(mesh=self, 
                                       value=numerix.arange(self.numberOfCells)).delete(exteriorCellIDs)
            
        return interiorCellIDs, exteriorCellIDs
       
    def _calcCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs, axis=-1)
        return (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _calcAdjacentCellIDs(self):
        mask = self.faceCellIDs[1].mask 
        return (self.faceCellIDs[0].filled(), 
                (mask * self.faceCellIDs[0].filled(0)  
                 + ~mask * self.faceCellIDs[1].filled(0))) 

    def _calcCellToCellIDs(self):    
        cellToCellIDs = numerix.take(self.faceCellIDs, self.cellFaceIDs, axis=1)
        cellToCellIDs = ((self._cellToFaceOrientations == 1) * cellToCellIDs[1]  
                         + (self._cellToFaceOrientations != 1) * cellToCellIDs[0]) 
        return cellToCellIDs 
     
    def _calcCellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self._maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        from fipy.variables.cellVariable import CellVariable
        cellIDs = CellVariable(mesh=self, value=cellIDs) 
        mask = self._cellToCellIDs.mask 
        return mask * cellIDs  + ~mask * self._cellToCellIDs.filled() 

    def _isOrthogonal(self):
        return False          

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
        self._faceNormals = self._calcFaceNormals()
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
        faceVertexIDs = self.faceVertexIDs.filled(-1)
        substitute = numerix.repeat(faceVertexIDs[numerix.newaxis, 0], 
                                    faceVertexIDs.shape[0], axis=0)
        mask = self.faceVertexIDs.maskArray
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
        return numerix.sqrtDot(cross, cross) / 2.
    
    def _calcFaceCenters(self):
        faceVertexCoords = numerix.take(self.vertexCoords, self.faceVertexIDs, axis=1)
    
        return faceVertexCoords.mean(axis=1).filled() 
     
    def _calcFaceNormals(self):
        faceVertexIDs = self.faceVertexIDs.filled(0)
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
        mask = faceCellCentersUp.maskArray
        faceCellCentersUp = mask * self.faceCenters + ~mask * faceCellCentersUp.filled()

        diff = faceCellCentersDown - faceCellCentersUp
        mag = numerix.sqrtDot(diff, diff)
        faceCellToCellNormals = diff / mag

        orientation = 1 - 2 * (numerix.dot(self._faceNormals, faceCellToCellNormals) < 0)
        return faceCellToCellNormals * orientation

    def _calcOrientedFaceNormals(self):
        return self._faceNormals
        
    def _calcCellVolumes(self):
        tmp = self._faceCenters[0] * self._faceAreas * self._faceNormals[0]
        tmp = numerix.take(tmp, self.cellFaceIDs, axis=-1) * self._cellToFaceOrientations
        cellVolumes = tmp.sum(0).filled() 
        cellVolumes.name = self.__class__.__name__ + ".cellVolumes" 
        return cellVolumes

    def _calcCellCenters(self):
        tmp = numerix.take(self._faceCenters, self.cellFaceIDs, axis=1)
        cellCenters = tmp.mean(axis=1).filled()
        cellCenters.name = self.__class__.__name__ + ".cellCenters" 
        return cellCenters 
        
    def _calcFaceToCellDistAndVec(self):
        tmp = self._faceCenters[...,numerix.newaxis,:] 
        # array -= masked_array screws up masking for on numpy 1.1

        tmp = tmp - numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        cellToFaceDistanceVectors = tmp
        faceToCellDistances = numerix.sqrt((tmp * tmp).sum(axis=0))
        faceToCellDistances.name = self.__class__.__name__ + ".faceToCellDistances"
        return faceToCellDistances, cellToFaceDistanceVectors

    def _calcCellDistAndVec(self):
        tmp = numerix.take(self._cellCenters, self.faceCellIDs, axis=1)
        tmp = tmp[...,1,:] - tmp[...,0,:]
        cellDistanceVectors = (tmp.mask * self.cellToFaceDistanceVectors[:,0].filled()  
                               + ~tmp.mask * tmp.filled()) 
        cellDistances = cellDistanceVectors.mag 
        cellDistances.name = self.__class__.__name__ + ".cellDistances"
        return cellDistances, cellDistanceVectors

    def _calcFaceTangents(self):
        faceVertexCoord = numerix.take(self.vertexCoords,  
                                       self.faceVertexIDs[0],  
                                       axis=1) 
        tmp = (self.faceCenters - faceVertexCoord).filled()
        faceTangents1 = tmp / numerix.sqrtDot(tmp, tmp)
        tmp = numerix.cross(faceTangents1, self._faceNormals, axis=0)
        faceTangents2 = tmp / numerix.sqrtDot(tmp, tmp)
        faceTangents1.name = self.__class__.__name__ + ".faceTangents1" 
        faceTangents2.name = self.__class__.__name__ + ".faceTangents2" 
        return faceTangents1, faceTangents2
        
    def _calcCellToCellDist(self):
        cellToCellDistances = numerix.take(self._cellDistances, self.cellFaceIDs, axis=-1) 
        cellToCellDistances.name = self.__class__.__name__ + ".cellToCellDistances" 
        return cellToCellDistances

    def _calcCellAreas(self):
        from fipy.tools.numerix import take
        return take(self._faceAreas, self.cellFaceIDs)
    
    def _calcCellNormals(self):
        cellNormals = numerix.take(self._faceNormals, self.cellFaceIDs, axis=1)
        cellFaceCellIDs = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        cellIDs = numerix.repeat(numerix.arange(self.numberOfCells)[numerix.newaxis,...], 
                                 self._maxFacesPerCell,
                                 axis=0)
        direction = (cellFaceCellIDs == cellIDs) * 2 - 1
        if self._maxFacesPerCell > 0:
            cellNormals = direction[numerix.newaxis, ...] * cellNormals 
        cellNormals.name = self.__class__.__name__ + ".cellNormals"       
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
        areaProjections = self._faceNormals * self._faceAreas
        areaProjections.name = self.__class__.__name__ + ".areaProjections"
        return areaProjections
        
    def _calcOrientedAreaProjections(self):
        return self._areaProjections

    def _calcFaceToCellDistanceRatio(self):
        dAP = self._cellDistances
        dFP = self._faceToCellDistances[0]
        
        faceToCellDistanceRatio = (dFP / dAP).filled() 
        faceToCellDistanceRatio.name = self.__class__.__name__ + ".faceToCellDistanceRatio" 
        return faceToCellDistanceRatio
       
    def _calcFaceAspectRatios(self):
        return self._scaledFaceAreas / self._cellDistances
    
    @property
    def _concatenatedClass(self):
        return Mesh

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
    
    @property
    def _concatenableMesh(self):
        return self

    def _translate(self, vector):
        newCoords = self.vertexCoords + vector
        newmesh = Mesh(vertexCoords=newCoords,
                       faceVertexIDs=self.faceVertexIDs,
                       cellFaceIDs=self.cellFaceIDs)
        return newmesh

    def _handleFaceConnection(self):
        self._orientedFaceNormals = self._calcOrientedFaceNormals() 
        self._cellToCellDistances = self._calcCellToCellDist()   
        self._setFaceDependentScaledValues()

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
                                 mask=self.mesh.cellFaceIDs.mask)
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
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)

        tmp = cellVertexIDs[:-1].masked(cellVertexIDs[:-1] == cellVertexIDs[1:])
        cellVertexIDs = tmp.append(cellVertexIDs[-1, numerix.newaxis], axis=0)
        cellVertexIDs = cellVertexIDs.sorted(axis=0, fill_value=-1)
        
        ## resize the array to remove extra masked values
        if cellVertexIDs.shape[-1] == 0:
            length = 0
        else:
            length = cellVertexIDs.mask.sum(axis=0).min()
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

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
