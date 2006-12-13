#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
 #                                    created: 11/10/03 {2:44:42 PM} 
 #                                last update: 5/15/06 {3:51:57 PM} 
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

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.tools.numerix import MA

from fipy.tools.dimensions.physicalField import PhysicalField

class Mesh:
    """
    Generic mesh class defining implementation-agnostic behavior.

    Make changes to mesh here first, then implement specific implementations in
    `pyMesh` and `numMesh`.

    Meshes contain cells, faces, and vertices.
    """

    def __init__(self):
	self.scale = {
	    'length': 1.,
	    'area': 1.,
	    'volume': 1.
	}
	
	self._calcTopology()
	self._calcGeometry()
    
    def __add__(self, other):
        """
        Either translate a `Mesh` or concatenate two `Mesh` objects.
        
            >>> from fipy.meshes.grid2D import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.getCellCenters()
            [[ 0.5, 0.5,]
             [ 1.5, 0.5,]
             [ 0.5, 1.5,]
             [ 1.5, 1.5,]] 1
             
        If a vector is added to a `Mesh`, a translated `Mesh` is returned
        
            >>> translatedMesh = baseMesh + (5, 10)
            >>> translatedMesh.getCellCenters()
            [[  5.5, 10.5,]
             [  6.5, 10.5,]
             [  5.5, 11.5,]
             [  6.5, 11.5,]]
             
        If a `Mesh` is added to a `Mesh`, a concatenation of the two 
        `Mesh` objects is returned
        
            >>> addedMesh = baseMesh + (baseMesh + (2, 0))
            >>> addedMesh.getCellCenters()
            [[ 0.5, 0.5,]
             [ 1.5, 0.5,]
             [ 0.5, 1.5,]
             [ 1.5, 1.5,]
             [ 2.5, 0.5,]
             [ 3.5, 0.5,]
             [ 2.5, 1.5,]
             [ 3.5, 1.5,]]
        
        The two `Mesh` objects must be properly aligned in order to concatenate them
        
            >>> addedMesh = baseMesh + (baseMesh + (3, 0))
            Traceback (most recent call last):
            ...
            MeshAdditionError: Vertices are not aligned

            >>> addedMesh = baseMesh + (baseMesh + (2, 2))
            Traceback (most recent call last):
            ...
            MeshAdditionError: Faces are not aligned

        No provision is made to avoid or consolidate overlapping `Mesh` objects
        
            >>> addedMesh = baseMesh + (baseMesh + (1, 0))
            >>> addedMesh.getCellCenters()
            [[ 0.5, 0.5,]
             [ 1.5, 0.5,]
             [ 0.5, 1.5,]
             [ 1.5, 1.5,]
             [ 1.5, 0.5,]
             [ 2.5, 0.5,]
             [ 1.5, 1.5,]
             [ 2.5, 1.5,]]
             
        Different `Mesh` classes can be concatenated
         
            >>> from fipy.meshes.tri2D import Tri2D
            >>> triMesh = Tri2D(dx = 1.0, dy = 1.0, nx = 2, ny = 1)
            >>> triMesh = triMesh + (2, 0)
            >>> triAddedMesh = baseMesh + triMesh
            >>> triAddedMesh.getCellCenters()
            [[ 0.5       , 0.5       ,]
             [ 1.5       , 0.5       ,]
             [ 0.5       , 1.5       ,]
             [ 1.5       , 1.5       ,]
             [ 2.83333333, 0.5       ,]
             [ 3.83333333, 0.5       ,]
             [ 2.5       , 0.83333333,]
             [ 3.5       , 0.83333333,]
             [ 2.16666667, 0.5       ,]
             [ 3.16666667, 0.5       ,]
             [ 2.5       , 0.16666667,]
             [ 3.5       , 0.16666667,]]

        but their faces must still align properly
        
            >>> triMesh = Tri2D(dx = 1.0, dy = 2.0, nx = 2, ny = 1)
            >>> triMesh = triMesh + (2, 0)
            >>> triAddedMesh = baseMesh + triMesh
            Traceback (most recent call last):
            ...
            MeshAdditionError: Faces are not aligned

        `Mesh` concatenation is not limited to 2D meshes
        
            >>> from fipy.meshes.grid3D import Grid3D
            >>> threeDBaseMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                         nx = 2, ny = 2, nz = 2)
            >>> threeDSecondMesh = Grid3D(dx = 1.0, dy = 1.0, dz = 1.0, 
            ...                           nx = 1, ny = 1, nz = 1)
            >>> threeDAddedMesh = threeDBaseMesh + (threeDSecondMesh + (2, 0, 0))
            >>> threeDAddedMesh.getCellCenters()
            [[ 0.5, 0.5, 0.5,]
             [ 1.5, 0.5, 0.5,]
             [ 0.5, 1.5, 0.5,]
             [ 1.5, 1.5, 0.5,]
             [ 0.5, 0.5, 1.5,]
             [ 1.5, 0.5, 1.5,]
             [ 0.5, 1.5, 1.5,]
             [ 1.5, 1.5, 1.5,]
             [ 2.5, 0.5, 0.5,]]

        but the different `Mesh` objects must, of course, have the same 
        dimensionality.
        
            >>> InvalidMesh = threeDBaseMesh + baseMesh
            Traceback (most recent call last):
            ...
            MeshAdditionError: Dimensions do not match
        """
        pass
        
    def __mul__(self, factor):
        """
        Dilate a `Mesh` by `factor`.
        
            >>> from fipy.meshes.grid2D import Grid2D
            >>> baseMesh = Grid2D(dx = 1.0, dy = 1.0, nx = 2, ny = 2)
            >>> print baseMesh.getCellCenters()
            [[ 0.5, 0.5,]
             [ 1.5, 0.5,]
             [ 0.5, 1.5,]
             [ 1.5, 1.5,]] 1
             
        The `factor` can be a scalar
        
            >>> dilatedMesh = baseMesh * 3
            >>> dilatedMesh.getCellCenters()
            [[ 1.5, 1.5,]
             [ 4.5, 1.5,]
             [ 1.5, 4.5,]
             [ 4.5, 4.5,]]
             
        or a vector
        
            >>> dilatedMesh = baseMesh * (3, 2)
            >>> dilatedMesh.getCellCenters()
            [[ 1.5, 1. ,]
             [ 4.5, 1. ,]
             [ 1.5, 3. ,]
             [ 4.5, 3. ,]]
        
        but the vector must have the same dimensionality as the `Mesh`
        
            >>> dilatedMesh = baseMesh * (3, 2, 1)
            Traceback (most recent call last):
            ...
            ValueError: frames are not aligned

        """
        pass
        
    def __repr__(self):
        return "%s()" % self.__class__.__name__
        
    """topology methods"""
    
    def _calcTopology(self):
	self._calcInteriorAndExteriorFaceIDs()
	self._calcInteriorAndExteriorCellIDs()
	self._calcCellToFaceOrientations()
	self._calcAdjacentCellIDs()
	self._calcCellToCellIDs()
        self._calcCellToCellIDsFilled()
       
    """calc topology methods"""
	
    def _calcInteriorAndExteriorFaceIDs(self):
	pass

    def _calcExteriorCellIDs(self):
	pass
	
    def _calcInteriorCellIDs(self):
        pass
##	self.interiorCellIDs = list(sets.Set(range(self.numberOfCells)) - sets.Set(self.exteriorCellIDs))
##        onesWhereInterior = numerix.zeros(self.numberOfCells)
##        numerix.put(onesWhereInterior, self.exteriorCells, numerix.zeros((len(self.exteriorCellIDs))))
##        self.interiorCellIDs = numerix.nonzero(onesWhereInterior)
##        self.interiorCellIDs = (0,0)
        
    def _calcInteriorAndExteriorCellIDs(self):
	self._calcExteriorCellIDs()
	self._calcInteriorCellIDs()

    def _calcCellToFaceOrientations(self):
	pass

    def _calcAdjacentCellIDs(self):
	pass

    def _calcCellToCellIDs(self):
	pass

    def _calcCellToCellIDsFilled(self):
        N = self.getNumberOfCells()
        M = self._getMaxFacesPerCell()
        cellIDs = numerix.reshape(numerix.repeat(numerix.arange(N), M), (N, M))
        cellToCellIDs = self._getCellToCellIDs()
        self.cellToCellIDsFilled = MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)

    
    """get topology methods"""

    def _getFaceVertexIDs(self):
        return self.faceVertexIDs

    def _getCellFaceIDs(self):
        return self.cellFaceIDs

    def getExteriorFaces(self):
	pass

    def getInteriorFaces(self):
        pass
	
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
    
    def _getNumberOfVertices(self):
        return len(self.vertexCoords[:,0])
	
    def _getAdjacentCellIDs(self):
        return self.adjacentCellIDs

    def getDim(self):
        return self.dim

    def _getCellsByID(self, ids = None):
	pass
	    
    def getCells(self, filter = None, ids = None, **args):
	"""Return `Cell` objects of `Mesh`."""
	cells = self._getCellsByID(ids)
	
        if filter is not None:
            cells = [cell for cell in cells if filter(cell, **args)]

        return cells
        
    def _getFaces(self):
        pass
    
    def getFaces(self, filter = None, **args):
	"""Return `Face` objects of `Mesh`."""
	faces = self._getFaces()
	
        if filter is not None:
	    return [face for face in faces if filter(face, **args)]

        return faces

    def _getMaxFacesPerCell(self):
	pass

    def _getNumberOfFaces(self):
	return self.numberOfFaces

    def _getCellToCellIDs(self):
        return self.cellToCellIDs

    def _getCellToCellIDsFilled(self):
        return self.cellToCellIDsFilled
	
    """geometry methods"""
    
    def _calcGeometry(self):
	self._calcFaceAreas()
	self._calcFaceNormals()
	self._calcOrientedFaceNormals()
	self._calcCellVolumes()
	self._calcCellCenters()
	self._calcFaceToCellDistances()
	self._calcCellDistances()        
	self._calcFaceTangents()
	self._calcCellToCellDistances()
	self._calcScaledGeometry()
        self._calcCellAreas()
       
    """calc geometry methods"""
    
    def _calcFaceAreas(self):
	pass
	
    def _calcFaceNormals(self):
	pass
	
    def _calcOrientedFaceNormals(self):
	pass
	
    def _calcCellVolumes(self):
	pass
	
    def _calcCellCenters(self):
	pass
	
    def _calcFaceToCellDistances(self):
	pass

    def _calcCellDistances(self):
	pass
        
    def _calcAreaProjections(self):
	pass

    def _calcOrientedAreaProjections(self):
	pass

    def _calcFaceTangents(self):
	pass

    def _calcFaceToCellDistanceRatio(self):
	pass

    def _calcFaceAspectRatios(self):
	self.faceAspectRatios = self._getFaceAreas() / self._getCellDistances()

    def _calcCellToCellDistances(self):
	pass

    def _calcCellAreas(self):
        from fipy.tools.numerix import take
        self.cellAreas =  take(self._getFaceAreas(), self.cellFaceIDs)
    
    """get geometry methods"""
        
    def _getFaceAreas(self):
        return self.scaledFaceAreas

    def _getFaceNormals(self):
        return self.faceNormals
	
    def getCellVolumes(self):
	return self.scaledCellVolumes

    def getCellCenters(self):
	return self.scaledCellCenters

    def _getFaceToCellDistances(self):
        return self.scaledFaceToCellDistances

    def _getCellDistances(self):
        return self.scaledCellDistances

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
	return self.scaledCellToCellDistances

    def _getCellNormals(self):
        return self.cellNormals

    def _getCellAreas(self):
        return self.cellAreas

    def _getCellAreaProjections(self):
        return self.cellNormals * self._getCellAreas()[..., numerix.NewAxis]

    """scaling"""

    def setScale(self, value = 1.):
        self.scale['length'] = PhysicalField(value = value)
        if self.scale['length'].getUnit().isDimensionless():
            self.scale['length'] = 1
	self._calcHigherOrderScalings()
	self._calcScaledGeometry()

    def _calcHigherOrderScalings(self):
	self.scale['area'] = self.scale['length']**2
	self.scale['volume'] = self.scale['length']**3

    def _calcScaledGeometry(self):
	self.scaledFaceAreas = self.scale['area'] * self.faceAreas
	self.scaledCellVolumes = self.scale['volume'] * self.cellVolumes
	self.scaledCellCenters = self.scale['length'] * self.cellCenters
	
	self.scaledFaceToCellDistances = self.scale['length'] * self.faceToCellDistances
	self.scaledCellDistances = self.scale['length'] * self.cellDistances
	self.scaledCellToCellDistances = self.scale['length'] * self.cellToCellDistances
	
	self._calcAreaProjections()
	self._calcOrientedAreaProjections()
	self._calcFaceToCellDistanceRatio()
	self._calcFaceAspectRatios()
	
    """point to cell distances"""
    
    def _getPointToCellDistances(self, point):
	tmp = self.getCellCenters() - PhysicalField(point)
	from fipy.tools import numerix
	return numerix.sqrtDot(tmp, tmp)

    def getNearestCell(self, point):
	return self._getCellsByID([self._getNearestCellID(point)])[0]

    def _getNearestCellID(self, point):
        try:
            tmp = self.getCellCenters() - point
        except TypeError:
            tmp = self.getCellCenters() - PhysicalField(point)
        i = numerix.argmin(numerix.add.reduce((tmp * tmp), axis = 1))
	return i    

## pickling

##    self.__getinitargs__(self):
##        return (self.vertexCoords, self.faceVertexIDs, self.cellFaceIDs)
    

## ##     def __getstate__(self):
## ##         dict = {
## ##             'vertexCoords' : self.vertexCoords,            
## ##             'faceVertexIDs' : self.faceVertexIDs,
## ##             'cellFaceIDs' : self.cellFaceIDs }
## ##         return dict
## ## 
## ##     def __setstate__(self, dict):
## ##         self.__init__(dict['vertexCoords'], dict['faceVertexIDs'], dict['cellFaceIDs'])
        
                      
    
def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
