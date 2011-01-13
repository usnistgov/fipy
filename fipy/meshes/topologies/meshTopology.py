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

from abstractMeshTopology import AbstractMeshTopology

class MeshTopology(AbstractMeshTopology):
    def __init__(self, cellFaceIDs, faceCellIDs, numCells, maxFacesPerCell, 
                 mesh):
        self.cellFaceIDs = cellFaceIDs
        self.faceCellIDs = faceCellIDs

        (self._interiorFaces,
        self._exteriorFaces) = self._calcInteriorAndExteriorFaceIDs(mesh)
        (self._interiorCellIDs,
        self._exteriorCellIDs) = self._calcInteriorAndExteriorCellIDs(numCells)
        self._cellToFaceOrientations = self._calcCellToFaceOrientations()
        self._adjacentCellIDs = self._calcAdjacentCellIDs()
        self._cellToCellIDs = self._calcCellToCellIDs()
        self._cellToCellIDsFilled = self._calcCellToCellIDsFilled(numCells,
                                                                 maxFacesPerCell)

    def _getCellFaceIDs(self):
        return self._cellFaceIDs

    def _setCellFaceIDs(self, newVal):
        self._cellFaceIDs = newVal

    cellFaceIDs = property(_getCellFaceIDs, _setCellFaceIDs)
    
    def _getInteriorFaces(self):
        return self._interiorFaces

    def _getExteriorFaces(self):
        return self._exteriorFaces

    def _getInteriorCellIDs(self):
        return self._interiorCellIDs

    def _getExteriorCellIDs(self):
        return self._exteriorCellIDs

    def _getCellToFaceOrientations(self):
        return self._cellToFaceOrientations

    def _getAdjacentCellIDs(self):
        return self._adjacentCellIDs

    def _getCellToCellIDs(self):
        return self._cellToCellIDs

    def _getCellToCellIDsFilled(self):
        return self._cellToCellIDsFilled
    
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces)
    interiorCellIDs = property(_getInteriorCellIDs)
    exteriorCellIDs = property(_getExteriorCellIDs)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
      
    def _calcInteriorAndExteriorFaceIDs(self, mesh):
        from fipy.variables.faceVariable import FaceVariable
        mask = MA.getmask(self.faceCellIDs[1])
        exteriorFaces = FaceVariable(mesh=mesh, 
                                     value=mask)
        interiorFaces = FaceVariable(mesh=mesh, 
                                     value=numerix.logical_not(mask))
        return interiorFaces, exteriorFaces
           
    def _calcInteriorAndExteriorCellIDs(self, numCells):
        try:
            import sets
            exteriorCellIDs = sets.Set(self.faceCellIDs[0, self.exteriorFaces.getValue()])
            interiorCellIDs = list(sets.Set(range(numCells)) - self.exteriorCellIDs)
            exteriorCellIDs = list(self.exteriorCellIDs)
        except:
            exteriorCellIDs = self.faceCellIDs[0, self.exteriorFaces.getValue()]
            tmp = numerix.zeros(numCells)
            numerix.put(tmp, exteriorCellIDs, numerix.ones(len(exteriorCellIDs)))
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
        cellToCellIDs = MA.where(self.cellToFaceOrientations == 1, 
                                 cellToCellIDs[1], cellToCellIDs[0])
        return cellToCellIDs 
     
    def _calcCellToCellIDsFilled(self, numCells, maxFacesPerCell):
        N = numCells
        M = maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        return MA.where(MA.getmaskarray(self.cellToCellIDs), cellIDs, 
                        self.cellToCellIDs)
     
