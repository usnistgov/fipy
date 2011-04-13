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

from fipy.variables.cellVariable import CellVariable

from abstractMeshTopology import AbstractMeshTopology

class MeshTopology(AbstractMeshTopology):
    """
    Covers all non-uniform Mesh topologies.
    """

    def __init__(self, cellFaceIDs, faceCellIDs, numCells, maxFacesPerCell, 
                 mesh):
        self.mesh = mesh
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
        exteriorFaces = self.faceCellIDs[1].maskArray
        interiorFaces = ~exteriorFaces

        return interiorFaces, exteriorFaces
           
    def _calcInteriorAndExteriorCellIDs(self, numCells):
        ids = numerix.take(self.faceCellIDs[0], self.exteriorFaces, axis=-1).filled().sorted()
        extras = numerix.array([True] * (len(ids) - len(ids[:-1])), dtype=bool)
        exteriorCellIDs = ids[(ids[:-1] != ids[1:]).append(extras)]
        
        interiorCellIDs = CellVariable(mesh=self.mesh, value=numerix.arange(numCells)).delete(exteriorCellIDs)
    
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
        cellToCellIDs = ((self.cellToFaceOrientations == 1) * cellToCellIDs[1] 
                         + (self.cellToFaceOrientations != 1) * cellToCellIDs[0])
        return cellToCellIDs 
     
    def _calcCellToCellIDsFilled(self, numCells, maxFacesPerCell):
        N = numCells
        M = maxFacesPerCell
        cellIDs = CellVariable(mesh=self.mesh, value=numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0))
        mask = self.cellToCellIDs.mask
        return mask * cellIDs  + ~mask * self._cellToCellIDs.filled()
