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

"""TODO: prettify"""

class UniformMeshTopology1D(AbstractMeshTopology):

    def __init__(self, facesLeft, facesRight, numberOfCells, numberOfFaces, 
                       mesh):
        self._exteriorFaces = facesLeft | facesRight
        self.numberOfCells = numberOfCells
        self.numberOfFaces = numberOfFaces
        self.mesh = mesh

    def _getInteriorFaces(self):
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[numerix.arange(self.numberOfFaces-2) + 1] = True
        return interiorFaces
            
    def _getCellToFaceOrientations(self):
        orientations = numerix.ones((2, self.numberOfCells))
        if self.numberOfCells > 0:
            orientations[0] *= -1
            orientations[0,0] = 1
        return orientations

    def _getAdjacentCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = numerix.array((c1 - 1, c1))
        if self.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,-1] = ids[0,-1]
        return ids[0], ids[1]

    def _getCellToCellIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        ids = MA.array((c1 - 1, c1 + 1))
        if self.numberOfCells > 0:
            ids[0,0] = MA.masked
            ids[1,-1] = MA.masked
        return ids
        
    def _getCellToCellIDsFilled(self):
        ids = self.cellToCellIDs.filled()
        if self.numberOfCells > 0:
            ids[0,0] = 0
            ids[1,-1] = self.numberOfCells - 1
        return ids          

    def _getExteriorFaces(self):
        return self._exteriorFaces

    def _setExteriorFaces(self, e):
        self.exteriorFaces = e
    
    """Properties conforming to the MeshTopology interface."""
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces, _setExteriorFaces)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
                                                                              
