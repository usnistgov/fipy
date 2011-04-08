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

class UniformMeshTopology3D(AbstractMeshTopology):

    def __init__(self, nx, ny, nz,
                       numberOfCells,
                       maxFacesPerCell,
                       XYFaceIDs,
                       XZFaceIDs,
                       YZFaceIDs,
                       faceCellIDs,
                       cellFaceIDs,
                       mesh):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.numberOfCells = numberOfCells
        self.maxFacesPerCell = maxFacesPerCell
        self.XYFaceIDs = XYFaceIDs
        self.XZFaceIDs = XZFaceIDs
        self.YZFaceIDs = YZFaceIDs
        self.faceCellIDs = faceCellIDs
        self.cellFaceIDs = cellFaceIDs
        self.mesh = mesh

    def _getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        XYids = self.XYFaceIDs
        XZids = self.XZFaceIDs
        YZids = self.YZFaceIDs
        
        exteriorIDs = numerix.concatenate((numerix.ravel(XYids[...,      0].swapaxes(0,1)), 
                                           numerix.ravel(XYids[...,     -1].swapaxes(0,1)),
                                           numerix.ravel(XZids[...,  0,...]), 
                                           numerix.ravel(XZids[..., -1,...]),
                                           numerix.ravel(YZids[ 0,     ...]), 
                                           numerix.ravel(YZids[-1,     ...])))
                                                     
        from fipy.variables.faceVariable import FaceVariable
        exteriorFaces = FaceVariable(mesh=self.mesh, value=False)
        exteriorFaces[exteriorIDs] = True
        return exteriorFaces

    def _getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells
        """
        XYids = self.XYFaceIDs
        XZids = self.XZFaceIDs
        YZids = self.YZFaceIDs
        
        interiorIDs = numerix.concatenate((numerix.ravel(XYids[ ...     ,1:-1]),
                                           numerix.ravel(XZids[ ...,1:-1, ...]),
                                           numerix.ravel(YZids[1:-1,      ...].swapaxes(0,1))))
                                                     
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    def _getCellToFaceOrientations(self):
        tmp = numerix.take(self.faceCellIDs[0], self.cellFaceIDs)
        return (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _getAdjacentCellIDs(self):
        faceCellIDs = self.faceCellIDs
        mask = faceCellIDs.getMask()
        faceCellIDs = faceCellIDs.filled()
        return ((mask[0] * faceCellIDs[1] + ~mask[0] * faceCellIDs[0]),
                (mask[1] * faceCellIDs[0] + ~mask[1] * faceCellIDs[1]))

    def _getCellToCellIDs(self):
        ids = MA.zeros((6, self.nx, self.ny, self.nz), 'l')
        indices = numerix.indices((self.nx, self.ny, self.nz))
        ids[0] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx - 1
        ids[1] = indices[0] + (indices[1] + indices[2] * self.ny) * self.nx + 1
        ids[2] = indices[0] + (indices[1] + indices[2] * self.ny - self.nz) * self.nx
        ids[3] = indices[0] + (indices[1] + indices[2] * self.ny + self.nz) * self.nx
        ids[4] = indices[0] + (indices[1] + (indices[2] - 1) * self.ny) * self.nx
        ids[5] = indices[0] + (indices[1] + (indices[2] + 1) * self.ny) * self.nx
        
        ids[0, 0,    ...] = MA.masked
        ids[1,-1,    ...] = MA.masked
        ids[2,..., 0,...] = MA.masked
        ids[3,...,-1,...] = MA.masked
        ids[4,...,     0] = MA.masked
        ids[5,...,    -1] = MA.masked

        return CellVariable(mesh=self.mesh, 
                            value=MA.reshape(ids.swapaxes(1,3), (6, self.numberOfCells)))
        
    def _getCellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self.maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self.cellToCellIDs
        mask = cellToCellIDs.getMask()
        return mask * cellIDs + ~mask * cellToCellIDs.filled()

    """Properties conforming to the MeshTopology interface."""
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
                   
