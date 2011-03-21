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
from fipy.tools import inline
from fipy.tools.numerix import MA

from abstractMeshTopology import AbstractMeshTopology  

class UniformMeshTopology2D(AbstractMeshTopology):

    def __init__(self, nx, ny,
                       numberOfFaces, numberOfCells,
                       numberOfVerticalColumns,
                       numberOfHorizontalFaces,
                       numberOfHorizontalRows,
                       maxFacesPerCell,
                       mesh):
        self.nx = nx
        self.ny = ny
        self.numberOfFaces = numberOfFaces
        self.numberOfCells = numberOfCells
        self.numberOfVerticalColumns = numberOfVerticalColumns
        self.numberOfHorizontalFaces = numberOfHorizontalFaces
        self.numberOfHorizontalRows = numberOfHorizontalRows
        self.maxFacesPerCell = maxFacesPerCell
        self.mesh = mesh

    def _getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        exteriorIDs = numerix.concatenate((numerix.arange(0, self.nx),
                                           numerix.arange(0, self.nx) + self.nx * self.ny,
                                           numerix.arange(0, self.ny) * self.numberOfVerticalColumns + self.numberOfHorizontalFaces,
                                           numerix.arange(0, self.ny) * self.numberOfVerticalColumns + self.numberOfHorizontalFaces + self.nx))
                       
        from fipy.variables.faceVariable import FaceVariable
        exteriorFaces = FaceVariable(mesh=self.mesh, value=False)
        exteriorFaces[exteriorIDs] = True
        return exteriorFaces

    def _getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        Hids = numerix.arange(0, self.numberOfHorizontalFaces)
        Hids = numerix.reshape(Hids, (self.numberOfHorizontalRows, self.nx))
        Hids = Hids[1:-1,...]
        
        Vids = numerix.arange(self.numberOfHorizontalFaces, self.numberOfFaces)
        Vids = numerix.reshape(Vids, (self.ny, self.numberOfVerticalColumns))
        Vids = Vids[...,1:-1]
        
        interiorIDs = numerix.concatenate((numerix.reshape(Hids, (self.nx * (self.ny - 1),)), 
                                           numerix.reshape(Vids, ((self.nx - 1) * self.ny,))))
                                           
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    def _getCellToFaceOrientations(self):
        cellFaceOrientations = numerix.ones((4, self.numberOfCells))
        if self.numberOfCells > 0:
            cellFaceOrientations[0, self.nx:] = -1
            cellFaceOrientations[3, :] = -1
            cellFaceOrientations[3, ::self.nx] = 1
        return cellFaceOrientations

    if inline.doInline:
        def _getAdjacentCellIDs(self):
            faceCellIDs0 =  numerix.zeros(self.numberOfFaces)
            faceCellIDs1 =  numerix.zeros(self.numberOfFaces)

            inline._runInline("""
                int ID = j * ni + i;

                faceCellIDs0[ID] = ID - ni;
                faceCellIDs1[ID] = ID;

                faceCellIDs0[ID + Nhor + j] = ID - 1;
                faceCellIDs1[ID + Nhor + j] = ID;

                if (j == 0) {
                    faceCellIDs0[ID] = ID;
                }

                if (j == nj - 1) {
                    faceCellIDs0[ID + ni] = ID;
                    faceCellIDs1[ID + ni] = ID;
                }

                if (i == 0) {
                    faceCellIDs0[ID + Nhor + j] = ID;
                }

                if ( i == ni - 1 ) {
                    faceCellIDs0[ID + Nhor + j + 1] = ID;
                    faceCellIDs1[ID + Nhor + j + 1] = ID;
                }
                
            """,
            Nhor=self.numberOfHorizontalFaces,
            faceCellIDs0=faceCellIDs0,
            faceCellIDs1=faceCellIDs1,
            ni=self.nx,
            nj=self.ny)

            return (faceCellIDs0, faceCellIDs1)
    else:
        def _getAdjacentCellIDs(self):
            Hids = numerix.zeros((self.numberOfHorizontalRows, self.nx, 2))
            indices = numerix.indices((self.numberOfHorizontalRows, self.nx))
            
            Hids[...,1] = indices[1] + indices[0] * self.nx
            Hids[...,0] = Hids[...,1] - self.nx
            
            if self.numberOfHorizontalRows > 0:
                Hids[0,...,0] = Hids[0,...,1]
                Hids[0,...,1] = Hids[0,...,0]
                Hids[-1,...,1] = Hids[-1,...,0]
          
            Vids = numerix.zeros((self.ny, self.numberOfVerticalColumns, 2))
            indices = numerix.indices((self.ny, self.numberOfVerticalColumns))
            Vids[...,1] = indices[1] + indices[0] * self.nx
            Vids[...,0] = Vids[...,1] - 1
            
            if self.numberOfVerticalColumns > 0:
                Vids[...,0,0] = Vids[...,0,1]
                Vids[...,0,1] = Vids[...,0,0]
                Vids[...,-1,1] = Vids[...,-1,0]

            faceCellIDs =  numerix.concatenate((numerix.reshape(Hids, (self.numberOfHorizontalFaces, 2)), 
                                                numerix.reshape(Vids, (self.numberOfFaces - self.numberOfHorizontalFaces, 2))))

            return (faceCellIDs[:,0], faceCellIDs[:,1])

    def _getCellToCellIDs(self):
        ids = MA.zeros((4, self.nx, self.ny), 'l')
        indices = numerix.indices((self.nx, self.ny))
        ids[0] = indices[0] + (indices[1] - 1) * self.nx
        ids[1] = (indices[0] + 1) + indices[1] * self.nx
        ids[2] = indices[0] + (indices[1] + 1) * self.nx
        ids[3] = (indices[0] - 1) + indices[1] * self.nx
        
        if self.ny > 0:
            ids[0,..., 0] = MA.masked
            ids[2,...,-1] = MA.masked
        if self.nx > 0:
            ids[1,-1,...] = MA.masked
            ids[3, 0,...] = MA.masked
        
        return MA.reshape(ids.swapaxes(1,2), (4, self.numberOfCells))
        
    def _getCellToCellIDsFilled(self):
        N = self.numberOfCells
        M = self.maxFacesPerCell
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self.cellToCellIDs
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)

    """Properties conforming to the MeshTopology interface."""
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
                
