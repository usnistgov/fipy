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

"""TODO: prettify"""

class UniformMeshTopology1D(AbstractMeshTopology):

    def __init__(self, mesh):
        self._exteriorFaces = mesh.facesLeft | mesh.facesRight
        self.mesh = mesh

    def _getInteriorFaces(self):
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[numerix.arange(self.mesh.numberOfFaces-2) + 1] = True
        return interiorFaces
            
    def _getCellToFaceOrientations(self):
        orientations = numerix.ones((2, self.mesh.numberOfCells))
        if self.mesh.numberOfCells > 0:
            orientations[0] *= -1
            orientations[0,0] = 1
        return orientations

    def _getAdjacentCellIDs(self):
        c1 = numerix.arange(self.mesh.numberOfFaces)
        ids = numerix.array((c1 - 1, c1))
        if self.mesh.numberOfFaces > 0:
            ids[0,0] = ids[1,0]
            ids[1,-1] = ids[0,-1]
        return ids[0], ids[1]

    def _getCellToCellIDs(self):
        c1 = numerix.arange(self.mesh.numberOfCells)
        ids = MA.array((c1 - 1, c1 + 1))
        if self.mesh.numberOfCells > 0:
            ids[0,0] = MA.masked
            ids[1,-1] = MA.masked
        return ids
        
    def _getCellToCellIDsFilled(self):
        ids = self._getCellToCellIDs().filled()
        if self.mesh.numberOfCells > 0:
            ids[0,0] = 0
            ids[1,-1] = self.mesh.numberOfCells - 1
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

class UniformMeshTopology2D(AbstractMeshTopology):

    def __init__(self, mesh):
        self.mesh = mesh

    def _getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        exteriorIDs = numerix.concatenate((numerix.arange(0, self.mesh.nx),
                                           numerix.arange(0, self.mesh.nx) + self.mesh.nx * self.mesh.ny,
                                           numerix.arange(0, self.mesh.ny) * self.mesh.numberOfVerticalColumns + self.mesh.numberOfHorizontalFaces,
                                           numerix.arange(0, self.mesh.ny) * self.mesh.numberOfVerticalColumns + self.mesh.numberOfHorizontalFaces + self.mesh.nx))
                       
        from fipy.variables.faceVariable import FaceVariable
        exteriorFaces = FaceVariable(mesh=self.mesh, value=False)
        exteriorFaces[exteriorIDs] = True
        return exteriorFaces

    def _getInteriorFaces(self):
        """
        Return only the faces that have two neighboring cells.
        """
        Hids = numerix.arange(0, self.mesh.numberOfHorizontalFaces)
        Hids = numerix.reshape(Hids, (self.mesh.numberOfHorizontalRows, self.mesh.nx))
        Hids = Hids[1:-1,...]
        
        Vids = numerix.arange(self.mesh.numberOfHorizontalFaces, self.mesh.numberOfFaces)
        Vids = numerix.reshape(Vids, (self.mesh.ny, self.mesh.numberOfVerticalColumns))
        Vids = Vids[...,1:-1]
        
        interiorIDs = numerix.concatenate((numerix.reshape(Hids, (self.mesh.nx * (self.mesh.ny - 1),)), 
                                           numerix.reshape(Vids, ((self.mesh.nx - 1) * self.mesh.ny,))))
                                           
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    def _getCellToFaceOrientations(self):
        cellFaceOrientations = numerix.ones((4, self.mesh.numberOfCells))
        if self.mesh.numberOfCells > 0:
            cellFaceOrientations[0, self.mesh.nx:] = -1
            cellFaceOrientations[3, :] = -1
            cellFaceOrientations[3, ::self.mesh.nx] = 1
        return cellFaceOrientations

    def _getAdjacentCellIDs(self):
        return inline._optionalInline(self._getAdjacentCellIDsIn, self._getAdjacentCellIDsPy)
    
    def _getAdjacentCellIDsIn(self):
        faceCellIDs0 =  numerix.zeros(self.mesh.numberOfFaces)
        faceCellIDs1 =  numerix.zeros(self.mesh.numberOfFaces)

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
        Nhor=self.mesh.numberOfHorizontalFaces,
        faceCellIDs0=faceCellIDs0,
        faceCellIDs1=faceCellIDs1,
        ni=self.mesh.nx,
        nj=self.mesh.ny)

        return (faceCellIDs0, faceCellIDs1)

    def _getAdjacentCellIDsPy(self):
        Hids = numerix.zeros((self.mesh.numberOfHorizontalRows, self.mesh.nx, 2))
        indices = numerix.indices((self.mesh.numberOfHorizontalRows, self.mesh.nx))
        
        Hids[...,1] = indices[1] + indices[0] * self.mesh.nx
        Hids[...,0] = Hids[...,1] - self.mesh.nx
        
        if self.mesh.numberOfHorizontalRows > 0:
            Hids[0,...,0] = Hids[0,...,1]
            Hids[0,...,1] = Hids[0,...,0]
            Hids[-1,...,1] = Hids[-1,...,0]
      
        Vids = numerix.zeros((self.mesh.ny, self.mesh.numberOfVerticalColumns, 2))
        indices = numerix.indices((self.mesh.ny, self.mesh.numberOfVerticalColumns))
        Vids[...,1] = indices[1] + indices[0] * self.mesh.nx
        Vids[...,0] = Vids[...,1] - 1
        
        if self.mesh.numberOfVerticalColumns > 0:
            Vids[...,0,0] = Vids[...,0,1]
            Vids[...,0,1] = Vids[...,0,0]
            Vids[...,-1,1] = Vids[...,-1,0]

        faceCellIDs =  numerix.concatenate((numerix.reshape(Hids, (self.mesh.numberOfHorizontalFaces, 2)), 
                                            numerix.reshape(Vids, (self.mesh.numberOfFaces - self.mesh.numberOfHorizontalFaces, 2))))

        return (faceCellIDs[:,0], faceCellIDs[:,1])


    def _getCellToCellIDs(self):
        ids = MA.zeros((4, self.mesh.nx, self.mesh.ny), 'l')
        indices = numerix.indices((self.mesh.nx, self.mesh.ny))
        ids[0] = indices[0] + (indices[1] - 1) * self.mesh.nx
        ids[1] = (indices[0] + 1) + indices[1] * self.mesh.nx
        ids[2] = indices[0] + (indices[1] + 1) * self.mesh.nx
        ids[3] = (indices[0] - 1) + indices[1] * self.mesh.nx
        
        if self.mesh.ny > 0:
            ids[0,..., 0] = MA.masked
            ids[2,...,-1] = MA.masked
        if self.mesh.nx > 0:
            ids[1,-1,...] = MA.masked
            ids[3, 0,...] = MA.masked
        
        return MA.reshape(ids.swapaxes(1,2), (4, self.mesh.numberOfCells))
        
    def _getCellToCellIDsFilled(self):
        N = self.mesh.numberOfCells
        M = self.mesh._getMaxFacesPerCell()
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self._getCellToCellIDs()
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)

    """Properties conforming to the MeshTopology interface."""
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
       
class UniformMeshTopology3D(AbstractMeshTopology):

    def __init__(self, mesh):
        self.mesh = mesh

    def _getExteriorFaces(self):
        """
        Return only the faces that have one neighboring cell.
        """
        XYids = self.mesh._getXYFaceIDs()
        XZids = self.mesh._getXZFaceIDs()
        YZids = self.mesh._getYZFaceIDs()
        
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
        XYids = self.mesh._getXYFaceIDs()
        XZids = self.mesh._getXZFaceIDs()
        YZids = self.mesh._getYZFaceIDs()
        
        interiorIDs = numerix.concatenate((numerix.ravel(XYids[ ...     ,1:-1]),
                                           numerix.ravel(XZids[ ...,1:-1, ...]),
                                           numerix.ravel(YZids[1:-1,      ...].swapaxes(0,1))))
                                                     
        from fipy.variables.faceVariable import FaceVariable
        interiorFaces = FaceVariable(mesh=self.mesh, value=False)
        interiorFaces[interiorIDs] = True
        return interiorFaces

    def _getCellToFaceOrientations(self):
        tmp = numerix.take(self.mesh.faceCellIDs[0], self.mesh.cellFaceIDs)
        return (tmp == MA.indices(tmp.shape)[-1]) * 2 - 1

    def _getAdjacentCellIDs(self):
        faceCellIDs = self.mesh.faceCellIDs
        return (MA.where(MA.getmaskarray(faceCellIDs[0]), faceCellIDs[1], faceCellIDs[0]).filled(),
                MA.where(MA.getmaskarray(faceCellIDs[1]), faceCellIDs[0], faceCellIDs[1]).filled())

    def _getCellToCellIDs(self):
        ids = MA.zeros((6, self.mesh.nx, self.mesh.ny, self.mesh.nz), 'l')
        indices = numerix.indices((self.mesh.nx, self.mesh.ny, self.mesh.nz))
        ids[0] = indices[0] + (indices[1] + indices[2] * self.mesh.ny) * self.mesh.nx - 1
        ids[1] = indices[0] + (indices[1] + indices[2] * self.mesh.ny) * self.mesh.nx + 1
        ids[2] = indices[0] + (indices[1] + indices[2] * self.mesh.ny - self.mesh.nz) * self.mesh.nx
        ids[3] = indices[0] + (indices[1] + indices[2] * self.mesh.ny + self.mesh.nz) * self.mesh.nx
        ids[4] = indices[0] + (indices[1] + (indices[2] - 1) * self.mesh.ny) * self.mesh.nx
        ids[5] = indices[0] + (indices[1] + (indices[2] + 1) * self.mesh.ny) * self.mesh.nx
        
        ids[0, 0,    ...] = MA.masked
        ids[1,-1,    ...] = MA.masked
        ids[2,..., 0,...] = MA.masked
        ids[3,...,-1,...] = MA.masked
        ids[4,...,     0] = MA.masked
        ids[5,...,    -1] = MA.masked

        return MA.reshape(ids.swapaxes(1,3), (6, self.mesh.numberOfCells))
        
    def _getCellToCellIDsFilled(self):
        N = self.mesh.numberOfCells
        M = self.mesh._getMaxFacesPerCell()
        cellIDs = numerix.repeat(numerix.arange(N)[numerix.newaxis, ...], M, axis=0)
        cellToCellIDs = self._getCellToCellIDs()
        return MA.where(MA.getmaskarray(cellToCellIDs), cellIDs, cellToCellIDs)     

    """Properties conforming to the MeshTopology interface."""
    interiorFaces = property(_getInteriorFaces)
    exteriorFaces = property(_getExteriorFaces)
    cellToFaceOrientations = property(_getCellToFaceOrientations)
    adjacentCellIDs = property(_getAdjacentCellIDs)
    cellToCellIDs = property(_getCellToCellIDs)
    cellToCellIDsFilled = property(_getCellToCellIDsFilled)
                                                                                                 
