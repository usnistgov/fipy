#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
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
 
from abstractGridBuilder import AbstractGridBuilder

from fipy.tools import inline  
from fipy.tools import numerix
from fipy.tools import vector
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.meshes.builders.utilityClasses import (UniformNumPts,
                                                 DOffsets,
                                                 UniformOrigin,
                                                 NonuniformNumPts)

class Grid2DBuilder(AbstractGridBuilder):

    def buildGridData(self, *args, **kwargs):
        # call super for side-effects
        super(Grid2DBuilder, self).buildGridData(*args, **kwargs)

        self.numberOfVerticalColumns = self.spatialDict["numVerticalCols"]
        self.numberOfHorizontalRows = self.spatialDict["numHorizontalRows"]

    @property
    def _specificGridData(self):
        return [self.numberOfHorizontalRows,
                self.numberOfVerticalColumns,
                self.numberOfHorizontalFaces]  
     
    @staticmethod
    def createVertices(nx, ny, dx, dy, numVerts, numVertCols):
        x = AbstractGridBuilder.calcVertexCoordinates(dx, nx)
        x = numerix.resize(x, (numVerts,))
            
        y = AbstractGridBuilder.calcVertexCoordinates(dy, ny)
        y = numerix.repeat(y, numVertCols)
        
        return numerix.array((x, y))
    
    @staticmethod
    def createFaces(nx, numVerts, numVertCols):
        """
        v1, v2 refer to the vertices.
        Horizontal faces are first

        Ugly return to avoid side-effects.
        """
        v1 = numerix.arange(numVerts)
        v2 = v1 + 1

        horizontalFaces = vector.prune(numerix.array((v1, v2)), numVertCols, nx, axis=1)

        v1 = numerix.arange(numVerts - numVertCols)
        v2 = v1 + numVertCols
        verticalFaces =  numerix.array((v1, v2))

        ## The cell normals must point out of the cell.
        ## The left and bottom faces have only one neighboring cell,
        ## in the 2nd neighbor position (there is nothing in the 1st).
        ## 
        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[0,:nx] = tmp[1,:nx]
        horizontalFaces[1,:nx] = tmp[0,:nx]

        numberOfHorizontalFaces = horizontalFaces.shape[-1]

        tmp = verticalFaces.copy()
        verticalFaces[0, :] = tmp[1, :]
        verticalFaces[1, :] = tmp[0, :]
        if numVertCols > 0:
            verticalFaces[0, ::numVertCols] = tmp[0, ::numVertCols]
            verticalFaces[1, ::numVertCols] = tmp[1,::numVertCols]

        return (numerix.concatenate((horizontalFaces, verticalFaces), axis=1),
                numberOfHorizontalFaces)
           
    @staticmethod
    def createCells(nx, ny, numFaces, numHorizFaces, numVertCols):
        """
        cells = (f1, f2, f3, f4) going anticlock wise.
        f1 etc. refer to the faces
        """
        return inline._optionalInline(Grid2DBuilder._createCellsIn,
                                      Grid2DBuilder._createCellsPy,
                                      nx, ny, numFaces, numHorizFaces,
                                      numVertCols)

    
    @staticmethod
    def _createCellsPy(nx, ny, numFaces, numHorizFaces, numVertCols):
        cellFaceIDs = numerix.zeros((4, nx * ny))
        faceIDs = numerix.arange(numFaces)
        if numFaces > 0:
            cellFaceIDs[0,:] = faceIDs[:numHorizFaces - nx]
            cellFaceIDs[2,:] = cellFaceIDs[0,:] + nx
            cellFaceIDs[1,:] = vector.prune(faceIDs[numHorizFaces:], 
                                            numVertCols)
            cellFaceIDs[3,:] = cellFaceIDs[1,:] - 1
        return cellFaceIDs

    @staticmethod
    def _createCellsIn(nx, ny, numFaces, numHorizFaces, numVertCols):
        cellFaceIDs = numerix.zeros((4, nx * ny))
        
        inline._runInline("""
            int ID = j * ni + i;
            int NCELLS = ni * nj;
            cellFaceIDs[ID + 0 * NCELLS] = ID;
            cellFaceIDs[ID + 2 * NCELLS] = cellFaceIDs[ID + 0 * NCELLS] + ni;
            cellFaceIDs[ID + 3 * NCELLS] = horizontalFaces + ID + j;
            cellFaceIDs[ID + 1 * NCELLS] = cellFaceIDs[ID + 3 * NCELLS] + 1;
        """,
        horizontalFaces=numHorizFaces,
        cellFaceIDs=cellFaceIDs,
        ni=nx,
        nj=ny)

        return cellFaceIDs

    def _packOverlap(self, first, second):
        return {'left': 0, 'right': 0, 'bottom': first, 'top': second}  

    def _packOffset(self, arg):
        return (0, arg)

class NonuniformGrid2DBuilder(Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = NonuniformNumPts

        super(NonuniformGrid2DBuilder, self).__init__()

    def buildGridData(self, *args, **kwargs):
        # call super for side-effects
        super(NonuniformGrid2DBuilder, self).buildGridData(*args, **kwargs)

        (self.offsets, 
         self.ds) = DOffsets.calcDOffsets(self.ds, self.ns, self.offset)

        self.vertices = Grid2DBuilder.createVertices(self.ns[0], self.ns[1],
                                        self.ds[0], self.ds[1],
                                        self.numberOfVertices, 
                                        self.numberOfVerticalColumns) \
                          + ((self.offsets[0],), (self.offsets[1],)) 

        (self.faces,
         self.numberOfHorizontalFaces) = Grid2DBuilder.createFaces(self.ns[0],
                                          self.numberOfVertices,
                                          self.numberOfVerticalColumns)
        self.numberOfFaces = len(self.faces[0])
        self.cells = Grid2DBuilder.createCells(self.ns[0], self.ns[1],
                                               self.numberOfFaces,
                                               self.numberOfHorizontalFaces,
                                               self.numberOfVerticalColumns)

    @property
    def _specificGridData(self):
        return super(NonuniformGrid2DBuilder, self)._specificGridData \
                 + [self.vertices,
                    self.faces,
                    self.cells,
                    self.offsets]
               


class UniformGrid2DBuilder(Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = UniformNumPts

        super(UniformGrid2DBuilder, self).__init__()

    def buildGridData(self, ds, ns, overlap, communicator, origin):
        # call super for side-effects
        super(UniformGrid2DBuilder, self).buildGridData(ds, ns, overlap,
                                                        communicator)
        
        self.origin = UniformOrigin.calcOrigin(origin, 
                                               self.offset, self.ds, self.scale)
           
        self.numberOfHorizontalFaces = self.ns[0] * self.numberOfHorizontalRows
        self.numberOfVerticalFaces = self.numberOfVerticalColumns * self.ns[1]
        self.numberOfFaces = self.numberOfHorizontalFaces \
                               + self.numberOfVerticalFaces

    @property
    def _specificGridData(self):
        return super(UniformGrid2DBuilder, self)._specificGridData \
                + [self.numberOfVerticalFaces,
                   self.origin]
                   

    
                                  
