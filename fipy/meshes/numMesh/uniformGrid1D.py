#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "uniformGrid1D.py"
 #                                    created: 2/22/06 {11:32:04 AM}
 #                                last update: 11/6/07 {11:16:35 AM} 
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
 #  
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

"""
1D Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools.numerix import MA

from fipy.meshes.numMesh.grid1D import Grid1D
from fipy.meshes.meshIterator import FaceIterator
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import numerix

from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable
from fipy.variables.vertexVariable import _VertexVariable

class UniformGrid1D(Grid1D):
    """
    Creates a 1D grid mesh.
    
        >>> mesh = UniformGrid1D(nx = 3)
        >>> print mesh.getCellCenters()
        [[ 0.5  1.5  2.5]]
         
    """
    def __init__(self, dx = 1., nx = 1, origin = (0,)):
        self.dim = 1
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.getUnit())
        self.dx /= scale
        
        self.origin = PhysicalField(value = origin)
        self.origin /= scale
        
        self.nx = int(nx)
        
        self.numberOfVertices = self.nx + 1
        self.numberOfFaces = self.nx + 1
        self.numberOfCells = self.nx
        
        self.exteriorFaces = self.getFacesLeft() + self.getFacesRight()
        
        self.scale = {
            'length': 1.,
            'area': 1.,
            'volume': 1.
        }
        
        self.setScale(value = scale)
        
    def _translate(self, vector):
        return UniformGrid1D(dx = self.dx, nx = self.nx, origin = self.origin + vector)

    def __mul__(self, factor):
        return UniformGrid1D(dx = self.dx * factor, nx = self.nx, origin = self.origin * factor)

    def _getConcatenableMesh(self):
        from fipy.meshes.numMesh.mesh1D import Mesh1D
        return Mesh1D(vertexCoords = self.getVertexCoords(), 
                      faceVertexIDs = self._createFaces(), 
                      cellFaceIDs = self._createCells())
                      
    def _concatenate(self, other, smallNumber):
        """
        Following test was added due to a bug in adding Meshes.
        
            >>> a = UniformGrid1D(nx=10) + 10
            >>> print a.getCellCenters()[0,0]
            10.5
            >>> b = 10 + UniformGrid1D(nx=10)
            >>> print b.getCellCenters()[0,0]
            10.5
            >>> c =  UniformGrid1D(nx=10) + (UniformGrid1D(nx=10) + 10)
            >>> print c.getCellCenters()[0,-1]
            19.5
            
        """
        return self._getConcatenableMesh()._concatenate(other = other, smallNumber = smallNumber)
        
##     get topology methods

##         from common/mesh
        
    def _getCellFaceIDs(self):
        return CellVariable(mesh=self, 
                            value=MA.array(self._createCells()), 
                            elementshape=(2,))
     
    def getInteriorFaces(self):
        return FaceIterator(mesh=self, 
                            ids=numerix.arange(self.numberOfFaces-2) + 1)
            
    def _getCellFaceOrientations(self):
        orientations = CellVariable(mesh=self, value=1., elementshape=(2,))
        orientations[0] *= -1
        orientations[0,0] = 1
        return orientations
        
    def _getAdjacentCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = numerix.array((c1 - 1, c1))
        ids[0,0] = ids[1,0]
        ids[1,-1] = ids[0,-1]
        return (FaceVariable(mesh=self, value=ids[0]), 
                FaceVariable(mesh=self, value=ids[1]))

    def _getCellToCellIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        ids = MA.array((c1 - 1, c1 + 1))
        ids[0,0] = MA.masked
        ids[1,-1] = MA.masked
        return CellVariable(mesh=self, value=ids)

    def _getCellToCellIDsFilled(self):
        ids = self._getCellToCellIDs().filled().copy()
        ids[0,0] = 0
        ids[1,-1] = self.numberOfCells - 1
        return ids
        
    def _getMaxFacesPerCell(self):
        return 2
        
##         from numMesh/mesh

    def getVertexCoords(self):
        return self.getFaceCenters()

    def getFaceCellIDs(self):
        c1 = numerix.arange(self.numberOfFaces)
        ids = MA.array((c1 - 1, c1))
        ids[0,0] = ids[1,0]
        ids[1,0] = MA.masked
        ids[1,-1] = MA.masked
        return FaceVariable(mesh=self, value=ids, elementshape=(2,))

##     get geometry methods
        
##         from common/mesh
        
    def _getFaceAreas(self):
        return FaceVariable(mesh=self, value=1.)

    def _getFaceNormals(self):
        faceNormals = FaceVariable(mesh=self, value=1., rank=1)
        # The left-most face has neighboring cells None and the left-most cell.
        # We must reverse the normal to make fluxes work correctly.
        faceNormals[...,0] *= -1
        return faceNormals

    def getCellVolumes(self):
        return CellVariable(mesh=self, value=self.dx)

    def getCellCenters(self):
        return CellVariable(mesh=self,
                            value=((numerix.arange(self.numberOfCells)[numerix.NewAxis, ...] + 0.5) 
                                   * self.dx + self.origin) * self.scale['length'],
                            rank=1)

    def _getCellDistances(self):
        distances = FaceVariable(mesh=self, value=0.)
        distances[1:-1] = self.dx
        distances[0] = self.dx / 2.
        distances[-1] = self.dx / 2.
        return distances

    def _getFaceToCellDistanceRatio(self):
        distances = FaceVariable(mesh=self, value=0.5)
        distances[0] = 1
        distances[-1] = 1
        return distances

    def _getOrientedAreaProjections(self):
        return self._getAreaProjections()

    def _getAreaProjections(self):
        return self._getFaceNormals()

    def _getOrientedFaceNormals(self):
        return self._getFaceNormals()

    def _getFaceTangents1(self):
        return FaceVariable(mesh=self, value=0., rank=1)

    def _getFaceTangents2(self):
        return FaceVariable(mesh=self, value=0., rank=1)
     
    def _getFaceAspectRatios(self):
        return 1. / self._getCellDistances()
    
    def _getCellToCellDistances(self):
        distances = MA.zeros((2, self.numberOfCells), 'd')
        distances[:] = self.dx
        distances[0,0] = self.dx / 2.
        distances[1,-1] = self.dx / 2.
        return CellVariable(mesh=self, value=distances, elementshape=(2,))

    def _getCellNormals(self):
        normals = CellVariable(mesh=self, value=1., elementshape=(1,2))
        normals[:,0] = -1
        return normals

    def _getCellAreas(self):
        return CellVariable(mesh=self, value=1, elementshape=(2,))

    def _getCellAreaProjections(self):
        return CellVariable(mesh=self,
                            value=MA.array(self._getCellNormals().getValue()),
                            elementshape=(1,2))

##         from numMesh/mesh

    def getFaceCenters(self):
        return FaceVariable(mesh=self,
                            value=numerix.arange(self.numberOfFaces)[numerix.NewAxis, ...] * self.dx + self.origin,
                            rank=1)

    def _getCellVertexIDs(self):
        c1 = numerix.arange(self.numberOfCells)
        return CellVariable(mesh=self, value=(c1 + 1, c1), elementshape=(2,))


##     scaling
    
    def _calcScaledGeometry(self):
        pass

    def _test(self):
        """
        These tests are not useful as documentation, but are here to ensure
        everything works as expected. The following was broken, now fixed.

            >>> from fipy import *
            >>> mesh = Grid1D(nx=3., dx=1.)
            >>> var = CellVariable(mesh=mesh)
            >>> DiffusionTerm().solve(var)

        """

def _test():
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
