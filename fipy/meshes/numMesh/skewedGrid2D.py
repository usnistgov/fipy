#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "SkewedGrid2D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 1/3/07 {3:10:37 PM} 
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

__docformat__ = 'restructuredtext'
 
from fipy.tools import numerix
from fipy.tools.numerix import random
from fipy.meshes.numMesh.mesh2D import Mesh2D
from fipy.meshes.numMesh.face import Face
import fipy.tools.vector as vector
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.meshes.meshIterator import FaceIterator

class SkewedGrid2D(Mesh2D):
    """
    Creates a 2D grid mesh with horizontal faces numbered first and then
    vertical faces.  The points are skewed by a random amount (between `rand`
    and `-rand`) in the X and Y directions.
    """
    def __init__(self, dx = 1., dy = 1., nx = None, ny = 1, rand = 0):
        self.nx = nx
        self.ny = ny
        
        self.dx = PhysicalField(value = dx)
        scale = PhysicalField(value = 1, unit = self.dx.getUnit())
        self.dx /= scale
        
        self.dy = PhysicalField(value = dy)
        if self.dy.getUnit().isDimensionless():
            self.dy = dy
        else:
            self.dy /= scale
        
        self.numberOfVertices = (self.nx + 1) * (self.ny + 1)
        
        vertices = self._createVertices()
        changedVertices = numerix.zeros(vertices.shape)
        changedVertices = changedVertices.astype(numerix.Float)
        for i in range(len(vertices)):
            if((i % (nx+1)) != 0 and (i % (nx+1)) != nx and (i / nx+1) != 0 and (i / nx+1) != ny):
                changedVertices[i, 0] = vertices[i, 0] + (rand * ((random.random() * 2) - 1))
                changedVertices[i, 1] = vertices[i, 1] + (rand * ((random.random() * 2) - 1))
            else:
                changedVertices[i, 0] = vertices[i, 0]
                changedVertices[i, 1] = vertices[i, 1]
        faces = self._createFaces()
        cells = self._createCells()
        Mesh2D.__init__(self, changedVertices, faces, cells)
        
        self.setScale(value = scale)
        
    def _createVertices(self):
        x = numerix.arange(self.nx + 1) * self.dx
        y = numerix.arange(self.ny + 1) * self.dy
        x = numerix.resize(x, (self.numberOfVertices,))
        y = numerix.repeat(y, self.nx + 1)
        return numerix.transpose(numerix.array((x, y)))
    
    def _createFaces(self):
        """
        v1, v2 refer to the cells.
        Horizontel faces are first
        """
        v1 = numerix.arange(self.numberOfVertices)
        v2 = v1 + 1
        horizontalFaces = vector.prune(numerix.transpose(numerix.array((v1, v2))), self.nx + 1, self.nx)
        v1 = numerix.arange(self.numberOfVertices - (self.nx + 1))
        v2 = v1 + self.nx + 1
        verticalFaces =  numerix.transpose(numerix.array((v1, v2)))

        ## reverse some of the face orientations to obtain the correct normals

        tmp = horizontalFaces.copy()
        horizontalFaces[:self.nx, 0] = tmp[:self.nx, 1]
        horizontalFaces[:self.nx, 1] = tmp[:self.nx, 0]

        tmp = verticalFaces.copy()
        verticalFaces[:, 0] = tmp[:, 1]
        verticalFaces[:, 1] = tmp[:, 0]
        verticalFaces[::(self.nx + 1), 0] = tmp[::(self.nx + 1), 0]
        verticalFaces[::(self.nx + 1), 1] = tmp[::(self.nx + 1), 1]

        return numerix.concatenate((horizontalFaces, verticalFaces))

    def _createCells(self):
        """
        cells = (f1, f2, f3, f4) going anticlock wise.
        f1 etx refer to the faces
        """
        self.numberOfHorizontalFaces = self.nx * (self.ny + 1)
        self.numberOfFaces = self.numberOfHorizontalFaces + self.ny * (self.nx + 1)
        f1 = numerix.arange(self.numberOfHorizontalFaces - self.nx)
        f3 = f1 + self.nx
        f2 = vector.prune(numerix.arange(self.numberOfHorizontalFaces,  self.numberOfFaces), self.nx + 1)
        f4 = f2 - 1
        return numerix.transpose(numerix.array((f1, f2, f3, f4)))

    def getFacesLeft(self):
        """Return list of faces on left boundary of Grid2D.
        """
        return FaceIterator(mesh = self, ids = numerix.arange(self.numberOfHorizontalFaces, self.numberOfFaces, self.nx + 1))
    
    def getFacesRight(self):
        """Return list of faces on right boundary of Grid2D.
        """
        return FaceIterator(mesh = self, ids = numerix.arange(self.numberOfHorizontalFaces + self.nx, self.numberOfFaces, self.nx + 1))
        
    def getFacesTop(self):
        """Return list of faces on top boundary of Grid2D.
        """
        return FaceIterator(mesh = self, ids = numerix.arange(self.numberOfHorizontalFaces - self.nx, self.numberOfHorizontalFaces))
        
    def getFacesBottom(self):
        """Return list of faces on bottom boundary of Grid2D.
        """
        return FaceIterator(mesh = self, ids = numerix.arange(self.nx))
        
    def getScale(self):
        return self.scale['length']
        
    def getPhysicalShape(self):
        """Return physical dimensions of Grid2D.
        """
        return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale()))

    def _getMeshSpacing(self):
        return numerix.array((self.dx,self.dy))
    
    def getShape(self):
        return (self.nx, self.ny)
    
## pickling

    def __getstate__(self):
        return {
            'dx' : self.dx,            
            'dy' : self.dy,
            'nx' : self.nx,
            'ny' : self.ny
        }

    def __setstate__(self, dict):
        self.__init__(dx = dict['dx'], dy = dict['dy'], nx = dict['nx'], ny = dict['ny'])

