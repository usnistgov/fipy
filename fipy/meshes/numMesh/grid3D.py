#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid2D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 12/7/04 {4:27:59 PM} 
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
3D rectangular-prism Mesh

X axis runs from left to right.
Y axis runs from bottom to top.
Z axis runs from front to back.

Numbering System:

Vertices: Numbered in the usual way. X coordinate changes most quickly, then Y, then Z.

Cells: Same numbering system as vertices.

Faces: XY faces numbered first, then XZ faces, then YZ faces. Within each subcategory, it is numbered in the usual way.
"""

import Numeric

from fipy.meshes.numMesh.mesh import Mesh
from fipy.meshes.numMesh.face import Face
import fipy.tools.vector as vector
from fipy.tools.dimensions.physicalField import PhysicalField


class Grid3D(Mesh):
    """
    Creates a 3D grid mesh.
    """
    def __init__(self, dx = 1., dy = 1., dz = 1., nx = 1, ny = 1, nz = 1):
        self.nx = nx
        self.ny = ny
	self.nz = nz
        
	self.dx = PhysicalField(value = dx)
	scale = PhysicalField(value = 1, unit = self.dx.getUnit())
	self.dx /= scale
	
	self.dy = PhysicalField(value = dy)
	if self.dy.getUnit().isDimensionless():
	    self.dy = dy
	else:
	    self.dy /= scale

	self.dz = PhysicalField(value = dz)
	if self.dz.getUnit().isDimensionless():
	    self.dz = dz
	else:
	    self.dz /= scale
	
        self.numberOfVertices = (self.nx + 1) * (self.ny + 1) * (self.nz + 1)
        
	vertices = self.createVertices()
        faces = self.createFaces()
        cells = self.createCells()
        Mesh.__init__(self, vertices, faces, cells)
	
	self.setScale(value = scale)
        
    def createVertices(self):
        x = Numeric.arange(self.nx + 1) * self.dx
        y = Numeric.arange(self.ny + 1) * self.dy
        x = Numeric.resize(x, (self.numberOfVertices,))
        y = Numeric.repeat(y, self.nx + 1)
        y = Numeric.resize(y, (self.numberOfVertices,))
        z = Numeric.arange(self.nz + 1) * self.dz
        z = Numeric.repeat(z, (self.nx + 1) * (self.ny + 1))
        return Numeric.transpose(Numeric.array((x, y, z)))
    
    def createFaces(self):
        """
        XY faces are first, then XZ faces, then YZ faces
        """
        ## do the XY faces
        v1 = Numeric.arange((self.nx + 1) * (self.ny))
        v1 = vector.prune(v1, self.nx + 1, self.nx)
        v1 = self.repeatWithOffset(v1, (self.nx + 1) * (self.ny + 1), self.nz + 1) 
        v2 = v1 + 1
        v3 = v1 + (self.nx + 2)
        v4 = v1 + (self.nx + 1)
        XYFaces = Numeric.transpose(Numeric.array((v1, v2, v3, v4)))

        ## do the XZ faces
        v1 = Numeric.arange((self.nx + 1) * (self.ny + 1))
        v1 = vector.prune(v1, self.nx + 1, self.nx)
        v1 = self.repeatWithOffset(v1, (self.nx + 1) * (self.ny + 1), self.nz)
        v2 = v1 + 1
        v3 = v1 + ((self.nx + 1)*(self.ny + 1)) + 1
        v4 = v1 + ((self.nx + 1)*(self.ny + 1))
        XZFaces = Numeric.transpose(Numeric.array((v1, v2, v3, v4)))
        
        ## do the YZ faces
        v1 = Numeric.arange((self.nx + 1) * self.ny)
        v1 = self.repeatWithOffset(v1, (self.nx + 1) * (self.ny + 1), self.nz)
        v2 = v1 + (self.nx + 1)
        v3 = v1 + ((self.nx + 1)*(self.ny + 1)) + (self.nx + 1)                                  
        v4 = v1 + ((self.nx + 1)*(self.ny + 1))
        YZFaces = Numeric.transpose(Numeric.array((v1, v2, v3, v4)))                           

        ## reverse some of the face orientations to obtain the correct normals
        ##tmp = horizontalFaces.copy()
        ##horizontalFaces[:self.nx, 0] = tmp[:self.nx, 1]
        ##horizontalFaces[:self.nx, 1] = tmp[:self.nx, 0]
        ##tmp = verticalFaces.copy()
        ##verticalFaces[:, 0] = tmp[:, 1]
        ##verticalFaces[:, 1] = tmp[:, 0]
        ##verticalFaces[::(self.nx + 1), 0] = tmp[::(self.nx + 1), 0]
        ##verticalFaces[::(self.nx + 1), 1] = tmp[::(self.nx + 1), 1]

        self.numberOfXYFaces = (self.nx * self.ny * (self.nz + 1))
        self.numberOfXZFaces = (self.nx * (self.ny + 1) * self.nz)
        self.numberOfYZFaces = ((self.nx + 1) * self.ny * self.nz)
        self.totalNumberOfFaces = self.numberOfXYFaces + self.numberOfXZFaces + self.numberOfYZFaces
        
        return Numeric.concatenate((XYFaces, XZFaces, YZFaces))
    
    def createCells(self):
        """
        cells = (front face, back face, left face, right face, bottom face, top face)
        front and back faces are YZ faces
        left and right faces are XZ faces
        top and bottom faces are XY faces
        """
        self.numberOfCells = self.nx * self.ny * self.nz
        
        ## front and back faces
        frontFaces = Numeric.arange(self.numberOfYZFaces)
        frontFaces = vector.prune(frontFaces, self.nx + 1, self.nx)
        frontFaces = frontFaces + self.numberOfXYFaces + self.numberOfXZFaces
        backFaces = frontFaces + 1

        ## left and right faces
        leftFaces = Numeric.arange(self.nx * self.ny)
        leftFaces = self.repeatWithOffset(leftFaces, self.nx * (self.ny + 1), self.nz) 
        leftFaces = Numeric.ravel(leftFaces)
        leftFaces = leftFaces + self.numberOfXYFaces
        rightFaces = leftFaces + self.nx

        ## bottom and top faces
        bottomFaces = Numeric.arange(self.nx * self.ny * self.nz)
        topFaces = bottomFaces + (self.nx * self.ny)
        
        return Numeric.transpose(Numeric.array((frontFaces, backFaces, leftFaces, rightFaces, bottomFaces, topFaces)))

    def getFacesBottom(self):
	"""Return list of faces on bottom boundary of Grid3D.
	"""
	return [Face(self, id) for id in self.repeatWithOffset(Numeric.arange(self.numberOfXYFaces, self.numberOfXYFaces + self.nx), self.nx * (self.ny + 1), self.nz)]
	
    def getFacesTop(self):
	"""Return list of faces on top boundary of Grid3D.
	"""
	return [Face(self, id) for id in self.repeatWithOffset(Numeric.arange(self.numberOfXYFaces + (self.nx * self.ny), self.numberOfXYFaces + (self.nx * self.ny) + self.nx), self.nx * (self.ny + 1), self.nz)]
	
    def getFacesBack(self):
	"""Return list of faces on back boundary of Grid3D.
	"""
	return [Face(self, id) for id in Numeric.arange(self.numberOfXYFaces - (self.nx * self.ny), self.numberOfXYFaces)]
	
    def getFacesFront(self):
	"""Return list of faces on front boundary of Grid3D.
	"""
	return [Face(self, id) for id in Numeric.arange(self.nx * self.ny)]

    def getFacesLeft(self):
        """Return list of faces on left boundary of Grid3D.
        """
        return [Face(self, id) for id in Numeric.arange(self.numberOfXYFaces + self.numberOfXZFaces, self.totalNumberOfFaces, self.nx + 1)]

    def getFacesRight(self):
        """Return list of faces on right boundary of Grid3D.
        """
        return [Face(self, id) for id in Numeric.arange(self.numberOfXYFaces + self.numberOfXZFaces + self.nx, self.totalNumberOfFaces, self.nx + 1)]
    
    def getScale(self):
	return self.scale['length']
	
    def getPhysicalShape(self):
	"""Return physical dimensions of Grid2D.
	"""
	return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale(), self.nz * self.dz * self.getScale()))

    def getMeshSpacing(self):
	return Numeric.array((self.dx, self.dy, self.dz))
    
    def getShape(self):
        return (self.nx, self.ny, self.nz)

    def repeatWithOffset(self, array, offset, reps):
        a = Numeric.fromfunction(lambda rnum, x: array + (offset * rnum), (reps, Numeric.size(array)))
        return Numeric.ravel(a)

    def calcFaceAreas(self):
        XYFaceAreas = Numeric.ones(self.numberOfXYFaces)
        XYFaceAreas = XYFaceAreas * self.dx * self.dy
        XZFaceAreas = Numeric.ones(self.numberOfXZFaces)
        XZFaceAreas = XZFaceAreas * self.dx * self.dz        
        YZFaceAreas = Numeric.ones(self.numberOfYZFaces)
        YZFaceAreas = YZFaceAreas * self.dy * self.dz
        self.faceAreas =  Numeric.concatenate((XYFaceAreas, XZFaceAreas, YZFaceAreas))

    def calcFaceNormals(self):
        XYFaceNormals = Numeric.zeros((self.numberOfXYFaces, 3))
        XYFaceNormals[(self.nx * self.ny):, 2] = 1
        XYFaceNormals[:(self.nx * self.ny), 2] = -1
        XZFaceNormals = Numeric.zeros((self.numberOfXZFaces, 3))
        xzd = Numeric.arange(self.numberOfXZFaces)
        xzd = xzd % (self.nx * (self.ny + 1))
        xzd = (xzd < self.nx)
        xzd = 1 - (2 * xzd)
        XZFaceNormals[:, 1] = xzd
        YZFaceNormals = Numeric.zeros((self.numberOfYZFaces, 3))
        YZFaceNormals[:, 0] = 1
        YZFaceNormals[::self.nx + 1, 0] = -1
        self.faceNormals = Numeric.concatenate((XYFaceNormals, XZFaceNormals, YZFaceNormals))
        
    def calcFaceTangents(self):
        ## need to see whether order matters.
        faceTangents1 = Numeric.zeros((self.totalNumberOfFaces, 3))
        faceTangents1 = faceTangents1.astype(Numeric.Float)
        faceTangents2 = Numeric.zeros((self.totalNumberOfFaces, 3))
        faceTangents2 = faceTangents2.astype(Numeric.Float)
        ## XY faces
        faceTangents1[:self.numberOfXYFaces, 0] = 1.
        faceTangents2[:self.numberOfXYFaces, 1] = 1.
        ## XZ faces
        faceTangents1[self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces, 0] = 1.
        faceTangents2[self.numberOfXYFaces:self.numberOfXYFaces + self.numberOfXZFaces, 2] = 1.
        ## YZ faces
        faceTangents1[self.numberOfXYFaces + self.numberOfXZFaces:, 1] = 1.
        faceTangents2[self.numberOfXYFaces + self.numberOfXZFaces:, 2] = 1.
        self.faceTangents1 = faceTangents1
        self.faceTangents2 = faceTangents2

    def calcHigherOrderScalings(self):
        self.scale['area'] = self.scale['length']**2
	self.scale['volume'] = self.scale['length']**3
        
## pickling

    def __getstate__(self):
        dict = {
            'dx' : self.dx,            
            'dy' : self.dy,
            'dz' : self.dz,
            'nx' : self.nx,
            'ny' : self.ny,
            'nz' : self.nz
            }
        return dict

    def __setstate__(self, dict):
        self.__init__(dx = dict['dx'], dy = dict['dy'], dz = dict['dz'], nx = dict['nx'], ny = dict['ny'], nz = dict['nz'])

