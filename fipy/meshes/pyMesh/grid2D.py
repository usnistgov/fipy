#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "grid2D.py"
 #
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
 # ###################################################################
 ##

"""2D rectangular Mesh
"""
__docformat__ = 'restructuredtext'

from fipy.tools import numerix

from fipy.meshes.pyMesh.mesh import Mesh
from fipy.meshes.pyMesh.vertex import Vertex
from fipy.meshes.pyMesh.face2D import Face2D
from fipy.meshes.pyMesh.cell import Cell
from fipy.tools.dimensions.physicalField import PhysicalField

class Grid2D(Mesh):
    """2D rectangular Mesh
    
    Numbering system

	nx=5
	
	ny=3

	Cells::
	
	    *************************************
	    *      *      *       *      *      *
	    * 10   * 11   * 12    * 13   * 14   *
	    *************************************
	    *      *      *       *      *      *
	    * 5    * 6    * 7     * 8    * 9    *
	    *************************************
	    *      *      *       *      *      *
	    * 0    * 1    * 2     * 3    * 4    *
	    *************************************

	Faces (before reordering)::
	    
	    ***15******16*****17******18****19***
	    *      *      *       *      *      *
	    32    33      34      35     36     37
	    ***10******11*****12******13*****14**
	    *      *      *       *      *      *
	    26     27     28      29     30     31
	    ***5*******6******7*******8******9***
	    *      *      *       *      *      *
	    20     21     22      23     24     25
	    ***0*******1******2*******3******4***

    
	Faces (after reordering)::
	    
	    ***27******28*****29******30****31***
	    *      *      *       *      *      *
	    34     18     19      20     21     37
	    ***5*******6******7*******8******9***
	    *      *      *       *      *      *
	    33     14     15      16     17     36
	    ***0*******1******2*******3******4***
	    *      *      *       *      *      *
	    32     10     11      12     13     35
	    ***22******23*****24******25*****26**
		
	Vertices::

	    18*****19*****20******21*****22****23
	    *      *      *       *      *      *
	    *      *      *       *      *      *
	    12*****13*****14******15*****16****17
	    *      *      *       *      *      *
	    *      *      *       *      *      *
	    6******7******8*******9******10****11
	    *      *      *       *      *      *
	    *      *      *       *      *      *
	    0******1******2*******3******4******5
    """
    
    def __init__(self, dx, dy, nx, ny):
	"""Grid2D is initialized by caller
	
	:Parameters:
	    
	  - `dx`: dimension of each cell in **x** direction
	  - `dy`: dimension of each cell in **y** direction
	  - `nx`: number of cells in **x** direction
	  - `ny`: number of cells in **y** direction
	"""
        self.nx=nx
        self.ny=ny
        self.dx=PhysicalField(value = dx)
        self.dy=PhysicalField(value = dy)
	self.scale = PhysicalField(value = 1, unit = self.dx.getUnit())
	self.dx /= self.scale
	self.dy /= self.scale
	
	vertices = self._createVertices()
	rowFaces,colFaces = self._createFaces(vertices)
	cells = self._createCells(rowFaces,colFaces)
	faces,interiorFaces = self._reorderFaces(rowFaces,colFaces)
	Mesh.__init__(self, cells, faces, interiorFaces, vertices)
		
    def _createVertices(self):
	"""Return list of `Vertex` objects
	"""
	vertices = ()
        ny=self.ny
        nx=self.nx
        dx=self.dx
        dy=self.dy
	for j in range(ny+1):
	    for	i in range(nx+1):
		vertices += (Vertex(numerix.array([i * dx, j * dy],'d')),)
## 		vertices += (Vertex(PhysicalField(value = [i * dx, j * dy])),)
        return vertices	
		    
    def _createFaces(self, vertices):
	"""Return 2-`tuple` of `Face` objects bounded by `vertices`. 
	
	First `tuple` are the `Face` objects that separate rows of `Cell` objects.  
	Second `tuple` are the `Face` objects that separate columns of `Cell 
	objects.  These initial lists are layed out for efficiency of composing 
	and indexing into the lists to compose `Cell` objects.  They will 
	subsequently be reordered for efficiency of computations.
	"""
        nx=self.nx
        ny=self.ny

	id = 0
	rowFaces = ()
	for j in range(ny+1):
	    oneRow = ()
	    for i in range(nx):
		oneRow += (Face2D((vertices[i + j * (nx + 1)],vertices[i + 1 + j * (nx + 1)]),id),)
		id += 1
	    rowFaces += (oneRow,)
	colFaces = []
	for j in range(ny):
	    oneCol = ()
	    for i in range(nx+1):
		oneCol += (Face2D((vertices[i + j * (nx + 1)],vertices[i + (j + 1) * (nx + 1)]),id),)
		id += 1
	    colFaces += (oneCol,)
	return (rowFaces,colFaces)
	
    def _reorderFaces(self,rowFaces,colFaces):
	"""Return a `tuple` of `Face` objects ordered for best efficiency.
	
	Composed from `rowFaces` and `colFaces` such that all interior faces
	are listed contiguously, rows then columns, followed by all boundary
	faces, rows then columns.
	"""
	interiorFaces = ()

	for rowFace in rowFaces[1:-1]:
	    interiorFaces += rowFace
	for colFace in colFaces:
	    interiorFaces += colFace[1:-1]
	    
	faces = interiorFaces
	faces += rowFaces[0] + rowFaces[-1]
	for colFace in colFaces:
	    faces += (colFace[0],)
	for colFace in colFaces:
	    faces += (colFace[-1],)
	
	id = 0
	for face in faces:
	    face._setID(id)
	    id += 1

	return (faces, interiorFaces)
	
    def _createCells(self,rowFaces,colFaces):
	"""Return list of `Cell` objects.
	"""
	nx=self.nx
	ny=self.ny
	cells = ()
	for j in range(ny):
	    for i in range(nx):
                id = j * nx + i
		cells += (
		    Cell(
			faces = (rowFaces[j][i],
				rowFaces[j+1][i],
				colFaces[j][i],
				colFaces[j][i+1]),
			faceOrientations = (-1,1,1,-1),
			id = id
			),
		    ) 
		    
	return cells

    def _createInteriorFaces(self,faces):
	"""Return list of faces that are not on boundary of Grid2D.
	"""
        interiorFaces = ()
        for face in faces:
            if len(face.getCells()) == 2:
                interiorFaces += (face,)
        return interiorFaces

    def getFacesLeft(self):
        """
        Return list of faces on left boundary of Grid2D with the
        x-axis running from left to right.
        """
	nx=self.nx
	ny=self.ny
	start = len(self.interiorFaces) + 2 * nx
	return self.faces[start:start + ny]
	
    def getFacesRight(self):
        """
        Return list of faces on right boundary of Grid2D with the
        x-axis running from left to right.
        """
	nx=self.nx
	ny=self.ny
	start = len(self.interiorFaces) + 2 * nx + ny
	return self.faces[start:start + ny]
	
    def getFacesTop(self):
        """
        Return list of faces on top boundary of Grid2D with the
        y-axis running from bottom to top.
        """
	nx=self.nx
	start = len(self.interiorFaces) + nx
	return self.faces[start:start + nx]
	
    def getFacesBottom(self):
        """
        Return list of faces on bottom boundary of Grid2D with the
        y-axis running from bottom to top.
        """
	nx=self.nx
	start = len(self.interiorFaces)
	return self.faces[start:start + nx]
	    
    def getShape(self):
	"""Return cell dimensions `Grid2D`.
	"""
        return (self.nx,self.ny)
        
    def getPhysicalShape(self):
	"""Return physical dimensions of Grid2D.
	"""
	return PhysicalField(value = (self.nx * self.dx * self.getScale(), self.ny * self.dy * self.getScale()))

    def _getMaxFacesPerCell(self):
        return 4

    def _getFaceAreas(self):
	return Mesh._getFaceAreas(self) * self.getScale()

    def getCellVolumes(self):
        if self.getScale() is 1:
            return Mesh.getCellVolumes(self)
        else:
            return Mesh.getCellVolumes(self) * self.getScale() * self.getScale()

    def getCellCenters(self):
        if self.getScale() is 1:
            return Mesh.getCellCenters(self)
        else:
            return Mesh.getCellCenters(self) * self.getScale()

    def _getCellDistances(self):
        if self.getScale() is 1:
            return Mesh._getCellDistances(self)
        else:
            return Mesh._getCellDistances(self) * self.getScale()
        
    def _getFaceToCellDistances(self):
        if self.getScale() is 1:
            return Mesh._getFaceToCellDistances(self)
        else:
            return Mesh._getFaceToCellDistances(self) * self.getScale()

    def _getMeshSpacing(self):
	return PhysicalField(value = ((self.dx * self.getScale(),),(self.dy * self.getScale(),)))
