#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PFM - Python-based phase field solver
 # 
 #  FILE: "grid2D.py"
 #                                    created: 11/10/03 {3:30:42 PM} 
 #                                last update: 11/24/03 {7:29:24 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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

"""2D rectangular Mesh
"""

from mesh import Mesh
from vertex import Vertex
from face2D import Face2D
from cell import Cell
import Numeric

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
	
	Arguments:
	    
	    'dx' -- dimension of each cell in **x** direction

	    'dy' -- dimension of each cell in **y** direction

	    'nx' -- number of cells in **x** direction

	    'ny' -- number of cells in **y** direction
	"""
        self.nx=nx
        self.ny=ny
        self.dx=dx
        self.dy=dy
	vertices = self.createVertices()
	rowFaces,colFaces = self.createFaces(vertices)
	cells = self.createCells(rowFaces,colFaces)
	faces,interiorFaces = self.reorderFaces(rowFaces,colFaces)
#         interiorFaces = self.createInteriorFaces(faces)
	
	Mesh.__init__(self, cells, faces, interiorFaces, vertices)
		
    def createVertices(self):
	"""Return list of Vertices
	"""
	vertices = ()
        ny=self.ny
        nx=self.nx
        dx=self.dx
        dy=self.dy
	for j in range(ny+1):
	    for	i in range(nx+1):
		vertices += (Vertex(Numeric.array([i * dx, j * dy],'d')),)
        return vertices	
		    
    def createFaces(self, vertices):
	"""Return list of Faces bounded by vertices.
	"""
        nx=self.nx
        ny=self.ny
	faces = ()
# 	# create interior row faces
# 	for j in range(ny-1):
# 	    for i in range(nx):
# 		faces += (Face2D((vertices[i + (j + 1) * (nx + 1)],vertices[i + 1 + (j + 1) * (nx + 1)]),id),)
# 		id += 1
# 	
# 	# create interior column faces
# 	for j in range(ny):
# 	    for i in range(nx-1):
# 		faces += (Face2D((vertices[i + 1 + j * (nx + 1)],vertices[i + 1 + (j + 1) * (nx + 1)]),id),)
# 		id += 1
# 
# 	# create exterior row faces
# 	for i in range(nx):
# 	    faces += (Face2D((vertices[i],vertices[i + 1]),id),)
# 	    id += 1
# 	for i in range(nx):
# 	    faces += (Face2D((vertices[i + ny * (nx + 1)],vertices[i + 1 + ny * (nx + 1)]),id),)
# 	    id += 1
# 	
# 	# create exterior column faces
# 	for j in range(ny):
# 	    faces += (Face2D((vertices[j * (nx + 1)],vertices[(j + 1) * (nx + 1)]),id),)
# 	    id += 1
# 	for j in range(ny):
# 	    faces += (Face2D((vertices[nx + j * (nx + 1)],vertices[nx + (j + 1) * (nx + 1)]),id),)
# 	    id += 1

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
	
    def reorderFaces(self,rowFaces,colFaces):
	def flatten(T):
	    """Stolen (with small modification) from 
	    <http://mail.python.org/pipermail/tutor/2001-April/005025.html>
	    """
	    def isTuple(x): return type(x) == type(())
	    
	    if not isTuple(T): return (T,)
	    elif len(T) == 0: return ()
	    else: return flatten(T[0]) + flatten(T[1:]) 

# 	print "rowFaces: ", rowFaces[1:-1]
# 	print "colFaces: ", colFaces[:][1:-1]
	
	faces = ()
	
	faces = rowFaces[1:-1]
	for colFace in colFaces:
	    faces += colFace[1:-1]
	faces = flatten(faces)
	interiorFaces = faces
	
# 	print "interiorFaces: ", interiorFaces
	
	faces += rowFaces[0] + rowFaces[-1]
	for colFace in colFaces:
	    faces += (colFace[0],)
	for colFace in colFaces:
	    faces += (colFace[-1],)
	
# 	print "faces: ", faces
	
	id = 0
	for face in faces:
	    face.setId(id)
	    id += 1
	    
# 	print "faces: ", faces
# 	print "interiorFaces: ", interiorFaces
	
	
	return (faces, interiorFaces)
	
    def createCells(self,rowFaces,colFaces):
	"""Return list of Cells.
	"""
	nx=self.nx
	ny=self.ny
	cells = ()
	for j in range(ny):
	    for i in range(nx):
                id = j * nx + i
		cells += (
		    Cell(
		    (rowFaces[j][i],
		    rowFaces[j+1][i],
		    colFaces[j][i],
		    colFaces[j][i+1]),id
		    ),
		    )                
	return cells

#     def createCells(self, faces):
# 	"""Return list of Cells.
# 	"""
#         nx=self.nx
#         ny=self.ny
# 	cells = ()
# 	for j in range(ny):
# 	    for i in range(nx):
#                 id = j * nx + i
# 		cells += (
#                     Cell(
#                     (faces[i + j * nx],
#                     faces[i + (j+1) * nx],
#                     faces[nx * (ny + 1) + i + j * (nx + 1)],
#                     faces[nx * (ny + 1) + i + 1 + j * (nx + 1)]),id
#                     ),
#                     )                
# 	return cells

    def createInteriorFaces(self,faces):
	"""Return list of faces that are not on boundary of Grid2D.
	"""
        interiorFaces = ()
        for face in faces:
            if len(face.getCells()) == 2:
                interiorFaces += (face,)
        return interiorFaces

#     def getFacesLeft(self):
# 	"""Return list of faces on left boundary of Grid2D.
# 	"""
#         nx=self.nx
#         ny=self.ny
#         facesLeft = ()
#         for i in range(ny):
#             facesLeft += (self.faces[nx * (ny + 1) + (nx + 1) * i],)
#         return facesLeft

    def getFacesLeft(self):
	"""Return list of faces on left boundary of Grid2D.
	"""
	nx=self.nx
	ny=self.ny
	numInt = len(self.interiorFaces)
	return self.faces[numInt + 2 * nx:numInt + 2 * nx + ny]
	
#     def getFacesRight(self):
# 	"""Return list of faces on right boundary of Grid2D.
# 	"""
#         nx=self.nx
#         ny=self.ny
#         facesRight = ()
#         for i in range( ny):
#             facesRight += (self.faces[nx * (ny + 1) + (nx + 1) * (i + 1) - 1],)
#         return facesRight

    def getFacesRight(self):
	"""Return list of faces on right boundary of Grid2D.
	"""
	nx=self.nx
	ny=self.ny
	numInt = len(self.interiorFaces)
	return self.faces[numInt + 2 * nx + ny:numInt + 2 * nx + 2 * ny]
	
#     def getFacesTop(self):
# 	"""Return list of faces on top boundary of Grid2D.
# 	"""
#         nx=self.nx
#         ny=self.ny
#         return self.faces[nx*ny:nx*ny+nx]

    def getFacesTop(self):
	"""Return list of faces on top boundary of Grid2D.
	"""
	nx=self.nx
	ny=self.ny
	numInt = len(self.interiorFaces)
	return self.faces[numInt + nx:numInt + 2 * nx]
	
#     def getFacesBottom(self):
# 	"""Return list of faces on bottom boundary of Grid2D.
# 	"""
#         nx=self.nx
#         return self.faces[0:nx]
                       
    def getFacesBottom(self):
	"""Return list of faces on bottom boundary of Grid2D.
	"""
	nx=self.nx
	ny=self.ny
	numInt = len(self.interiorFaces)
	return self.faces[numInt:numInt + nx]
	    
    def getShape(self):
	"""Return cell dimensions 'Grid2D'.
	"""
        return (self.nx,self.ny)
        
    def makeGridData(self,array):
	"""Return 'array' data mapped onto cell geometry of 'Grid2D'.
	"""
        return Numeric.reshape(array,self.getShape())

    def getPhysicalShape(self):
	"""Return physical dimensions of Grid2D.
	"""
        return (self.nx*self.dx,self.ny*self.dy)
